# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] St.Jude Children's Research Hospital, Memphis, TN
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# ******************************************************************************

import os
import threading
import time

from emtools.utils import Timer, Pretty, Process
from emtools.jobs import Pipeline
from emtools.pwx import BatchManager

from pyworkflow import SCIPION_DEBUG_NOCLEAN, BETA
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.protocol import STEPS_SERIAL
from pwem.objects import Float, SetOfMovies, SetOfMicrographs
from . import ProtMotionCorr
from .protocol_base import ProtMotionCorrBase

from .. import Plugin



class ProtMotionCorrTasks(ProtMotionCorr):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'tasks'
    _devStatus = BETA
    stepsExecutionMode = STEPS_SERIAL

    def __init__(self, **kwargs):
        ProtMotionCorr.__init__(self, **kwargs)
        # Disable parallelization options just take into account GPUs
        self.numberOfMpi.set(0)
        self.numberOfThreads.set(0)
        self.allowMpi = False
        self.allowThreads = False

    # We are not using the steps mechanism for parallelism from Scipion
    def _stepsCheck(self):
        pass

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        ProtMotionCorr._defineAlignmentParams(self, form)
        self._defineStreamingParams(form)
        # Make default 1 minute for sleeping when no new input movies
        form.getParam('streamingSleepOnWait').setDefault(60)

    # --------------------------- STEPS functions -----------------------------
    def _convertInputStep(self):
        # ProtMotionCorr keeps extra/DONE only for the classic
        # ProtAlignMovies workflow. Tasks uses logical persisted state
        # and must not recreate that filesystem checkpoint.
        ProtMotionCorrBase._convertInputStep(self)

    def _insertAllSteps(self):
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        self._insertFunctionStep(self._convertInputStep, needsGPU=False)
        # Make the following step to always run, despite finished
        # this may be useful when new input items (from streaming)
        # and need to continue
        self._insertFunctionStep(self._processAllMoviesStep, Pretty.now(), needsGPU=True)

    def _processAllMoviesStep(self, when):
        self.lock = threading.Lock()
        self.processed = 0
        self.registered = 0
        self._batchFailures = []
        self._failedBatchIds = set()

        self.error(f">>> {Pretty.now()}: ----------------- "
                   f"Start processing movies----------- ")
        self.program = Plugin.getProgram()
        self.command = self._getCmd()
        self._firstTimeOutput = True

        waitSecs = self.streamingSleepOnWait.get()
        moviesIter = self._iterLogicalInputMovies(waitSecs)
        batchMgr = BatchManager(self.streamingBatchSize.get(), moviesIter,
                                self._getTmpPath())

        mc = Pipeline()
        g = mc.addGenerator(batchMgr.generate)
        gpus = self.getGpuList()
        outputQueue = None
        self.debug(f"GPUS: {gpus}")
        for gpu in gpus:
            p = mc.addProcessor(g.outputQueue, self._getMcProcessor(gpu),
                                outputQueue=outputQueue)
            outputQueue = p.outputQueue

        o1 = mc.addProcessor(outputQueue, self._getOutputProcessor())
        mc.run()

        # Worker failures are recorded inside Pipeline threads so emtools
        # can close all queues normally. Re-raise from the main Scipion
        # protocol thread so the step becomes FAILED instead of RUNNING.
        if self._batchFailures:
            raise self._batchFailures[0]

        # Mark the output as closed
        self._firstTimeOutput = False
        self._updateOutputSets([], pwobj.Set.STREAM_CLOSED)

    @staticmethod
    def _getOutputItemIds(outputSet):
        getUniqueValues = getattr(
            outputSet,
            'getUniqueValues',
            None,
        )
        if callable(getUniqueValues):
            try:
                return set(
                    getUniqueValues('id')
                )
            except Exception:
                pass

        iterItems = getattr(
            outputSet,
            'iterItems',
            None,
        )
        if callable(iterItems):
            return {
                item.getObjId()
                for item in iterItems()
            }

        return {
            item.getObjId()
            for item in outputSet
        }

    def _getRequiredOutputNamesForResume(self):
        outputNames = []

        if self._createOutputMovies():
            outputNames.append(
                'outputMovies'
            )

        if self._createOutputMicrographs():
            outputNames.append(
                'outputMicrographs'
            )

        if self._createOutputWeightedMicrographs():
            outputNames.append(
                'outputMicrographsDoseWeighted'
            )

        if self._doSplitEvenOdd():
            outputNames.extend([
                'outputMicrographsEven',
                'outputMicrographsOdd',
            ])

        return outputNames

    def _getPersistedOutputMovieIds(self):
        completedIds = None

        for outputName in self._getRequiredOutputNamesForResume():
            outputSet = getattr(
                self,
                outputName,
                None,
            )

            # A configured output that has not been published yet means
            # that no input Movie can be considered durably complete.
            if outputSet is None:
                return set()

            loadAllProperties = getattr(
                outputSet,
                'loadAllProperties',
                None,
            )
            if callable(loadAllProperties):
                loadAllProperties()

            outputIds = self._getOutputItemIds(
                outputSet
            )

            if completedIds is None:
                completedIds = outputIds
            else:
                completedIds.intersection_update(
                    outputIds
                )

        return (
            completedIds
            if completedIds is not None
            else set()
        )

    def _iterLogicalInputMovies(self, waitSecs):
        inputMovies = self.getInputMovies()
        scheduledIds = set(
            self._getPersistedOutputMovieIds()
        )

        while True:
            newMovies = []

            for movie in inputMovies.iterItems():
                objId = movie.getObjId()

                if objId in scheduledIds:
                    continue

                scheduledIds.add(objId)

                clone = getattr(
                    movie,
                    'clone',
                    None,
                )
                newMovies.append(
                    clone()
                    if callable(clone)
                    else movie
                )

            for movie in newMovies:
                yield movie

            if not inputMovies.isStreamOpen():
                break

            time.sleep(waitSecs)

            loadAllProperties = getattr(
                inputMovies,
                'loadAllProperties',
                None,
            )
            if callable(loadAllProperties):
                loadAllProperties()

    def _getLogicalOutputSet(self, outputName, setClass,
                             suffix='', fixSampling=True):
        outputSet = getattr(self, outputName, None)
        if outputSet is not None:
            loadAllProperties = getattr(
                outputSet,
                'loadAllProperties',
                None,
            )
            if callable(loadAllProperties):
                loadAllProperties()

            outputSet.enableAppend()
            return outputSet, False

        inputMovies = self.getInputMovies()
        outputSet = setClass.create(
            self._getPath(),
            template=outputName,
            suffix=suffix or None,
        )
        outputSet.copyInfo(inputMovies)

        if fixSampling:
            outputSet.setSamplingRate(
                inputMovies.getSamplingRate() * self._getBinFactor()
            )

        outputSet.setStreamState(pwobj.Set.STREAM_OPEN)
        outputSet.write()

        self._defineOutputs(**{outputName: outputSet})

        # _defineOutputs() may replace the provisional Set with the
        # canonical runtime representation (for example the PostgreSQL
        # Set exposed by ScipionWeb). Continue from that canonical object
        # so subsequent append/update/write operations target the persisted
        # output rather than the temporary Set used during creation.
        canonicalOutputSet = getattr(
            self,
            outputName,
            outputSet,
        )
        enableAppend = getattr(
            canonicalOutputSet,
            'enableAppend',
            None,
        )
        if callable(enableAppend):
            enableAppend()

        self._defineTransformRelation(
            self.inputMovies,
            canonicalOutputSet,
        )
        return canonicalOutputSet, True

    def _persistLogicalOutputSet(self, outputSet, streamMode):
        outputSet.setStreamState(streamMode)
        outputSet.write()
        self._store(outputSet)

    def _updateLogicalOutputMovieSet(self, newDone, streamMode):
        outputName = 'outputMovies'
        saveMovie = self.getAttributeValue('doSaveMovie', False)

        outputMovies = getattr(self, outputName, None)
        if outputMovies is None and not newDone:
            return

        outputMovies, created = self._getLogicalOutputSet(
            outputName,
            SetOfMovies,
            fixSampling=saveMovie,
        )

        if saveMovie:
            outputMovies.setGain(None)
            outputMovies.setDark(None)

        persistedOutputIds = self._getOutputItemIds(outputMovies)

        for movie in newDone:
            objId = movie.getObjId()
            if objId in persistedOutputIds:
                continue

            try:
                # Data-building failures for one Movie may be reported
                # locally, but persistence failures below must propagate.
                outMovie = self._createOutputMovie(movie)
                outMovie.copyObjId(movie)

                shifts = outMovie.getAlignment().getShifts()[0]
                if not shifts:
                    self.warning(
                        "Movie %s has empty alignment data; "
                        "it was not added to outputMovies."
                        % movie.getFileName()
                    )
                    continue
            except Exception as exc:
                self.error(
                    "Movie %s could not be prepared for outputMovies: %s"
                    % (movie.getFileName(), exc)
                )
                continue

            # Persistence is deliberately outside the exception boundary.
            # A failed append/update means the Movie is not durably
            # published and must not be treated as completed on resume.
            outputMovies.append(outMovie)
            outputMovies.update(outMovie)
            persistedOutputIds.add(objId)

        if created and newDone:
            self._storeSummary(newDone[0])
            if not saveMovie:
                outputMovies.setDim(
                    self.getInputMovies().getDim()
                )

        self._persistLogicalOutputSet(
            outputMovies,
            streamMode,
        )

    def _updateLogicalOutputMicSet(self, newDone, outputName,
                                   getOutputMicName, streamMode,
                                   suffix=''):
        outputMics = getattr(self, outputName, None)
        if outputMics is None and not newDone:
            return

        outputMics, _ = self._getLogicalOutputSet(
            outputName,
            SetOfMicrographs,
            suffix=suffix,
            fixSampling=True,
        )

        persistedOutputIds = self._getOutputItemIds(outputMics)

        for movie in newDone:
            objId = movie.getObjId()
            if objId in persistedOutputIds:
                continue

            mic = outputMics.ITEM_TYPE()
            mic.copyObjId(movie)
            mic.setMicName(movie.getMicName())

            micFn = self._getExtraPath(
                getOutputMicName(movie)
            )
            mic.setFileName(micFn)

            if not os.path.exists(micFn):
                self.warning(
                    "Micrograph %s was not generated; "
                    "it was not added to %s."
                    % (micFn, outputName)
                )
                continue

            try:
                self._preprocessOutputMicrograph(
                    mic,
                    movie,
                )
            except Exception as exc:
                self.error(
                    "Could not prepare %s for movie %s: %s"
                    % (
                        outputName,
                        movie.getFileName(),
                        exc,
                    )
                )
                continue

            outputMics.append(mic)
            outputMics.update(mic)
            persistedOutputIds.add(objId)

        self._persistLogicalOutputSet(
            outputMics,
            streamMode,
        )

    def _updateOutputSets(self, newDone, streamMode):
        if self._createOutputMovies():
            self._updateLogicalOutputMovieSet(
                newDone,
                streamMode,
            )

        if self._createOutputMicrographs():
            self._updateLogicalOutputMicSet(
                newDone,
                'outputMicrographs',
                self._getOutputMicName,
                streamMode,
            )

        if self._createOutputWeightedMicrographs():
            self._updateLogicalOutputMicSet(
                newDone,
                'outputMicrographsDoseWeighted',
                self._getOutputMicWtName,
                streamMode,
                suffix='dose-weighted',
            )

        if self._doSplitEvenOdd():
            self._updateLogicalOutputMicSet(
                newDone,
                'outputMicrographsEven',
                self._getOutputMicEvenName,
                streamMode,
                suffix='even',
            )
            self._updateLogicalOutputMicSet(
                newDone,
                'outputMicrographsOdd',
                self._getOutputMicOddName,
                streamMode,
                suffix='odd',
            )

    def _getMcProcessor(self, gpu):
        def _processBatch(batch):
            tries = 2
            while tries:
                tries -= 1
                try:
                    n = len(batch['items'])
                    t = Timer()

                    batch_path = batch['path']
                    batch_output = os.path.join(batch_path, 'output')

                    # The output folder may exists if re-trying
                    if not os.path.exists(batch_output):
                        Process.system(f"mkdir '{batch_output}'")

                    cmd = self.command.replace('-Gpu #', f'-Gpu {gpu}')
                    self.runJob(self.program, cmd, cwd=batch_path)

                    elapsed = t.getToc()
                    t.toc(f'Ran motioncor batch of {n} movies')
                    tries = 0  # Everything run OK, no more tries

                    with self.lock:
                        self.processed += n
                        thread_id = threading.get_ident()
                        self.debug(
                            f" {thread_id}: Processing "
                            f"{self.batch_str(batch)}, "
                            f"{elapsed}, Processed: {self.processed}"
                        )

                except Exception as e:
                    import traceback
                    traceback.print_exc()

                    if tries == 0:
                        self.error(
                            "ERROR: Motioncor has failed for batch %s after "
                            "the final retry. No more tries for this batch!!!"
                            % batch['id']
                        )

                        # IMPORTANT: do not raise from an emtools Pipeline
                        # worker thread. TaskGenerator.run() would terminate
                        # before notifyGeneratorEnds(), leaving downstream
                        # queues waiting forever. Record the failure and let
                        # the pipeline drain. The main Scipion step will
                        # re-raise it after mc.run().
                        with self.lock:
                            self._batchFailures.append(e)
                            self._failedBatchIds.add(batch['id'])
                        return batch

                    self.debug(
                        "ERROR: Motioncor has failed for batch %s. --> %s\n."
                        "Sleeping and re-trying in one minute."
                        % (batch['id'], str(e))
                    )
                    time.sleep(60)

            return batch

        return _processBatch

    def _getOutputProcessor(self):
        def _processBatchOutput(batch):
            try:
                return self._moveBatchOutput(batch)
            except Exception as exc:
                # emtools Pipeline runs processors in background threads.
                # Let the queue drain, but preserve the publication failure
                # so _processAllMoviesStep can raise it in the Scipion step.
                with self.lock:
                    self._batchFailures.append(exc)
                return batch

        return _processBatchOutput

    def _moveBatchOutput(self, batch):
        # Failed batches are forwarded only to let the emtools Pipeline
        # drain and close. Never publish outputs for them.
        if batch['id'] in self._failedBatchIds:
            self.debug(
                "Skipping output registration for failed batch %s."
                % batch['id']
            )
            return batch

        t = Timer()
        srcDir = batch['path']
        doClean = not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN)
        applyDose = self._createOutputWeightedMicrographs()
        saveUnweighted = self._doSaveUnweightedMic()
        usePatches = self.patchX != 0 or self.patchY != 0
        logSuffix = '%s-Full.log' % ('-Patch' if usePatches else '')
        newDone = []
        missing = {}

        def _moveToExtra(movie, src, dst):
            srcFn = os.path.join(srcDir, src)
            dstFn = self._getExtraPath(dst)
            if os.path.exists(srcFn):
                pwutils.moveFile(srcFn, dstFn)
                return True
            self.debug(f"Missing file: {srcFn}")
            missing[movie.getObjId()] = movie
            return False

        def _moveMovieFiles(movie):
            movieRoot = 'output/' + ProtMotionCorr._getMovieRoot(self, movie)
            self.debug(f"Moving output for movie: {movieRoot}")

            if applyDose:
                _moveToExtra(movie, movieRoot + '_DW.mrc', self._getOutputMicWtName(movie))

            if not applyDose or saveUnweighted:
                _moveToExtra(movie, movieRoot + '.mrc', self._getOutputMicName(movie))

            if self.splitEvenOdd:
                _moveToExtra(movie, movieRoot + '_EVN.mrc',
                             self._getOutputMicEvenName(movie))
                _moveToExtra(movie, movieRoot + '_ODD.mrc',
                             self._getOutputMicOddName(movie))

            _moveToExtra(movie, movieRoot + logSuffix, self._getMovieLogFile(movie))

        for movie in batch['items']:
            _moveMovieFiles(movie)
            if movie.getObjId() not in missing:
                newDone.append(movie)

        self.debug(f" Moving {self.batch_str(batch)}, "
                   f"newDone: {len(newDone)}, missing: {len(missing)}")

        if newDone:
            self._firstTimeOutput = not hasattr(self, 'outputMovies')
            self.debug(f">>> Updating outputs, newDone: {len(newDone)}, "
                       f"firstTimeOutput: {self._firstTimeOutput}")
            self._updateOutputSets(newDone, pwobj.Set.STREAM_OPEN)

        elapsed = t.getToc()
        t.toc('Registered outputs')

        with self.lock:
            self.registered += len(newDone)
            self.debug(f"OUTPUT: {self.batch_str(batch)}, "
                       f"{elapsed}, "
                       f"New done {len(newDone)}, "
                       f"Registered {self.registered}, "
                       f"Processed {self.processed}")
            for movie in missing.values():
                self.debug(f"FAILED: {movie.getFileName()}")

        # Clean batch folder if not in debug mode
        if doClean:
            Process.system('rm -rf %s' % batch['path'])

        t.toc(f"Moved output for batch {batch['id']}")

        return batch

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getCmd(self):
        """ Set return a command string that will be used for each batch. """
        inputMovies = self.getInputMovies()
        argsDict = self._getMcArgs()
        argsDict['-Gpu'] = '#'
        argsDict['-Serial'] = 1
        argsDict['-LogDir'] = "output/"

        # Get input format, but for the batch. In live streaming the
        # logical input may legitimately still be empty when Tasks starts.
        # Refresh it until the first Movie arrives instead of dereferencing
        # None or relying on a physical backing file.
        firstMovie = inputMovies.getFirstItem()
        while firstMovie is None and inputMovies.isStreamOpen():
            time.sleep(self.streamingSleepOnWait.get())
            inputMovies.loadAllProperties()
            firstMovie = inputMovies.getFirstItem()

        if firstMovie is None:
            raise RuntimeError(
                "No movies available and the input stream is already "
                "closed; nothing to process."
            )

        ext = pwutils.getExt(firstMovie.getFileName()).lower()
        if ext in ['.mrc', '.mrcs']:
            inprefix = '-InMrc'
        elif ext in ['.tif', '.tiff']:
            inprefix = '-InTiff'
        elif ext in ['.eer']:
            inprefix = '-InEer'
        else:
            raise Exception(f"Unsupported format '{ext}' for batch processing "
                            f"in Motioncor protocol. ")

        argsDict[inprefix] = './'
        # In serial mode MotionCor scans the whole input directory. Without
        # a suffix filter, non-movie entries in the batch directory (notably
        # output/) can be treated as input stacks and cause header errors.
        argsDict['-InSuffix'] = ext
        argsDict['-OutMrc'] = 'output/'

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += self.extraParams2.get()

        return cmd

    def _createOutputWeightedMicrographs(self):
        # Keep every Tasks code path aligned with the actual MotionCor
        # command: requesting dose weighting in the form is not enough;
        # usable dose metadata must also be available.
        return ProtMotionCorr._useDoseWeightedOutput(self)

    def _setPlotInfo(self, movie, mic):
        # FIXME: For now not support PSD or Thumbnail
        if self._createOutputWeightedMicrographs():
            total, early, late = self.calcFrameMotion(movie)
            mic._rlnAccumMotionTotal = Float(total)
            mic._rlnAccumMotionEarly = Float(early)
            mic._rlnAccumMotionLate = Float(late)

    def _getMovieRoot(self, movie):
        return "mic_%06d" % movie.getObjId()

    def _getOutputMovieName(self, movie):
        """ Returns the name of the output movie.
        (relative to micFolder)
        """
        return "movie_%06d" % movie.getObjId()

    def _getOutputMicName(self, movie):
        """ Returns the name of the output micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '.mrc'

    def _getOutputMicWtName(self, movie):
        """ Returns the name of the output dose-weighted micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_DW.mrc'

    def _getOutputMicEvenName(self, movie):
        """ Returns the name of the output EVEN micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_EVN.mrc'

    def _getOutputMicOddName(self, movie):
        """ Returns the name of the output EVEN micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_ODD.mrc'

    def _getOutputMicThumbnail(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_thumbnail.png')

    def _getMovieLogFile(self, movie):
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s%s-Full.log' % (self._getMovieRoot(movie),
                                  '-Patch' if usePatches else '')

    def debug(self, msg):
        self.error(f"{Pretty.now()}: DEBUG >>> {msg}")

    def batch_str(self, batch):
        batch_ids = [m.getObjId() for m in batch['items']]
        return f"Batch {batch['index']}:{batch_ids}"
