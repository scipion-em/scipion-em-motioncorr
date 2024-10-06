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
from emtools.pwx import SetMonitor, BatchManager

from pyworkflow import SCIPION_DEBUG_NOCLEAN, BETA
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.protocol import STEPS_SERIAL
from pwem.objects import Float, SetOfMovies

from .. import Plugin
from .protocol_motioncorr import ProtMotionCorr


class ProtMotionCorrTasks(ProtMotionCorr):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'tasks'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtMotionCorr.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
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
    def _insertAllSteps(self):
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        self._insertFunctionStep(self._convertInputStep)
        # Make the following step to always run, despite finished
        # this may be useful when new input items (from streaming)
        # and need to continue
        self._insertFunctionStep(self._processAllMoviesStep, Pretty.now())

    def _processAllMoviesStep(self, when):
        self.lock = threading.Lock()
        self.processed = 0
        self.registered = 0

        self.error(f">>> {Pretty.now()}: ----------------- "
                   f"Start processing movies----------- ")
        self.program = Plugin.getProgram()
        self.command = self._getCmd()
        self._firstTimeOutput = True

        moviesMtr = SetMonitor(SetOfMovies,
                               self.getInputMovies().getFileName(),
                               blacklist=getattr(self, 'outputMovies', None))
        waitSecs = self.streamingSleepOnWait.get()
        moviesIter = moviesMtr.iterProtocolInput(self, 'movies', waitSecs=waitSecs)
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

        o1 = mc.addProcessor(outputQueue, self._moveBatchOutput)
        mc.run()
        # Mark the output as closed
        self._firstTimeOutput = False
        self._updateOutputSets([], pwobj.Set.STREAM_CLOSED)

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
                        self.debug(f" {thread_id}: Processing {self.batch_str(batch)}, "
                                   f"{elapsed}, Processed: {self.processed}")

                except Exception as e:
                    self.debug("ERROR: Motioncor has failed for batch %s. --> %s\n."
                               "Sleeping and re-trying in one minute." % (batch['id'], str(e)))
                    time.sleep(60)
                    import traceback
                    traceback.print_exc()
                    if tries == 0:
                        self.error("ERROR: No more tries for this batch!!!")

            return batch

        return _processBatch

    def _moveBatchOutput(self, batch):
        t = Timer()
        srcDir = batch['path']
        doClean = not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN)
        applyDose = self.doApplyDoseFilter
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

        # Get input format, but for the batch
        firstMovie = inputMovies.getFirstItem()
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
        argsDict['-OutMrc'] = 'output/'

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += self.extraParams2.get()

        return cmd

    def _setPlotInfo(self, movie, mic):
        # FIXME: For now not support PSD or Thumbnail
        if self.doApplyDoseFilter:
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
