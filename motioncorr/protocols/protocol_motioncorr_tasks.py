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

from collections import OrderedDict

from emtools.utils import Timer, Pipeline, Pretty
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
        self._insertFunctionStep(self._processAllMoviesStep)

    def _processAllMoviesStep(self):
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
        for gpu in gpus:
            p = mc.addProcessor(g.outputQueue, self._getMcProcessor(gpu),
                                outputQueue=outputQueue)
            outputQueue = p.outputQueue

        o1 = mc.addProcessor(outputQueue, self._moveBatchOutput)
        mc.run()
        # Mark the output as closed
        self._updateOutputSets([], pwobj.Set.STREAM_CLOSED)

    def _getMcProcessor(self, gpu):
        def _processBatch(batch):
            tries = 2
            while tries:
                tries -= 1
                try:
                    n = len(batch['items'])
                    t = Timer()

                    cmd = self.command.replace('-Gpu #', f'-Gpu {gpu}')
                    self.runJob(self.program, cmd, cwd=batch['path'])

                    elapsed = t.getToc()
                    t.toc(f'Ran motioncor batch of {n} movies')

                    with self.lock:
                        self.processed += n
                        batch_ids = [m.getObjId() for m in batch['items']]
                        self.error(f"Batch {batch['index']}:{batch_ids}, "
                                   f"{elapsed}, Processed: {self.processed}")

                except Exception as e:
                    self.error("ERROR: Motioncor2 has failed for batch %s. --> %s\n"
                               % (batch['id'], str(e)))
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

        def _moveToExtra(src, dst):
            srcFn = os.path.join(srcDir, src)
            dstFn = self._getExtraPath(dst)
            self.error(f"Moving {srcFn} -> {dstFn}"
                       f"\n\texists source: {os.path.exists(srcFn)}")
            if os.path.exists(srcFn):
                pwutils.moveFile(srcFn, dstFn)
                return True
            return False

        def _moveMovieFiles(movie):
            movieRoot = 'output_' + ProtMotionCorr._getMovieRoot(self, movie)
            self.error(f"Moving output for movie: {movieRoot}")

            if applyDose:
                a = _moveToExtra(movieRoot + '_DW.mrc', self._getOutputMicWtName(movie))

            if not applyDose or saveUnweighted:
                b = _moveToExtra(movieRoot + '.mrc', self._getOutputMicName(movie))

            if self.splitEvenOdd:
                _moveToExtra(movieRoot + '_EVN.mrc',
                             self._getOutputMicEvenName(movie))
                _moveToExtra(movieRoot + '_ODD.mrc',
                             self._getOutputMicOddName(movie))

            _moveToExtra(movieRoot + logSuffix, self._getMovieLogFile(movie))
            if a or b:
                newDone.append(movie)

        for movie in batch['items']:
            _moveMovieFiles(movie)

        if newDone:
            self._firstTimeOutput = not hasattr(self, 'outputMovies')
            self.error(f">>> Updating outputs, newDone: {len(newDone)}, "
                       f"firstTimeOutput: {self._firstTimeOutput}")
            self._updateOutputSets(newDone, pwobj.Set.STREAM_OPEN)

        elapsed = t.getToc()
        t.toc('Registered outputs')

        with self.lock:
            self.registered += len(newDone)
            batch_ids = [m.getObjId() for m in batch['items']]
            self.error(f"OUTPUT: Batch {batch['index']}:{batch_ids}, "
                       f"{elapsed}, "
                       f"New done {len(newDone)}, "
                       f"Registered {self.registered}, "
                       f"Processed {self.processed}")

        # Clean batch folder if not in debug mode
        if doClean:
            os.system('rm -rf %s' % batch['path'])

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

        # FIXME: Duplicated from ProtMotionCorr._processMovie
        if self.isEER:
            argsDict['-EerSampling'] = self.eerSampling.get() + 1
            argsDict['-FmIntFile'] = "../../extra/FmIntFile.txt"

        if self.doApplyDoseFilter:
            preExp, dose = self._getCorrectedDose(inputMovies)
            argsDict.update({'-FmDose': dose,
                             '-InitDose': preExp if preExp > 0.001 else 0})

        if inputMovies.getGain():
            argsDict.update({'-Gain': f'"{inputMovies.getGain()}"',
                             '-RotGain': self.gainRot.get(),
                             '-FlipGain': self.gainFlip.get()})

        if inputMovies.getDark():
            argsDict['-Dark'] = inputMovies.getDark()

        # Get input format, but for the batch
        firstMovie = inputMovies.getFirstItem()
        ext = pwutils.getExt(firstMovie.getFileName()).lower()
        if ext in ['.mrc', '.mrcs']:
            inprefix = '-InMrc'
        elif ext in ['.tif', '.tiff']:
            inprefix = '-InTiff'
        else:
            raise Exception(f"Unsupported format '{ext}' for batch processing "
                            f"in Motioncor protocol. ")

        argsDict[inprefix] = './'
        argsDict['-OutMrc'] = 'output_'

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
