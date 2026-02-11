# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# **************************************************************************
import logging
import os
from os.path import basename, join, abspath
import pyworkflow.utils as pwutils
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam, BooleanParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils import cyanStr, makePath, removeBaseExt, Message
from tomo.objects import TiltImageM
# from tomo.protocols import ProtTsCorrectMotion
from .. import Plugin
from .protocol_base import ProtMotionCorrBase


logger = logging.getLogger(__name__)


class ProtTsMotionCorr(ProtMotionCorrBase):#, ProtTsCorrectMotion):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'align tilt-series movies'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.doSaveMovie = False
        self.doApplyDoseFilter = False
        self.tsDict = None
        self.patchStr = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, addEvenOddParam = True):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTiltSeriesM', PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      important=True,
                      label='Input Tilt-Series (movies)',
                      help='Select input tilt-series movies that you want'
                           'to correct for beam-induced motion. ')

        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN',
                             help='Frames range to ALIGN on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie.')
        line.addParam('alignFrame0', IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', IntParam, default=0,
                      label='to')
        group.addParam('useAlignToSum', BooleanParam, default=True,
                       label='Use ALIGN frames range to SUM?',
                       help="If *Yes*, the same frame range will be used to "
                            "ALIGN and to SUM. If *No*, you can selected a "
                            "different range for SUM (must be a subset).")
        line = group.addLine('Frames to SUM', condition="not useAlignToSum",
                             help='Frames range to SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to sum, it means that you will sum '
                                  'until the last frame of the movie.')
        line.addParam('sumFrame0', IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', IntParam, default=0,
                      label='to')
        group.addParam('binFactor', FloatParam, default=1.,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        line = group.addLine('Crop offsets (px)',
                             expertLevel=LEVEL_ADVANCED)
        line.addParam('cropOffsetX', IntParam, default=0, label='X')
        line.addParam('cropOffsetY', IntParam, default=0, label='Y')

        line = group.addLine('Crop dimensions (px)',
                             expertLevel=LEVEL_ADVANCED,
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', IntParam, default=0, label='X')
        line.addParam('cropDimY', IntParam, default=0, label='Y')

        if self.evenOddCapable and addEvenOddParam:
            form.addParam('splitEvenOdd', BooleanParam,
                          default=False,
                          label='Split & sum odd/even frames?',
                          help='(Used for denoising data preparation). If set to Yes, 2 additional movies/tilt '
                               'series will be generated, one generated from the even frames and the other from the '
                               'odd ones using the same alignment for the whole stack of frames.')

        self._defineCommonParams(form, allowDW=False)
        # Patch alignment is not recommended
        form.getParam('patchX').setDefault(0)
        form.getParam('patchY').setDefault(0)
        form.addParallelSection(threads=4, mpi=1)


    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self._convertInputStep, needsGPU=False)
        for tsId, tsM in self.tsMDict.items():
            makePath(self.getTsResultsPath(tsId))
            for counter, tiM in enumerate(tsM.iterItems()):
                tiM = tiM.clone()
                self._insertFunctionStep(self.processAngularStackStep,
                                         tsId,
                                         tiM,
                                         counter + 1,
                                         needsGPU=True)


        # self._insertFunctionStep(self.runDataExtraction, needsGPU=False)
        # self._insertFunctionStep(self.prepareTrainingStep, needsGPU=False)
        # self._insertFunctionStep(self.trainingStep, needsGPU=True)
        # self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsMDict = {tsM.getTsId(): tsM.clone() for tsM in self.inputTiltSeriesM.get().iterItems()}
        usePatches = self.patchX != 0 or self.patchY != 0
        self.patchStr = '-Patch' if usePatches else ''

    #     inputTs = self._getInputTs()
    #     acq = inputTs.getAcquisition()
    #     gain, dark = self.getGainAndDark()
    #     self.__basicArgs = [
    #         acq.getDoseInitial(), acq.getDosePerFrame(), gain, dark]


    def processAngularStackStep(self, tsId: str, tiM: TiltImageM, counter: int):
        tiMFileName = tiM.getFileName()
        logger.info(cyanStr(f'tsId = {tsId} - processing angular stack {basename(tiMFileName)}'))
        argsDict = self._getMcArgs(acqOrder=tiM.getAcquisitionOrder())
        argsDict['-OutMrc'] = f'"{abspath(self.getOutTsFName(tsId, counter))}"'
        argsDict['-LogDir'] = f'"{abspath(self.getTsResultsPath(tsId))}"'
        params = self._getInputFormat(tiMFileName, absPath=True)
        params += ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        params += ' ' + self.extraParams2.get()

        try:
            self.runJob(self.program, params, env=Plugin.getEnviron())

        except Exception as e:
            self.error(f"ERROR: Motioncor has failed for {tiMFileName} --> {str(e)}\n")
            import traceback
            traceback.print_exc()

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        tiltImageM = self._tsDict.getTi(tsId, tiltImageId)
        workingFolder = self._getTmpPath(self._getTiltImageMRoot(tiltImageM))
        pwutils.makePath(workingFolder)
        self._processTiltImageM(workingFolder, tiltImageM, *args)

        if self._doSplitEvenOdd():
            baseName = self._getTiltImageMRoot(tiltImageM)

            evenName = (os.path.abspath(self._getExtraPath(baseName + '_EVN.mrc')))
            oddName = (os.path.abspath(self._getExtraPath(baseName + '_ODD.mrc')))

            # Store the corresponding tsImM to use its data later in the even/odd TS
            tiltImageM.setOddEven([oddName, evenName])

        tiFn, _ = self._getOutputTiltImagePaths(tiltImageM)
        if not os.path.exists(tiFn):
            raise FileNotFoundError(f"Expected output file '{tiFn}' not produced!")

        if not pwutils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pwutils.cleanPath(workingFolder)

    # def _processTiltImageM(self, workingFolder, tiltImageM, *args):
    #     outputFn, _ = self._getOutputTiltImagePaths(tiltImageM)
    #
    #     def _getPath(path):
    #         """ shortcut to get relative path from workingFolder. """
    #         return os.path.relpath(path, workingFolder)
    #
    #     self.info(f"workingFolder: {workingFolder}")
    #     self.info(f"outputFn: {outputFn}")
    #
    #     argsDict = self._getMcArgs(acqOrder=tiltImageM.getAcquisitionOrder())
    #     argsDict['-OutMrc'] = f'"{_getPath(outputFn)}"'
    #
    #     tiFn = tiltImageM.getFileName()
    #     inputFn = os.path.abspath(tiFn)
    #
    #     self.info(f"inputFn: {tiFn}")
    #
    #     params = self._getInputFormat(inputFn, absPath=True)
    #     params += ' '.join(['%s %s' % (k, v)
    #                         for k, v in argsDict.items()])
    #
    #     params += ' ' + self.extraParams2.get()
    #
    #     try:
    #         self.runJob(Plugin.getProgram(), params,
    #                     cwd=workingFolder,
    #                     env=Plugin.getEnviron())
    #
    #         # Move output log to extra dir
    #         logFn = os.path.join(workingFolder, self._getMovieLogFile(tiltImageM))
    #         logFnExtra = self._getExtraPath(self._getMovieLogFile(tiltImageM))
    #         pwutils.moveFile(logFn, logFnExtra)
    #
    #     except Exception as e:
    #         self.error(f"ERROR: Motioncor has failed for {tiFn} --> {str(e)}\n")
    #         import traceback
    #         traceback.print_exc()



    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'TiltSeries'):
            summary.append('Aligned %d tilt series movies using Motioncor.'
                           % self.getInputMovies().getSize())
        else:
            summary.append('Output is not ready')

        return summary

    def _validate(self):
        errors = ProtMotionCorrBase._validate(self)

        return errors

    # --------------------------- UTILS functions -----------------------------
    def getInputMovies(self):
        return self.inputTiltSeriesM.get()

    def _getMovieLogFile(self, tiMFileName: str) -> str:
        return f'{removeBaseExt(tiMFileName)}{self.patchStr}-Full.log'

    @staticmethod
    def _createOutputWeightedTS():
        return False

    @staticmethod
    def _getOutputName():
        return 'TiltSeries'

    def getTsResultsPath(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def getTiMLogDir(self, tsId: str, counter: int) -> str:
        return self._getTmpPath(f'{tsId}_{counter:02d}')

    def getOutTsFName(self, tsId: str, counter: int) -> str:
        return join(self.getTsResultsPath(tsId), f'{tsId}_{counter:02d}.mrc')


