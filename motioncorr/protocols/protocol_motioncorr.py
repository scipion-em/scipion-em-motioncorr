# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Vahid Abrishami (vabrishami@cnb.csic.es) [2]
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [3]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [4]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [3] Department of Anatomy and Cell Biology, McGill University
# * [4] MRC Laboratory of Molecular Biology (MRC-LMB)
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
import logging
import time
import traceback
from collections import Counter
from enum import Enum
from math import ceil, sqrt
from os.path import exists, basename, join, abspath
from threading import Thread
import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.plotter import Plotter
from pwem.objects import Image, Float, SetOfMovies, SetOfMicrographs, Movie
from pwem.protocols import ProtAlignMovies
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.utils import cyanStr, makePath, Message, redStr, removeBaseExt
from .. import Plugin
from .protocol_base import ProtMotionCorrBase
from relion.convert.convert31 import OpticsGroups
from ..convert import parseMovieAlignment2

logger = logging.getLogger(__name__)


class MotionCorrOutputs(Enum):
    movies = SetOfMovies()
    micrographs = SetOfMicrographs()
    micrographsDW = SetOfMicrographs()
    micrographsEven = SetOfMicrographs()
    micrographsOdd = SetOfMicrographs()


class ProtMotionCorrNewStreaming(ProtMotionCorrBase, ProtStreamingBase):  # , ProtAlignMovies):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'movie alignment NS'
    evenOddCapable = True
    _possibleOutputs = MotionCorrOutputs
    stepsExecutionMode = cons.STEPS_PARALLEL

    def __init__(self, **kwargs):
        # ProtAlignMovies.__init__(self, **kwargs)
        super().__init__(**kwargs)
        self.itemIdReadList = []
        self.sRate = None
        self.gain = None
        self.dark = None
        self.sRate = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported movies.')
        self._defineAlignmentParams(form)
        self._defineCommonParams(form)

    def _defineAlignmentParams(self, form):
        form.addHidden('doSaveAveMic', params.BooleanParam,
                       default=True)
        form.addHidden('useAlignToSum', params.BooleanParam,
                       default=True)

        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN and SUM',
                             help='Frames range to ALIGN and SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie. '
                                  'When using EER, this option is IGNORED!')
        line.addParam('alignFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', params.IntParam, default=0,
                      label='to')

        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        line = group.addLine('Crop offsets (px)',
                             expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X',
                      expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y',
                      expertLevel=cons.LEVEL_ADVANCED)

        line = group.addLine('Crop dimensions (px)',
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.',
                             expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropDimX', params.IntParam, default=0, label='X',
                      expertLevel=cons.LEVEL_ADVANCED)
        line.addParam('cropDimY', params.IntParam, default=0, label='Y',
                      expertLevel=cons.LEVEL_ADVANCED)

        if self.evenOddCapable:
            form.addParam('splitEvenOdd', params.BooleanParam,
                          default=False,
                          label='Split & sum odd/even frames?',
                          help='Generate odd and even sums using odd and even frames '
                               'respectively when this option is enabled.')

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Save aligned movie?")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Compute PSD?",
                      help="If Yes, the protocol will compute for each "
                           "aligned micrograph the PSD using EMAN2.")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      expertLevel=cons.LEVEL_ADVANCED,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail with EMAN2 and keep it with the '
                           'micrograph object for visualization purposes. ')

        form.addParam('extraProtocolParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional protocol parameters',
                      help="Here you can provide some extra parameters for the "
                           "protocol, not the underlying motioncor program."
                           "You can provide many options separated by space. "
                           "\n\n*Options:* \n\n"
                           "--dont_use_worker_thread \n"
                           " Now by default we use a separate thread to compute"
                           " PSD and thumbnail (if is required). This allows "
                           " more effective use of the GPU card, but requires "
                           " an extra CPU. Use this option (NOT RECOMMENDED) if "
                           " you want to prevent this behaviour")

    # --------------------------- STEPS functions -----------------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        self._initialize()
        inMoviesSet = self.getInputMovies()
        self.sRate = inMoviesSet.getSamplingRate()
        self.readingOutput()

        while True:
            with self._lock:
                inIds = set(inMoviesSet.getUniqueValues('id'))

            # In the if statement below, Counter is used because in the objId comparison the order doesn’t matter
            # but duplicates do. With a direct comparison, the closing step may not be inserted because of the order:
            # ['id_a', 'id_b'] != ['id_b', 'id_a'], but they are the same with Counter.
            if not inMoviesSet.isStreamOpen() and Counter(self.itemIdReadList) == Counter(inIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self._closeOutputSet,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedIds = inIds - set(self.itemIdReadList)
            moviesToProcessDict = {objId: movie.clone() for movie in inMoviesSet.iterItems()
                                   if (objId := movie.getObjId()) in nonProcessedIds}
            for objId, movie in moviesToProcessDict.items():
                cInPId = self._insertFunctionStep(self.convertInputStep,
                                                  prerequisites=[],
                                                  needsGPU=False)
                pMovPid = self._insertFunctionStep(self.processMovieStep,
                                                   objId,
                                                   movie,
                                                   prerequisites=cInPId,
                                                   needsGPU=True)
                cOutId = self._insertFunctionStep(self.createOutputStep,
                                                  prerequisites=pMovPid,
                                                  needsGPU=False)
                closeSetStepDeps.append(cOutId)
                logger.info(cyanStr(f"Steps created for objId = {objId} - {movie.getFileName()}"))
                self.itemIdReadList.append(objId)

            time.sleep(10)
            if inMoviesSet.isStreamOpen():
                with self._lock:
                    inMoviesSet.loadAllProperties()  # refresh status for the streaming

    def readingOutput(self) -> None:
        movies = getattr(self, self._possibleOutputs.movies.name, None)
        outputList = []
        mics = getattr(self, self._possibleOutputs.micrographs.name, None)
        micsDW = getattr(self, self._possibleOutputs.micrographsDW.name, None)
        micsEven = getattr(self, self._possibleOutputs.micrographsEven.name, None)
        micsOdd = getattr(self, self._possibleOutputs.micrographsOdd.name, None)
        if self.splitEvenOdd:
            outputList.append(micsEven)
            outputList.append(micsOdd)
        if self.doApplyDoseFilter:
            outputList.append(micsDW)
        else:
            outputList.append(mics)

        if None not in outputList:
            for item in movies:
                self.itemIdReadList.append(item.getObjId())
            self.info(cyanStr(f'Ids processed: {self.itemIdReadList}'))
        else:
            self.info(cyanStr('No movies have been processed yet'))

    def _initialize(self):
        inputMovies = self.getInputMovies()
        self.gain = inputMovies.getGain()
        self.dark = inputMovies.getDark()
        self.sRate = inputMovies.getSamplingRate()

    def convertInputStep(self):
        # Convert the gain and dark images if they haven't been converted yet or if they weren't provided
        if self.dark or self.gain and not exists(self._getExtraPath('DONE')):
            super()._convertInputStep()

    def processMovieStep(self, objId: int, movie: Movie):
        # movieFolder = self._getOutputMovieFolder(objId)
        # movieFn = movie.getFileName()
        # movieName = basename(movieFn)
        # makePath(movieFolder)
        # createAbsLink(abspath(movieFn), join(movieFolder, movieName))

        # if movieName.endswith('bz2'):
        #     newMovieName = movieName.replace('.bz2', '')
        #     # We assume that if compressed the name ends with .mrc.bz2
        #     if not exists(newMovieName):
        #         self.runJob('bzip2', '-d -f %s' % movieName, cwd=movieFolder)
        # 
        # elif movieName.endswith('tbz'):
        #     newMovieName = movieName.replace('.tbz', '.mrc')
        #     # We assume that if compressed the name ends with .tbz
        #     if not exists(newMovieName):
        #         self.runJob('tar', 'jxf %s' % movieName, cwd=movieFolder)

        # elif movieName.endswith('.txt'):
        #     # Support a list of frame as a simple .txt file containing
        #     # all the frames in a raw list, we could use a xmd as well,
        #     # but a plain text was choose to simply its generation
        #     movieTxt = join(movieFolder, movieName)
        #     with open(movieTxt) as f:
        #         movieOrigin = basename(movieFn)
        #         newMovieName = movieName.replace('.txt', '.mrcs')
        #         ih = emlib.image.ImageHandler()
        #         for i, line in enumerate(f):
        #             if line.strip():
        #                 inputFrame = join(movieOrigin, line.strip())
        #                 ih.convert(inputFrame,
        #                            (i + 1, join(movieFolder, newMovieName)))
        # else:
        #     newMovieName = movieName

        # Just store the original name in case it is needed in _processMovie
        # movie._originalFileName = pwobj.String(objDoStore=False)
        # movie._originalFileName.set(movie.getFileName())
        # Now set the new filename (either linked or converted)
        # movie.setFileName(join(movieFolder, newMovieName))
        movieFName = movie.getFileName()
        logger.info(cyanStr(f"Processing movie: {movieFName}"))

        # self._processMovie(objId, movie)
        outputMicFn = self._getResultMicFn(movie)

        argsDict = self._getMcArgs()
        argsDict['-OutMrc'] = f'"{outputMicFn}"'
        argsDict['-LogDir'] = f'"{self._getExtraPath()}"'

        args = self._getInputFormat(movieFName, absPath=True)
        args += ' '.join(['%s %s' % (k, v)
                          for k, v in argsDict.items()])
        args += ' ' + self.extraParams2.get()

        try:
            self.runJob(Plugin.getProgram(), args, env=Plugin.getEnviron())
        except Exception as e:
            logger.error(redStr(f"ERROR: Motioncor has failed for {movie.getFileName()} "
                                f"--> {str(e)}\n"))
            traceback.print_exc()

        # we need to move shifts log to extra dir before parsing
        try:
            self._saveAlignmentPlots(outputMicFn, self.sRate)
            outMicFn = self._getExtraPath(self._getMicFn(movie))

            if self.doComputePSD:
                self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))

            if self._doComputeMicThumbnail():
                self.computeThumbnail(
                    outMicFn, outputFn=self._getOutputMicThumbnail(movie))
        except Exception as e:
            logger.error(redStr(f"ERROR: Extra work (i.e plots, PSD, thumbnail) "
                                f"has failed for {movie.getFileName()} --> {str(e)}\n"))
            traceback.print_exc()

            # self._moveOutput(movie)

            # def _extraWork():
            #     # we need to move shifts log to extra dir before parsing
            #     try:
            #         self._saveAlignmentPlots(movie, self.sRate)
            #         outMicFn = self._getExtraPath(self._getMicFn(movie))
            #
            #         if self.doComputePSD:
            #             self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))
            #
            #         if self._doComputeMicThumbnail():
            #             self.computeThumbnail(
            #                 outMicFn, outputFn=self._getOutputMicThumbnail(movie))
            #     except:
            #         self.error(f"ERROR: Extra work (i.e plots, PSD, thumbnail) "
            #                    f"has failed for {movie.getFileName()}\n")
            #
            # if self._useWorkerThread():
            #     thread = Thread(target=_extraWork)
            #     thread.start()
            # else:
            #     _extraWork()



        # if self._doMovieFolderCleanUp():
        #     self._cleanMovieFolder(movieFolder)

    def createOutputStep(self):
        pass

    # def _processMovie(self, objId: int, movie: Movie):
    #     # movieFolder = self._getOutputMovieFolder(objId)
    #     outputMicFn = self._getResultMicFn(movie)
    #
    #     argsDict = self._getMcArgs()
    #     argsDict['-OutMrc'] = f'"{outputMicFn}"'
    #     argsDict['-LogDir'] = './'
    #
    #     args = self._getInputFormat(movie.getFileName())
    #     args += ' '.join(['%s %s' % (k, v)
    #                       for k, v in argsDict.items()])
    #     args += ' ' + self.extraParams2.get()
    #
    #     try:
    #         self.runJob(Plugin.getProgram(), args, cwd=movieFolder,
    #                     env=Plugin.getEnviron())
    #         self._moveOutput(movie)
    #
    #         def _extraWork():
    #             # we need to move shifts log to extra dir before parsing
    #             try:
    #                 self._saveAlignmentPlots(movie, self.sRate)
    #                 outMicFn = self._getExtraPath(self._getMicFn(movie))
    #
    #                 if self.doComputePSD:
    #                     self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))
    #
    #                 if self._doComputeMicThumbnail():
    #                     self.computeThumbnail(
    #                         outMicFn, outputFn=self._getOutputMicThumbnail(movie))
    #             except:
    #                 self.error(f"ERROR: Extra work (i.e plots, PSD, thumbnail) "
    #                            f"has failed for {movie.getFileName()}\n")
    #
    #         if self._useWorkerThread():
    #             thread = Thread(target=_extraWork)
    #             thread.start()
    #         else:
    #             _extraWork()
    #
    #     except Exception as e:
    #         self.error(f"ERROR: Motioncor has failed for {movie.getFileName()} --> {str(e)}\n")
    #         import traceback
    #         traceback.print_exc()
    #
    # def _insertFinalSteps(self, deps):
    #     stepsId = []
    #     if self._useWorkerThread():
    #         stepsId.append(self._insertFunctionStep('waitForThreadStep',
    #                                                 prerequisites=deps,
    #                                                 needsGPU=False))
    #     return stepsId
    #
    # def waitForThreadStep(self):
    #     # Quick and dirty (maybe desperate) way to wait
    #     # if the PSD and thumbnail were computed in a thread
    #     # If running in streaming this will not be necessary
    #     time.sleep(10)  # wait 10 sec for the thread to finish

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'outputMicrographs') or \
                hasattr(self, 'outputMicrographsDoseWeighted'):
            summary.append('Aligned %d movies using Motioncor.'
                           % self.getInputMovies().getSize())
            if self.splitEvenOdd and self._createOutputWeightedMicrographs():
                summary.append('Even/odd outputs are dose-weighted!')
        else:
            summary.append('Output is not ready')

        return summary

    def _validate(self):
        errors = ProtMotionCorrBase._validate(self)

        # check eman2 plugin
        if self.doComputeMicThumbnail or self.doComputePSD:
            try:
                from pwem import Domain
                _ = Domain.importFromPlugin('eman2', doRaise=True)
            except:
                errors.append("EMAN2 plugin not found!\nComputing thumbnails "
                              "or PSD requires EMAN2 plugin and binaries installed.")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getOutputMovieFolder(self, objId: int) -> str:
        """ Create a Movie folder where to work with it. """
        return self._getTmpPath('movie_%06d' % objId)
    
    def _getCwdPath(self, movie, path):
        return join(self._getOutputMovieFolder(movie), path)

    def _getMovieLogFile(self, movieFName: str):
        usePatches = self.patchX != 0 or self.patchY != 0
        pattern = '-Patch' if usePatches else ''
        return self._getExtraPath(f'{removeBaseExt(movieFName)}-{pattern}-Full.log')

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movieFName: str):
        return self._getNameExt(movieFName, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd', 'png', extra=True)

    def getInputMovies(self):
        return self.inputMovies.get()

    def _computePSD(self, inputFn, outputFn, scaleFactor=6):
        """ Generate a thumbnail of the PSD with EMAN2"""
        args = f"{inputFn} {outputFn} "
        args += f"--process=math.realtofft --meanshrink {scaleFactor} "
        args += "--fixintscaling=sane"

        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2')
        from pyworkflow.utils.process import runJob
        runJob(self._log, eman2.Plugin.getProgram('e2proc2d.py'), args,
               env=eman2.Plugin.getEnviron())

        return outputFn

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)

    def _getMicFn(self, movie):
        if self.doApplyDoseFilter:
            return self._getOutputMicWtName(movie)
        else:
            return self._getOutputMicName(movie)

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = Image(location=self._getPsdCorr(movie))
        if self._doComputeMicThumbnail():
            mic.thumbnail = Image(location=self._getOutputMicThumbnail(movie))
        if self.doApplyDoseFilter:
            total, early, late = self.calcFrameMotion(movie)
            mic._rlnAccumMotionTotal = Float(total)
            mic._rlnAccumMotionEarly = Float(early)
            mic._rlnAccumMotionLate = Float(late)

    def _saveAlignmentPlots(self, movieFName: str, pixSize: float):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movieFName)
        first, _ = self._getFramesRange()
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movieFName))
        plotter.close()

    def _getMovieShifts(self, movieFName: str):
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts are in pixels irrespective of any binning.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movieFName))
        xShifts, yShifts = parseMovieAlignment2(logPath)

        return xShifts, yShifts

    # def _moveOutput(self, movie):
    #     """ Move output from tmp to extra folder. """
    #
    #     def _moveToExtra(fn):
    #         """ Move file from movies tmp folder to extra """
    #         if exists(self._getCwdPath(movie, fn)):
    #             moveFile(self._getCwdPath(movie, fn),
    #                              self._getExtraPath(fn))
    #
    #     _moveToExtra(self._getMovieLogFile(movie))
    #     _moveToExtra(self._getMicFn(movie))
    #     _moveToExtra(replaceExt(movie.getBaseName(), "star"))
    #
    #     if self._doSaveUnweightedMic():
    #         _moveToExtra(self._getOutputMicName(movie))
    #
    #     if self.splitEvenOdd:
    #         _moveToExtra(self._getOutputMicEvenName(movie))
    #         _moveToExtra(self._getOutputMicOddName(movie))
    #
    #     if self.doSaveMovie:
    #         outputMicFn = self._getCwdPath(movie, self._getOutputMicName(movie))
    #         outputMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
    #         movieFn = outputMicFn.replace('_aligned_mic.mrc',
    #                                       '_aligned_mic_Stk.mrc')
    #         moveFile(movieFn, outputMovieFn)

    def _updateOutputSet(self, outputName, outputSet,
                         state=pwobj.Set.STREAM_OPEN):
        """ Redefine this method to set EER attrs. """
        first = getattr(self, '_firstUpdate', True)

        if first and outputName == 'outputMovies':
            og = OpticsGroups.fromImages(outputSet)
            gain = self.getInputMovies().getGain()
            ogDict = {'rlnMicrographStartFrame': self.alignFrame0.get()}
            if self.isEER:
                ogDict.update({'rlnEERGrouping': self.eerGroup.get(),
                               'rlnEERUpsampling': self.eerSampling.get() + 1})
            if gain:
                ogDict['rlnMicrographGainName'] = gain
            og.updateAll(**ogDict)
            og.toImages(outputSet)

        ProtAlignMovies._updateOutputSet(self, outputName, outputSet,
                                         state=state)
        self._firstUpdate = False

    def _createOutputMicrographs(self):
        createWeighted = self._createOutputWeightedMicrographs()
        # To create the unweighted average micrographs
        # we only consider the 'doSaveUnweightedMic' flag if the
        # weighted ones should be created.
        return not createWeighted or self.doSaveUnweightedMic

    def _doSaveUnweightedMic(self):
        """ Wraps the logic for saving unweighted mics that needs to consider the doApplyDoseFilter to be true"""
        return self._createOutputWeightedMicrographs() and self.doSaveUnweightedMic

    def _createOutputWeightedMicrographs(self):
        return self.doApplyDoseFilter

    def _doComputeMicThumbnail(self):
        return self.doComputeMicThumbnail

    def _useWorkerThread(self):
        return '--dont_use_worker_thread' not in self.extraProtocolParams.get()

    def getSamplingRate(self):
        return self.getInputMovies().getSamplingRate()

    def calcFrameMotion(self, movie):
        # based on relion 3.1 motioncorr_runner.cpp
        shiftsX, shiftsY = self._getMovieShifts(movie)
        a0, aN = self._getFramesRange()
        nframes = aN - a0 + 1
        preExp, dose = self._getCorrectedDose()
        # when using EER, the hardware frames are grouped
        if self.isEER:
            dose *= self.eerGroup.get()
        cutoff = (4 - preExp) // dose  # early is <= 4e/A^2
        total, early, late = 0., 0., 0.
        x, y, xOld, yOld = 0., 0., 0., 0.
        pix = self.getSamplingRate()
        try:
            for frame in range(2, nframes + 1):  # start from the 2nd frame
                x, y = shiftsX[frame - 1], shiftsY[frame - 1]
                d = sqrt((x - xOld) * (x - xOld) + (y - yOld) * (y - yOld))
                total += d
                if frame <= cutoff:
                    early += d
                else:
                    late += d
                xOld = x
                yOld = y
            return list(map(lambda x: pix * x, [total, early, late]))
        except IndexError:
            self.error(f"Expected {nframes} frames, found less. "
                       f"Check movie {movie.getFileName()}")

    def _getResultMicFn(self, movie: Movie):
        return self._getExtraPath(f'{basename(movie.getFileName())}_aligned_mic.mrc')


def createGlobalAlignmentPlot(meanX, meanY, first, pixSize):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []

    def px_to_ang(px):
        y1, y2 = px.get_ylim()
        x1, x2 = px.get_xlim()
        ax_ang2.set_ylim(y1 * pixSize, y2 * pixSize)
        ax_ang.set_xlim(x1 * pixSize, x2 * pixSize)
        ax_ang.figure.canvas.draw()
        ax_ang2.figure.canvas.draw()

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax_px = figure.add_subplot(111)
    ax_px.grid()
    ax_px.set_xlabel('Shift x (px)')
    ax_px.set_ylabel('Shift y (px)')

    ax_ang = ax_px.twiny()
    ax_ang.set_xlabel('Shift x (A)')
    ax_ang2 = ax_px.twinx()
    ax_ang2.set_ylabel('Shift y (A)')

    i = first
    # The output and log files list the shifts relative to the first frame.
    # ROB unit seems to be pixels since sampling rate is only asked
    # by the program if dose filtering is required
    skipLabels = ceil(len(meanX) / 10.0)
    labelTick = 1

    for x, y in zip(meanX, meanY):
        sumMeanX.append(x)
        sumMeanY.append(y)
        if labelTick == 1:
            ax_px.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    # automatically update lim of ax_ang when lim of ax_px changes.
    ax_px.callbacks.connect("ylim_changed", px_to_ang)
    ax_px.callbacks.connect("xlim_changed", px_to_ang)

    ax_px.plot(sumMeanX, sumMeanY, color='b')
    ax_px.plot(sumMeanX, sumMeanY, 'yo')
    ax_px.plot(sumMeanX[0], sumMeanY[0], 'ro', markersize=10, linewidth=0.5)
    ax_px.set_title('Global frame alignment')

    plotter.tightLayout()

    return plotter


