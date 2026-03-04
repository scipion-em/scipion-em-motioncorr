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
from os.path import exists, abspath, basename
from typing import Tuple, List, Union, Optional, Any
import numpy as np
import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.gui.plotter import Plotter
from pwem.objects import SetOfMovies, SetOfMicrographs, Movie, Micrograph, MovieAlignment, Image
from pyworkflow.object import Set, CsvList, Float, Pointer
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.utils import cyanStr, Message, redStr, removeBaseExt, getExt, weakImport
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from .. import Plugin
from .protocol_base import ProtMotionCorrBase
from ..convert import parseMovieAlignment2

logger = logging.getLogger(__name__)
MC_EVEN_ODD_ATTRIBUTE = '_mcEvenOddMics'
EVEN_SUFFIX = '_EVN'
ODD_SUFFIX = '_ODD'
DW_SUFFIX = '_DW'
STK_SUFFIX = '_Stk'


class MotionCorrOutputs(Enum):
    movies = SetOfMovies()
    micrographs = SetOfMicrographs()
    micrographsDW = SetOfMicrographs()
    micrographsEven = SetOfMicrographs()
    micrographsOdd = SetOfMicrographs()


class ProtMotionCorrNewStreaming(ProtMotionCorrBase, ProtStreamingBase):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)

    New Streaming refers to the next-generation engine developed by the Scipion Team.
    While both the new and legacy streaming protocols will coexist for a transitional
    period, they remain fully compatible with each other.
    """

    _label = 'movie alignment New Streaming'
    _possibleOutputs = MotionCorrOutputs
    stepsExecutionMode = cons.STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.itemIdReadList = []
        self.sRate = None
        self.gain = None
        self.dark = None
        self.failedMovies = []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported movies.')
        self._defineAlignmentParams(form)
        self._defineCommonParams(form)

    @staticmethod
    def _defineAlignmentParams(form):
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

        form.addParam('splitEvenOdd', params.BooleanParam,
                      default=False,
                      label='Split & sum odd/even frames?',
                      help='Generate odd and even sums using odd and even frames '
                           'respectively when this option is enabled.')

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Save aligned movie?")

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
        outputsToCheck = self._getOutputsToCheck()

        while True:
            with self._lock:
                inIds = set(inMoviesSet.getUniqueValues('id'))

            # In the if statement below, Counter is used because in the objId comparison the order doesn’t matter
            # but duplicates do. With a direct comparison, the closing step may not be inserted because of the order:
            # ['id_a', 'id_b'] != ['id_b', 'id_a'], but they are the same with Counter.
            if not inMoviesSet.isStreamOpen() and Counter(self.itemIdReadList) == Counter(inIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetStep,
                                         outputsToCheck,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedIds = inIds - set(self.itemIdReadList)
            moviesToProcessDict = {objId: movie.clone() for movie in inMoviesSet.iterItems()
                                   if (objId := movie.getObjId()) in nonProcessedIds}
            for objId, movie in moviesToProcessDict.items():
                movieFName = movie.getFileName()
                cInPId = self._insertFunctionStep(self.convertInputStep,
                                                  movieFName,
                                                  prerequisites=[],
                                                  needsGPU=False)
                pMovPid = self._insertFunctionStep(self.processMovieStep,
                                                   movieFName,
                                                   prerequisites=cInPId,
                                                   needsGPU=True)
                cOutId = self._insertFunctionStep(self.createOutputStep,
                                                  movieFName,
                                                  movie,
                                                  prerequisites=pMovPid,
                                                  needsGPU=False)
                closeSetStepDeps.append(cOutId)
                logger.info(cyanStr(f"Steps created for objId = {objId} - {movie.getFileName()}"))
                self.itemIdReadList.append(objId)

            time.sleep(10)
            if inMoviesSet.isStreamOpen():
                with self._lock:
                    inMoviesSet.loadAllProperties()  # refresh status for the streaming

    def _initialize(self):
        inputMovies = self.getInputMovies()
        self.gain = inputMovies.getGain()
        self.dark = inputMovies.getDark()
        self.sRate = inputMovies.getSamplingRate()
        self.isEER = getExt(inputMovies.getFirstItem().getFileName()) == ".eer"

    def convertInputStep(self, movieFName: str):
        try:
            # Convert the gain and dark images if they haven't been converted yet or if they weren't provided
            if self.dark or self.gain and not exists(self._getExtraPath('DONE')):
                super()._convertInputStep()
        except Exception as e:
            self.failedMovies.append(movieFName)
            logger.error(redStr(f"ERROR: movie convert failed for {movieFName} with the exception {e}"))
            traceback.print_exc()

    def processMovieStep(self, movieFName: str):
        if movieFName in self.failedMovies:
            return

        logger.info(cyanStr(f"Processing movie: {movieFName}"))
        outputMicFn = self._getResultMicFn(movieFName)
        argsDict = self._getMcArgs()
        argsDict['-OutMrc'] = f'{outputMicFn}'
        argsDict['-LogDir'] = f'{self._getExtraPath()}'
        args = self._getInputFormat(movieFName, absPath=True)
        args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        args += ' ' + self.extraParams2.get()

        try:
            self.runJob(Plugin.getProgram(), args, env=Plugin.getEnviron())
        except Exception as e:
            self.failedMovies.append(movieFName)
            logger.error(redStr(f"ERROR: Motioncor failed for {movieFName} with the exception {e}"))
            traceback.print_exc()
            return

        try:
            self._saveAlignmentPlots(movieFName, self.sRate)
        except Exception as e:
            self.failedMovies.append(movieFName)
            logger.error(redStr(f"ERROR: Saving the alignment plots failed for {movieFName} "
                                f"with the exception {e}"))
            traceback.print_exc()

    @retry_on_sqlite_lock(log=logger)
    def createOutputStep(self, movieFName: str, inMovie: Movie):
        if movieFName in self.failedMovies:
            return
        try:
            if self.doSaveMovie.get():
                outMovieFn = self._getResultMicFn(movieFName, suffix=STK_SUFFIX)
                setMRCSamplingRate(outMovieFn, self.sRate)
            else:
                outMovieFn = movieFName
            with self._lock:
                # SET OF MOVIES ----------------------------------------------------------------
                outputMovies = self._getOutputMovies()
                outMovie = Movie()
                outMovie.copyInfo(inMovie)
                outMovie.setFileName(outMovieFn)
                outMovie.setMicName(basename(outMovieFn))
                # Movie alignment
                n = outMovie.getNumberOfFrames()
                alignment = self.getMovieAlignment(movieFName, n)
                outMovie.setAlignment(alignment)
                # Data persistence
                outputMovies.append(outMovie)
                outputMovies.update(outMovie)
                outputMovies.write()
                self._store(outputMovies)

                # SET OF MICROGRAPHS ----------------------------------------------------------
                if self.doApplyDoseFilter.get():
                    suffix = DW_SUFFIX
                    outputName = self._possibleOutputs.micrographsDW.name
                else:
                    suffix = ''
                    outputName = self._possibleOutputs.micrographs.name
                self._registerMics(movieFName, inMovie, outputName, suffix=suffix)
                if self.splitEvenOdd.get():
                    # Even
                    outputName = self._possibleOutputs.micrographsEven.name
                    self._registerMics(movieFName, inMovie, outputName, suffix=EVEN_SUFFIX)
                    # Odd
                    outputName = self._possibleOutputs.micrographsOdd.name
                    self._registerMics(movieFName, inMovie, outputName, suffix=ODD_SUFFIX)

                # Close explicitly the outputs (for streaming)
                self.closeOutputsForStreaming()

        except Exception as e:
            logger.error(
                redStr(f'Movie = {movieFName} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    def closeOutputSetStep(self, attrib: Union[List[str], str]):
        self._closeOutputSet()
        attribList = [attrib] if type(attrib) is str else attrib
        failedOutputList = []
        for attr in attribList:
            outTsSet = getattr(self, attr, None)
            if not outTsSet or (outTsSet and len(outTsSet) == 0):
                failedOutputList.append(attr)
        if failedOutputList:
            raise Exception(f'No output/s {failedOutputList} were generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'outputMicrographs') or \
                hasattr(self, 'outputMicrographsDoseWeighted'):
            summary.append('Aligned %d movies using Motioncor.'
                           % self.getInputMovies().getSize())
            if self.splitEvenOdd.get() and self.doApplyDoseFilter.get():
                summary.append('Even/odd outputs are dose-weighted!')
        else:
            summary.append('Output is not ready')

        return summary

    def _validate(self):
        return ProtMotionCorrBase._validate(self)

    # --------------------------- UTILS functions -----------------------------
    def readingOutput(self) -> None:
        movies = getattr(self, self._possibleOutputs.movies.name, None)
        outputList = []
        mics = getattr(self, self._possibleOutputs.micrographs.name, None)
        micsDW = getattr(self, self._possibleOutputs.micrographsDW.name, None)
        micsEven = getattr(self, self._possibleOutputs.micrographsEven.name, None)
        micsOdd = getattr(self, self._possibleOutputs.micrographsOdd.name, None)
        if self.splitEvenOdd.get():
            outputList.append(micsEven)
            outputList.append(micsOdd)
        if self.doApplyDoseFilter.get():
            outputList.append(micsDW)
        else:
            outputList.append(mics)

        if movies and None not in outputList:
            for item in movies:
                self.itemIdReadList.append(item.getObjId())
            self.info(cyanStr(f'Ids processed: {self.itemIdReadList}'))
        else:
            self.info(cyanStr('No movies have been processed yet'))

    def getInputMovies(self, asPointer: bool = False) -> Union[Pointer, SetOfMovies]:
        return self.inputMovies if asPointer else self.inputMovies.get()

    def _getResultMicFn(self, movieFName: str, suffix: str = '') -> str:
        bName = removeBaseExt(movieFName).replace('.mrc', '')
        return self._getExtraPath(f'{bName}_aligned_mic{suffix}.mrc')

    def setMicPlotInfo(self, mic: Micrograph, movieFName: str) -> None:
        mic.plotGlobal = Image(location=self._getPlotGlobal(movieFName))
        if self.doApplyDoseFilter.get():
            total, early, late = self.calcFrameMotion(movieFName)
            mic._rlnAccumMotionTotal = Float(total)
            mic._rlnAccumMotionEarly = Float(early)
            mic._rlnAccumMotionLate = Float(late)

    def setMicsEvenOdd(self, inMovieFName: str, alignedMic: Micrograph) -> None:
        if self.splitEvenOdd.get():
            setattr(alignedMic, MC_EVEN_ODD_ATTRIBUTE, CsvList(pType=str))
            alignedMic._mcEvenOddMics.set([
                self._getResultMicFn(inMovieFName, suffix=ODD_SUFFIX),
                self._getResultMicFn(inMovieFName, suffix=EVEN_SUFFIX),
            ])

    def _getMovieLogFile(self, movieFName: str) -> str:
        usePatches = self.patchX != 0 or self.patchY != 0
        pattern = '-Patch' if usePatches else ''
        bName = removeBaseExt(movieFName).replace('.mrc', '')
        return abspath(self._getExtraPath(f'{bName}{pattern}-Full.log'))

    def _getPlotGlobal(self, movieFName: str) -> str:
        return self._getExtraPath(f'{removeBaseExt(movieFName)}_global_shifts.png')

    def _saveAlignmentPlots(self, movieFName: str, pixSize: float) -> None:
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movieFName)
        first, _ = self._getFramesRange()
        plotter = self.createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movieFName))
        plotter.close()

    def _getMovieShifts(self, movieFName: str) -> Tuple[List[float], List[float]]:
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts are in pixels irrespective of any binning.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movieFName))
        xShifts, yShifts = parseMovieAlignment2(logPath)

        return xShifts, yShifts

    def _registerMics(self,
                      movieFName: str,
                      inMovie: Movie, outputName: str,
                      suffix: str = '') -> None:
        outMicSet = self._getOutputMics(outputName, suffix=suffix)
        outMic = Micrograph()
        outMic.copyInfo(inMovie)
        micFn = self._getResultMicFn(movieFName, suffix=suffix)
        setMRCSamplingRate(micFn, self.sRate)
        outMic.setFileName(micFn)
        self.setMicPlotInfo(outMic, movieFName)
        if suffix in [DW_SUFFIX, '']:
            self.setMicsEvenOdd(movieFName, outMic)
        # Data persistence
        outMicSet.append(outMic)
        outMicSet.update(outMic)
        outMicSet.write()
        self._store(outMicSet)

    def _getOutputMovies(self) -> SetOfMovies:
        attrName = self._possibleOutputs.movies.name
        outputMovies = getattr(self, attrName, None)
        if outputMovies:
            outputMovies.enableAppend()
        else:
            inputMoviesPointer = self.getInputMovies(asPointer=True)
            outputMovies = SetOfMovies.create(self._getPath(), template='movies')
            outputMovies.copyInfo(self.getInputMovies())
            with weakImport("relion"):
                from relion.convert import OpticsGroups
                og = OpticsGroups.fromImages(outputMovies)
                gain = self.getInputMovies().getGain()
                ogDict = {'rlnMicrographStartFrame': self.alignFrame0.get()}
                if self.isEER:
                    ogDict.update({'rlnEERGrouping': self.eerGroup.get(),
                                   'rlnEERUpsampling': self.eerSampling.get() + 1})
                if gain:
                    ogDict['rlnMicrographGainName'] = gain
                og.updateAll(**ogDict)
                og.toImages(outputMovies)
            outputMovies.setStreamState(Set.STREAM_OPEN)
            outputMovies.write()  # Write set properties, otherwise it may expose the set (sqlite) without properties.

            self._defineOutputs(**{attrName: outputMovies})
            self._defineSourceRelation(inputMoviesPointer, outputMovies)
        return outputMovies

    def getMovieAlignment(self, inMovieFName: str, nFrames: int) -> MovieAlignment:
        first, last = self._getFrameRange(nFrames, 'align')

        if self.doSaveMovie.get():  # Interpolated movie
            xShifts = np.zeros(nFrames).tolist()
            yShifts = xShifts
        else:
            logFile = self._getMovieLogFile(inMovieFName)
            xShifts, yShifts = parseMovieAlignment2(logFile)

        alignment = MovieAlignment(first=first,
                                   last=last,
                                   xshifts=xShifts,
                                   yshifts=yShifts)
        roiList = [self.getAttributeValue(s, 0) for s in
                   ['cropOffsetX', 'cropOffsetY', 'cropDimX', 'cropDimY']]
        alignment.setRoi(roiList)
        return alignment

    def _getOutputMics(self, outputName: str, suffix: str = '') -> SetOfMicrographs:
        outputMics = getattr(self, outputName, None)
        if outputMics:
            outputMics.enableAppend()
        else:
            inputMoviesPointer = self.getInputMovies(asPointer=True)
            outputMics = SetOfMicrographs.create(self._getPath(), template='movies', suffix=suffix)
            outputMics.copyInfo(self.getInputMovies())
            outputMics.setStreamState(Set.STREAM_OPEN)
            outputMics.write()  # Write set properties, otherwise it may expose the set (sqlite) without properties.

            self._defineOutputs(**{outputName: outputMics})
            self._defineSourceRelation(inputMoviesPointer, outputMics)
        return outputMics

    def _getOutputsToCheck(self) -> List[str]:
        outputsToCheck = [
            self._possibleOutputs.movies.name,
        ]
        if self.doApplyDoseFilter.get():
            outputsToCheck.append(self._possibleOutputs.micrographsDW.name)
        else:
            outputsToCheck.append(self._possibleOutputs.micrographs.name)
        if self.splitEvenOdd.get():
            outputsToCheck.append(self._possibleOutputs.micrographsEven.name)
            outputsToCheck.append(self._possibleOutputs.micrographsOdd.name)
        return outputsToCheck

    def closeOutputsForStreaming(self):
        # Close explicitly the outputs (for streaming)
        for outputName in self._possibleOutputs:
            output = getattr(self, outputName.name, None)
            if output:
                output.close()

    def _getFrameRange(self, n: int, prefix: str) -> Tuple[int, int]:
        """
        Params:
        :param n: Number of frames of the movies
        :param prefix: what range we want to consider, either 'align' or 'sum'
        :return: (i, f) initial and last frame range
        """
        # In case that the user select the same range for ALIGN and SUM
        # we also use the 'align' prefix
        if self.useAlignToSum.get():
            prefix = 'align'

        first = self.getAttributeValue('%sFrame0' % prefix)
        last = self.getAttributeValue('%sFrameN' % prefix)

        if first <= 1:
            first = 1

        if last <= 0:
            last = n

        return first, last

    def calcFrameMotion(self, movieFName: str) -> Optional[List[Any]]:
        # based on relion 3.1 motioncorr_runner.cpp
        shiftsX, shiftsY = self._getMovieShifts(movieFName)
        a0, aN = self._getFramesRange()
        nframes = aN - a0 + 1
        preExp, dose = self._getCorrectedDose()
        # when using EER, the hardware frames are grouped
        if self.isEER:
            dose *= self.eerGroup.get()
        cutoff = (4 - preExp) // dose  # early is <= 4e/A^2
        total, early, late = 0., 0., 0.
        x, y, xOld, yOld = 0., 0., 0., 0.
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
            return list(map(lambda x: self.sRate * x, [total, early, late]))
        except IndexError:
            logger.error(redStr(f"Expected {nframes} frames, found less. "
                                f"Check movie {movieFName}"))

    @staticmethod
    def createGlobalAlignmentPlot(meanX: List[float],
                                  meanY: List[float],
                                  first: int,
                                  pixSize: float) -> Plotter:
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
