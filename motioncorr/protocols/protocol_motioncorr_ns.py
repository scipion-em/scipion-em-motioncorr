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
    moviesFailed = SetOfMovies()


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
        # ProtStreamingBase keeps the generator step running while
        # processing steps consume additional executor threads.
        form.getParam('numberOfThreads').setDefault(2)

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
        self.readingOutput()
        outputsToCheck = self._getOutputsToCheck()

        prepareInputStepId = self._insertFunctionStep(
            self._convertInputStep,
            prerequisites=[],
            needsGPU=False,
        )

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
                pMovPid = self._insertFunctionStep(
                    self.processMovieStep,
                    movieFName,
                    prerequisites=prepareInputStepId,
                    needsGPU=True,
                )
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
        self.sRate = inputMovies.getSamplingRate() * self.binFactor.get()

        # This runs once, unconditionally, before the streaming loop
        # even starts. In live-acquisition streaming the input Set can
        # still be empty at launch (a normal race with the upstream
        # import-movies protocol) - getFirstItem() returns None for an
        # empty Set, and .getFileName() on None would crash the whole
        # protocol at t=0 instead of simply waiting for the first movie
        # like the rest of this generator already does.
        firstItem = inputMovies.getFirstItem()
        while firstItem is None and inputMovies.isStreamOpen():
            time.sleep(10)
            with self._lock:
                inputMovies.loadAllProperties()
            firstItem = inputMovies.getFirstItem()

        if firstItem is None:
            raise RuntimeError(
                "No movies available and the input stream is already "
                "closed; nothing to process."
            )

        self.isEER = getExt(firstItem.getFileName()) == ".eer"

        if (self.doApplyDoseFilter.get()
                and not self._hasValidDose()):
            raise RuntimeError(
                "Input movies do not contain usable dose information "
                "for dose weighting."
            )
    def processMovieStep(self, movieFName: str):
        if movieFName in self.failedMovies:
            return

        logger.info(cyanStr(f"Processing movie: {movieFName}"))

        try:
            # createOutputStep/the streaming generator has no exception
            # boundary of its own around this step, so building the
            # arguments (which can raise, e.g. _getInputFormat for an
            # unsupported extension) must be inside the try like the
            # runJob failure already is - otherwise a single movie
            # crashes the whole protocol instead of being routed to
            # failedMovies like every other failure in this function.
            outputMicFn = self._getResultMicFn(movieFName)
            argsDict = self._getMcArgs()
            argsDict['-OutMrc'] = f'{outputMicFn}'
            argsDict['-LogDir'] = f'{self._getExtraPath()}'
            args = self._getInputFormat(movieFName, absPath=True)
            args += ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
            args += ' ' + self.extraParams2.get()

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

    def createOutputStep(self, movieFName: str, inMovie: Movie):
        if movieFName in self.failedMovies:
            return

        if self.doSaveMovie.get():
            outMovieFn = self._getResultMicFn(movieFName, suffix=STK_SUFFIX)
            try:
                # Only the data-building call is isolated here.
                # Persistence failures below (append/update/write/_store)
                # must keep propagating instead of being converted into
                # per-movie processing failures.
                setMRCSamplingRate(outMovieFn, self.sRate)
            except Exception as e:
                self.failedMovies.append(movieFName)
                logger.error(redStr(f"ERROR: Patching the movie stack header failed for "
                                    f"{movieFName} with the exception {e}"))
                traceback.print_exc()
                return
        else:
            outMovieFn = movieFName

        with self._lock:
                # SET OF MOVIES ----------------------------------------------------------------
                outputMovies = self._getOutputMovies()
                if (len(outputMovies) == 0 or
                        inMovie.getObjId() not in outputMovies):
                    outMovie = Movie()
                    outMovie.copyInfo(inMovie)
                    outMovie.copyObjId(inMovie)
                    outMovie.setFileName(outMovieFn)
                    outMovie.setMicName(basename(outMovieFn))
                    # Movie alignment
                    n = outMovie.getNumberOfFrames()
                    try:
                        # Same isolation rationale as setMRCSamplingRate
                        # above - a movie whose alignment log fails to
                        # parse must not crash the whole protocol, but
                        # the persistence calls right after this must
                        # stay unprotected.
                        alignment = self.getMovieAlignment(movieFName, n)
                    except Exception as e:
                        self.failedMovies.append(movieFName)
                        logger.error(redStr(f"ERROR: Parsing the alignment for {movieFName} "
                                            f"failed with the exception {e}"))
                        traceback.print_exc()
                        return
                    outMovie.setAlignment(alignment)
                    # Data persistence
                    firstOutputMovie = outputMovies.getSize() == 0
                    outputMovies.append(outMovie)
                    if firstOutputMovie:
                        self._updateOutputMoviesOptics(outputMovies)
                    outputMovies.update(outMovie)
                    outputMovies.write()
                    self._store(outputMovies)

                # SET OF MICROGRAPHS ----------------------------------------------------------
                try:
                    # Match _getMcArgs: MotionCor2 was only asked to
                    # dose-weight (and therefore only writes a _DW
                    # output) when a usable dose was actually
                    # available - requesting the DW output here
                    # whenever doApplyDoseFilter is merely checked on
                    # the form would look for a file that was never
                    # produced.
                    if ProtMotionCorrBase._useDoseWeightedOutput(self):
                        suffix = DW_SUFFIX
                        outputName = self._possibleOutputs.micrographsDW.name
                    else:
                        suffix = ''
                        outputName = self._possibleOutputs.micrographs.name
                    if not self._registerMics(movieFName, inMovie, outputName, suffix=suffix):
                        return
                    if self.splitEvenOdd.get():
                        # Even
                        outputName = self._possibleOutputs.micrographsEven.name
                        if not self._registerMics(movieFName, inMovie, outputName, suffix=EVEN_SUFFIX):
                            return
                        # Odd
                        outputName = self._possibleOutputs.micrographsOdd.name
                        if not self._registerMics(movieFName, inMovie, outputName, suffix=ODD_SUFFIX):
                            return

                    # Close explicitly the outputs (for streaming)
                    self.closeOutputsForStreaming()

                except Exception as e:
                    logger.error(
                        redStr(f'Movie = {movieFName} -> Unable to register the output with exception {e}.'))
                    logger.error(traceback.format_exc())
                    raise

    def closeOutputSetStep(self, attrib: Union[List[str], str]):
        attribList = [attrib] if type(attrib) is str else attrib
        failedOutputList = []
        for attr in attribList:
            outputSet = getattr(self, attr, None)
            if not outputSet or len(outputSet) == 0:
                failedOutputList.append(attr)

        if self.failedMovies:
            if len(failedOutputList) == len(attribList):
                raise RuntimeError(
                    f"Motioncor failed for all {len(self.failedMovies)} movie(s). "
                    "No useful output was generated."
                )

            self._createFailedMoviesOutput()
            self.warning(
                f"Motioncor failed for {len(self.failedMovies)} movie(s). "
                f"They are available in '{self._possibleOutputs.moviesFailed.name}'."
            )

        if failedOutputList:
            raise RuntimeError(
                f"No output/s {failedOutputList} were generated. "
                "Please check the Output Log > run.stdout and run.stderr"
            )

        self._closeOutputSet()

    def _createFailedMoviesOutput(self):
        """Publish locally failed movies as a regular Scipion output Set."""
        inputMovies = self.getInputMovies()
        failedNames = set(self.failedMovies)
        outputName = self._possibleOutputs.moviesFailed.name

        failedMovies = self._createSetOfMovies(suffix='failed')
        failedMovies.copyInfo(inputMovies)

        for movie in inputMovies.iterItems():
            if movie.getFileName() in failedNames:
                failedMovies.append(movie.clone())

        failedMovies.setStreamState(Set.STREAM_CLOSED)
        failedMovies.write()
        self._defineOutputs(**{outputName: failedMovies})
        self._defineSourceRelation(self.getInputMovies(asPointer=True), failedMovies)

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
        errors = ProtMotionCorrBase._validate(self)
        self._validateThreads(errors)
        return errors

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
        if ProtMotionCorrBase._useDoseWeightedOutput(self):
            outputList.append(micsDW)
        else:
            outputList.append(mics)

        if movies and None not in outputList:
            completedIds = set(item.getObjId() for item in movies)
            for outputSet in outputList:
                completedIds.intersection_update(
                    item.getObjId() for item in outputSet
                )

            completedIds.intersection_update(
                self.getInputMovies().getUniqueValues('id')
            )

            for item in movies:
                if item.getObjId() in completedIds:
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
        if ProtMotionCorrBase._useDoseWeightedOutput(self):
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
                      suffix: str = '') -> bool:
        outMicSet = self._getOutputMics(outputName, suffix=suffix)
        if (len(outMicSet) > 0 and
                inMovie.getObjId() in outMicSet):
            return True

        outMic = Micrograph()
        outMic.copyInfo(inMovie)
        outMic.copyObjId(inMovie)
        micFn = self._getResultMicFn(movieFName, suffix=suffix)
        try:
            # Isolate data-building (reading the mic's own header,
            # parsing its alignment log for frame-motion stats) from
            # the persistence calls below. The external tool can
            # finish successfully without producing every expected
            # output file for a given movie (e.g. too few frames for
            # dose weighting) - a missing/unreadable file here must
            # not crash the whole protocol, but persistence stays
            # unprotected like everywhere else in this class.
            setMRCSamplingRate(micFn, self.sRate)
            outMic.setFileName(micFn)
            outMic.setSamplingRate(self.sRate)
            self.setMicPlotInfo(outMic, movieFName)
        except Exception as e:
            self.failedMovies.append(movieFName)
            logger.error(redStr(f"ERROR: Registering micrograph output failed for "
                                f"{movieFName} with the exception {e}"))
            traceback.print_exc()
            return False

        if suffix in [DW_SUFFIX, '']:
            self.setMicsEvenOdd(movieFName, outMic)
        # Data persistence
        outMicSet.append(outMic)
        outMicSet.update(outMic)
        outMicSet.write()
        self._store(outMicSet)
        return True

    def _getOutputMovies(self) -> SetOfMovies:
        attrName = self._possibleOutputs.movies.name
        outputMovies = getattr(self, attrName, None)
        if outputMovies is not None:
            outputMovies.enableAppend()
        else:
            inputMoviesPointer = self.getInputMovies(asPointer=True)
            outputMovies = SetOfMovies.create(self._getPath(), template='movies')
            outputMovies.copyInfo(self.getInputMovies())
            outputMovies.setSamplingRate(self.sRate)
            outputMovies.setStreamState(Set.STREAM_OPEN)
            outputMovies.write()  # Persist set properties before exposing the streaming output.

            self._defineOutputs(**{attrName: outputMovies})
            self._defineSourceRelation(inputMoviesPointer, outputMovies)
        return outputMovies

    def _updateOutputMoviesOptics(
            self,
            outputMovies: SetOfMovies,
    ) -> None:
        """Initialize Relion optics after the first output movie exists.

        The classic MotionCorr protocol builds optics only after the first
        output movie has populated the Set. Doing it while the streaming Set
        is still empty can leave required optics values as None.
        """
        with weakImport("relion"):
            from relion.convert import OpticsGroups

            og = OpticsGroups.fromImages(outputMovies)
            gain = self.getInputMovies().getGain()
            ogDict = {
                'rlnMicrographStartFrame': self.alignFrame0.get()
            }

            if self.isEER:
                ogDict.update({
                    'rlnEERGrouping': self.eerGroup.get(),
                    'rlnEERUpsampling': self.eerSampling.get() + 1,
                })

            if gain:
                ogDict['rlnMicrographGainName'] = gain

            og.updateAll(**ogDict)
            og.toImages(outputMovies)


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
        if outputMics is not None:
            outputMics.enableAppend()
        else:
            inputMoviesPointer = self.getInputMovies(asPointer=True)
            outputMics = SetOfMicrographs.create(self._getPath(), template='movies', suffix=suffix)
            outputMics.copyInfo(self.getInputMovies())
            outputMics.setSamplingRate(self.sRate)
            outputMics.setStreamState(Set.STREAM_OPEN)
            outputMics.write()  # Persist set properties before exposing the streaming output.

            self._defineOutputs(**{outputName: outputMics})
            self._defineSourceRelation(inputMoviesPointer, outputMics)
        return outputMics

    def _getOutputsToCheck(self) -> List[str]:
        outputsToCheck = [
            self._possibleOutputs.movies.name,
        ]
        if ProtMotionCorrBase._useDoseWeightedOutput(self):
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
            if output is not None:
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
        # dose can legitimately be 0.0 when the acquisition's dose per
        # frame is missing (_getCorrectedDose degrades to 0.0 instead
        # of crashing, and this is now a non-blocking _warnings()
        # notice rather than a hard _validate() error) - without a
        # known dose there is no way to tell when 4 e/A^2 was reached,
        # so treat every frame as "early" instead of raising
        # ZeroDivisionError.
        cutoff = (4 - preExp) // dose if dose else nframes  # early is <= 4e/A^2
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
