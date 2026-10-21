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

import os
import time
from math import ceil, sqrt
from os.path import exists
from threading import Thread

import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.plotter import Plotter
from pwem.objects import Image, Float, Movie
from pwem.protocols import ProtAlignMovies
from pyworkflow.utils import weakImport

from .. import Plugin, DEFECTS_FILE_EER
from .protocol_base import ProtMotionCorrBase

# Each movie's step returns one file here as its result file. On Continue,
# pyworkflow keeps a finished step finished only while its result files
# exist, so removing the ones of failed movies makes their steps run again.
STEP_RESULTS_DIR = 'STEP_RESULTS'
# Image extensions a movie root must not end with, see _getMovieRoot.
IMAGE_EXTENSIONS = ('.mrc', '.mrcs', '.tif', '.tiff', '.eer')
# Every output a movie is published to, when the run produces it.
OUTPUT_NAMES = ('outputMovies', 'outputMicrographs',
                'outputMicrographsDoseWeighted')


class ProtMotionCorr(ProtMotionCorrBase, ProtAlignMovies):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'movie alignment'
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        ProtMotionCorrBase.__init__(self, **kwargs)

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif', '.eer', '.gain'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
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

        self._defineCommonParams(form)

    # --------------------------- STEPS functions -----------------------------
    def _convertInputStep(self):
        # Legacy ProtAlignMovies.processMovieStep still writes per-movie
        # DONE markers under extra/DONE. Keep that compatibility local
        # to the classic protocol only; streaming protocols must not use
        # filesystem checkpoints for durable progress.
        pwutils.makePath(self._getExtraPath('DONE'))
        ProtMotionCorrBase._convertInputStep(self)

    def _processMovie(self, movie):
        inputMovies = self.getInputMovies()

        try:
            # processMovieStep (the pwem base class step calling this
            # hook) has no exception boundary of its own around it - a
            # single movie whose format/args fail to build here (e.g.
            # _getInputFormat raising ValueError for an unrecognized
            # extension) must not crash the whole protocol, so this
            # setup is inside the try like the runJob failure already
            # was.
            movieFolder = self._getOutputMovieFolder(movie)
            outputMicFn = self._getOutputMicName(movie)

            argsDict = self._getMcArgs()
            # runJob uses cwd=movieFolder; resolve a gain relative to the
            # project directory before entering this per-movie directory.
            gainArg = argsDict.get('-Gain')
            if gainArg:
                gainPath = str(gainArg).strip().strip('"')
                if not os.path.isabs(gainPath):
                    argsDict['-Gain'] = f'"{os.path.abspath(gainPath)}"'
            argsDict['-OutMrc'] = f'"{outputMicFn}"'
            argsDict['-LogDir'] = './'
            if self.isEER:
                argsDict.update({'-FmIntFile': "../../extra/FmIntFile.txt"})
            if self.defectFile.get():
                argsDict['-DefectFile'] = self.defectFile.get()
            elif self.defectMap.get():
                argsDict['-DefectMap'] = self.defectMap.get()
            elif exists(self._getExtraPath(DEFECTS_FILE_EER)):
                argsDict['-DefectFile'] = f"../../extra/{DEFECTS_FILE_EER}"
            args = self._getInputFormat(movie.getFileName())
            args += ' '.join(['%s %s' % (k, v)
                              for k, v in argsDict.items()])
            args += ' ' + self.extraParams2.get()

            self.runJob(Plugin.getProgram(), args, cwd=movieFolder,
                        env=Plugin.getEnviron())
            self._moveOutput(movie)

            def _extraWork():
                # we need to move shifts log to extra dir before parsing
                try:
                    self._saveAlignmentPlots(movie, inputMovies.getSamplingRate())
                    outMicFn = self._getExtraPath(self._getMicFn(movie))

                    if self.doComputePSD:
                        self._computePSD(outMicFn, outputFn=self._getPsdCorr(movie))

                    if self._doComputeMicThumbnail():
                        self.computeThumbnail(
                            outMicFn, outputFn=self._getOutputMicThumbnail(movie))
                except:
                    self.error(f"ERROR: Extra work (i.e plots, PSD, thumbnail) "
                               f"has failed for {movie.getFileName()}\n")

            if self._useWorkerThread():
                thread = Thread(target=_extraWork)
                thread.start()
            else:
                _extraWork()

        except Exception as e:
            self.error(f"ERROR: Motioncor has failed for {movie.getFileName()} --> {str(e)}\n")
            import traceback
            traceback.print_exc()
        
    def _insertFinalSteps(self, deps):
        stepsId = []
        if self._useWorkerThread():
            stepsId.append(self._insertFunctionStep('waitForThreadStep',
                                                    prerequisites=deps,
                                                    needsGPU=False))
        return stepsId

    # ------------------------- failed movies ----------------------------
    def _insertAllSteps(self):
        self._prepareFailedForRetry()
        ProtAlignMovies._insertAllSteps(self)

    def processMovieStep(self, movieDict, hasAlignment):
        ProtAlignMovies.processMovieStep(self, movieDict, hasAlignment)

        movie = Movie()
        movie.setAttributesFromDict(movieDict, setBasic=True,
                                    ignoreMissing=True)

        return self._stepResultFiles([movie.getObjId()])

    def _stepResultPath(self, movieId):
        return self._getExtraPath(STEP_RESULTS_DIR, '%06d.TXT' % movieId)

    def _stepResultFiles(self, movieIds):
        """Result files of a movie's step, whatever became of the movie:
        a failed one has its file removed when the next run starts."""
        pwutils.makePath(self._getExtraPath(STEP_RESULTS_DIR))
        paths = [self._stepResultPath(movieId) for movieId in movieIds]

        for path in paths:
            open(path, 'a').close()

        return paths

    def _publishedMovieIds(self):
        """Movies in every output there is: one only in some still lacks
        the others."""
        published = None

        for outputName in OUTPUT_NAMES:
            outputSet = getattr(self, outputName, None)

            if outputSet is not None:
                ids = set(outputSet.getIdSet())
                published = ids if published is None else published & ids

        return published or set()

    def _prepareFailedForRetry(self):
        """Have Continue align again the movies that failed.

        pwem counts a movie done by its DONE marker and its line in
        DONE/all.TXT, and pyworkflow does not run a finished step again:
        all three go for the movies that failed and are not published.
        """
        failedIds = set(self._readFailedList()) - self._publishedMovieIds()

        if not failedIds:
            return

        for movieId in failedIds:
            pwutils.cleanPath(self._stepResultPath(movieId))
            movie = Movie()
            movie.setObjId(movieId)
            pwutils.cleanPath(self._getMovieDone(movie))

        doneIds = [movieId for movieId in self._readDoneList()
                   if movieId not in failedIds]

        with open(self._getAllDone(), 'w') as handle:
            handle.writelines('%d\n' % movieId for movieId in doneIds)

    def _movieResultExists(self, movie):
        """Whether every micrograph this movie's outputs are built from
        exists: a movie is published to all of them or to none."""
        names = []

        if self._createOutputMicrographs():
            names.append(self._getOutputMicName(movie))
        if self._createOutputWeightedMicrographs():
            names.append(self._getOutputMicWtName(movie))
        if self._doSplitEvenOdd():
            names.extend([self._getOutputMicEvenName(movie),
                          self._getOutputMicOddName(movie)])

        return all(os.path.exists(self._getExtraPath(name)) for name in names)

    def _updateOutputSets(self, newDone, streamMode):
        failed = [movie for movie in newDone
                  if not self._movieResultExists(movie)]

        if failed:
            self.error("Movie(s) %s gave no aligned micrograph and are left "
                       "out; Continue aligns them again."
                       % [movie.getObjId() for movie in failed])
            self._writeFailedList(failed)

        self.__dict__.setdefault('_handledMovieIds', set()).update(
            movie.getObjId() for movie in newDone)
        published = [movie for movie in newDone if movie not in failed]

        if not published and streamMode != pwobj.Set.STREAM_CLOSED:
            return

        # pwem takes the first movie it is given for the output summary:
        # only while there is a movie to give, and no output yet.
        self._firstTimeOutput = (bool(published)
                                 and getattr(self, 'outputMovies', None) is None)
        ProtAlignMovies._updateOutputSets(self, published, streamMode)

    def _checkNewOutput(self):
        ProtAlignMovies._checkNewOutput(self)

        if getattr(self, 'finished', False):
            self._failIfEveryMovieFailed()

    def _failIfEveryMovieFailed(self):
        """A failed movie is left out; the run fails only when that
        happened to every movie."""
        if getattr(self, '_handledMovieIds', None) and not any(
                getattr(self, outputName, None) is not None
                and getattr(self, outputName).getSize()
                for outputName in OUTPUT_NAMES):
            raise Exception("Motioncor failed for every one of the %d "
                            "movies; see the errors above. Continue aligns "
                            "them again." % len(self._handledMovieIds))

    def waitForThreadStep(self):
        # Quick and dirty (maybe desperate) way to wait
        # if the PSD and thumbnail were computed in a thread
        # If running in streaming this will not be necessary
        time.sleep(10)  # wait 10 sec for the thread to finish

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
    def _getCwdPath(self, movie, path):
        return os.path.join(self._getOutputMovieFolder(movie), path)

    def _getMovieRoot(self, movie):
        """Preserve dots in movie names, especially acquisition timestamps.

        The pwem implementation assumes another short suffix is a second
        extension, truncating names such as May14_21.14.32_1.tif to May14_21.14.
        We only strip a real compression suffix when present.
        """
        root = self._getMovieFileRoot(movie)
        # A movie named like movie.mrc.eer keeps an image extension in its
        # root, and MotionCor cuts the output name it is given at it: asked
        # for movie.mrc_aligned_mic.mrc it writes movie.mrc, which is never
        # found again.
        root, extension = os.path.splitext(root)
        if extension.lower() not in IMAGE_EXTENSIONS:
            root += extension
        return root

    def _getMovieFileRoot(self, movie):
        """The movie's file name without its extension (or compression)."""
        filename = movie.getAttributeValue('_originalFileName',
                                           movie.getFileName())
        root, extension = os.path.splitext(os.path.basename(filename))
        if extension.lower() in ('.gz', '.bz2', '.zip'):
            root, _ = os.path.splitext(root)
        return root

    def _getMovieLogFile(self, movie):
        # MotionCor3 bases the log on the actual -InTiff/-InMrc/-InEer
        # filename, not on _originalFileName (which may differ).
        movieRoot = os.path.splitext(os.path.basename(movie.getFileName()))[0]
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s%s-Full.log' % (movieRoot,
                                  '-Patch' if usePatches else '')

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

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
        if self.doApplyDoseFilter and self._hasValidDose():
            total, early, late = self.calcFrameMotion(movie)
            mic._rlnAccumMotionTotal = Float(total)
            mic._rlnAccumMotionEarly = Float(early)
            mic._rlnAccumMotionLate = Float(late)

    def _saveAlignmentPlots(self, movie, pixSize):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFramesRange()
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _moveOutput(self, movie):
        """ Move output from tmp to extra folder. """
        def _moveToExtra(fn):
            """ Move file from movies tmp folder to extra """
            if os.path.exists(self._getCwdPath(movie, fn)):
                pwutils.moveFile(self._getCwdPath(movie, fn),
                                 self._getExtraPath(fn))

        _moveToExtra(self._getMovieLogFile(movie))
        _moveToExtra(self._getMicFn(movie))
        _moveToExtra(pwutils.replaceExt(movie.getBaseName(), "star"))

        if self._doSaveUnweightedMic():
            _moveToExtra(self._getOutputMicName(movie))

        if self.splitEvenOdd:
            _moveToExtra(self._getOutputMicEvenName(movie))
            _moveToExtra(self._getOutputMicOddName(movie))

        if self.doSaveMovie:
            outputMicFn = self._getCwdPath(movie, self._getOutputMicName(movie))
            outputMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
            movieFn = outputMicFn.replace('_aligned_mic.mrc',
                                          '_aligned_mic_Stk.mrc')
            pwutils.moveFile(movieFn, outputMovieFn)

    def _updateOutputSet(self, outputName, outputSet,
                         state=pwobj.Set.STREAM_OPEN):
        """ Redefine this method to set EER attrs. """
        first = getattr(self, '_firstUpdate', True)

        if first and outputName == 'outputMovies':
            with weakImport('relion'):
                from relion.convert.convert31 import OpticsGroups
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
