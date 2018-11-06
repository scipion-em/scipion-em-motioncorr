# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
from itertools import izip
from math import ceil
from threading import Thread

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.data import MovieAlignment
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.gui.plotter import Plotter

import motioncorr
from motioncorr.convert import *
from motioncorr.constants import *


class ProtMotionCorr(ProtAlignMovies):
    """ This protocol wraps motioncor2 movie alignment program developed at UCSF.

    Motioncor2 performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'movie alignment'
    CONVERT_TO_MRC = 'mrc'

    def __init__(self, **args):
        ProtAlignMovies.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def versionGE(self, version):
        """ Return True if current version of motioncor2 is greater
         or equal than the input argument.
         Params:
            version: string version (semantic version, e.g 1.0.1)
        """
        v1 = int(self._getVersion().replace('.', ''))
        v2 = int(version.replace('.', ''))

        if v1 < v2:
                return False

        return True

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=cons.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " Motioncor2 can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addHidden('doSaveAveMic', params.BooleanParam,
                       default=True)
        form.addHidden('useAlignToSum', params.BooleanParam,
                       default=True)

        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN and SUM',
                             help='Frames range to ALIGN and SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie.')
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

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save movie",
                      help="Save Aligned movie")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      label="Compute PSD (before/after)?",
                      help="If Yes, the protocol will compute for each movie "
                           "the average PSD before and after alignment, "
                           "for comparison")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail and keep it with the '
                           'micrograph object for visualization purposes. ')

        #form.addParam('computeAllFramesAvg', params.BooleanParam,
        #              expertLevel=cons.LEVEL_ADVANCED,
        #              default=False,
        #              label='Compute all frames average?',
        #              help='Computing all the frames average could provide a '
        #                   'sanity check about the microscope and the camera.')

        form.addParam('extraProtocolParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional protocol parameters',
                      help="Here you can provide some extra parameters for the "
                           "protocol, not the underlying motioncor2 program."
                           "You can provide many options separated by space. "
                           "\n\n*Options:* \n\n"
                           "--use_worker_thread \n"
                           " Use an extra thread to compute"
                           " PSD and thumbnail. This will allow a more effective"
                           " use of the GPU card, but requires an extra CPU. ")

        form.addSection(label="Motioncor2 params")
        form.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                      label='Apply dose filter',
                      help='Apply a dose-dependent filter to frames before '
                           'summing them. Pre-exposure and dose per frame '
                           'should  be specified during movies import.')

        line = form.addLine('Number of patches',
                            help='Number of patches to be used for patch based '
                                 'alignment. Set to *0 0* to do only global motion '
                                 'correction. \n')
        line.addParam('patchX', params.IntParam, default=5, label='X')
        line.addParam('patchY', params.IntParam, default=5, label='Y')

        if self.versionGE('1.0.1'): # Patch overlap was introduced in 1.0.1
            form.addParam('patchOverlap', params.IntParam, default=0,
                          label='Patches overlap (%)',
                          help='In versions > 1.0.1 it is possible to specify'
                               'the overlapping between patches. '
                               '\nFor example, overlap=20 means that '
                               'each patch will have a 20% overlapping \n'
                               'with its neighboring patches in each dimension.')

        form.addParam('group', params.IntParam, default='1',
                      label='Group N frames',
                      help='Group every specified number of frames by adding '
                           'them together. The alignment is then performed on '
                           'the summed frames. By default, no grouping is '
                           'performed.')

        form.addParam('tol', params.FloatParam, default='0.5',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Tolerance (px)',
                      help='Tolerance for iterative alignment, default *0.5px*.')

        form.addParam('doSaveUnweightedMic', params.BooleanParam, default=True,
                      condition='doSaveAveMic and doApplyDoseFilter',
                      label="Save unweighted micrographs?",
                      help="Yes by default, if you have selected to apply a "
                           "dose-dependent filter to the frames")

        group = form.addGroup('Magnification correction')
        group.addParam('doMagCor', params.BooleanParam, default=False,
                       label='Correct anisotropic magnification?',
                       help='Correct anisotropic magnification by '
                            'stretching image along the major axis, '
                            'the axis where the lower magnification is '
                            'detected.')
        group.addParam('useEst', params.BooleanParam, default=True,
                       label='Use previous estimation?',
                       condition='doMagCor',
                       help='Use previously calculated parameters of '
                            'magnification anisotropy (from magnification '
                            'distortion estimation protocol).')
        group.addParam('inputEst', params.PointerParam,
                       pointerClass='ProtMagDistEst',
                       condition='useEst and doMagCor',
                       label='Input protocol',
                       help='Select previously executed estimation protocol.')
        group.addParam('scaleMaj', params.FloatParam, default=1.0,
                       condition='not useEst and doMagCor',
                       label='Major scale factor',
                       help='Major scale factor.')
        group.addParam('scaleMin', params.FloatParam, default=1.0,
                       condition='not useEst and doMagCor',
                       label='Minor scale factor',
                       help='Minor scale factor.')
        group.addParam('angDist', params.FloatParam, default=0.0,
                       condition='not useEst and doMagCor',
                       label='Distortion angle (deg)',
                       help='Distortion angle, in degrees.')

        form.addParam('gainRot', params.EnumParam,
                      choices=['no rotation', '90 degrees',
                               '180 degrees', '270 degrees'],
                      label="Rotate gain reference:",
                      default=NO_ROTATION,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help="Rotate gain reference counter-clockwise.")

        form.addParam('gainFlip', params.EnumParam,
                      expertLevel=cons.LEVEL_ADVANCED,
                      choices=['no flip', 'upside down', 'left right'],
                      label="Flip gain reference:", default=NO_FLIP,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help="Flip gain reference after rotation.")

        form.addParam('defectFile', params.FileParam, allowsNull=True,
                      label='Camera defects file',
                      help='Defect file that stores entries of defects on camera.\n'
                           'Each entry corresponds to a rectangular region in image. '
                           'The pixels in such a region are replaced by '
                           'neighboring good pixel values. Each entry contains '
                           '4 integers x, y, w, h representing the x, y '
                           'coordinates, width, and height, respectively.')

        form.addParam('extraParams2', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""*Extra parameters:*\n
        -Bft\t\t100\t\tBFactor for alignment, in px^2.\n
        -Iter\t\t5\t\tMaximum iterations for iterative alignment.\n
        -MaskCent\t\t0 0\t\tCenter of subarea that will be used for alignment,
        \t\t\t\tdefault *0 0* corresponding to the frame center.\n
        -MaskSize\t\t1.0 1.0\t\tThe size of subarea that will be used for alignment,
        \t\t\t\tdefault *1.0 1.0* corresponding full size.\n
        -Align\t\t1\t\tGenerate aligned sum (1) or simple sum (0).\n
        -FmRef\t\t-1\t\tSpecify which frame to be the reference to which
        \t\t\t\tall other frames are aligned, by default (-1) the
        \t\t\t\tthe central frame is chosen. The central frame is
        \t\t\t\tat N/2 based upon zero indexing where N is the
        \t\t\t\tnumber of frames that will be summed, i.e., not
        \t\t\t\tincluding the frames thrown away.\n
        -Tilt\t\t0 0\t\tTilt angle range for a dose fractionated tomographic
        \t\t\t\ttilt series, e.g. *-60 60*\n

    Since version *1.1.0*:\n
        -GpuMemUsage\t\t0.5\t\tSpecify how much GPU memory is used to buffer movie frames.
        \t\t\t\tIt is recommended when running side by side processes in
        \t\t\t\tthe same card. By default is 50% (i. e 0.5).\n
        -InFmMotion\t\t1\t\tTakes into account of motion-induced blurring of
        \t\t\t\teach frame. It has shown resolution improvement
        \t\t\t\tin some test cases. By default this option is off.\n
        -Bft\t\t500 150\t\tSince version 1.1.0 this option can take two arguments.
        \t\t\t\tFirst one is used in global-motion measurement and the
        \t\t\t\tsecond one is for local-motion. (default 500 150).""")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def _processMovie(self, movie):
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getRelPath(self._getOutputMicName(movie),
                                       movieFolder)
        outputMovieFn = self._getRelPath(self._getOutputMovieName(movie),
                                         movieFolder)
        movieBaseName = pwutils.removeExt(movie.getFileName())
        aveMicFn = movieBaseName + '_uncorrected_avg.mrc'
        logFile = self._getRelPath(self._getMovieLogFile(movie),
                                   movieFolder)

        a0, aN = self._getRange(movie, 'align')
        program = self._getProgram()

        logFileBase = (logFile.replace('0-Full.log', '').replace(
            '0-Patch-Full.log', ''))
        # default values for motioncor2 are (1, 1)
        cropDimX = self.cropDimX.get() or 1
        cropDimY = self.cropDimY.get() or 1

        numbOfFrames = self._getNumberOfFrames(movie)

        if self.doApplyDoseFilter:
            preExp, dose = self._getCorrectedDose(inputMovies)
        else:
            preExp, dose = 0.0, 0.0

        # reset values = 1 to 0 (motioncor2 does it automatically,
        # but we need to keep this for consistency)
        if self.patchX.get() == 1:
            self.patchX.set(0)
        if self.patchY.get() == 1:
            self.patchY.set(0)

        argsDict = {'-OutMrc': '"%s"' % outputMicFn,
                    '-Patch': '%d %d' % (self.patchX, self.patchY),
                    '-MaskCent': '%d %d' % (self.cropOffsetX,
                                            self.cropOffsetY),
                    '-MaskSize': '%d %d' % (cropDimX, cropDimY),
                    '-FtBin': self.binFactor.get(),
                    '-Tol': self.tol.get(),
                    '-Group': self.group.get(),
                    '-FmDose': dose,
                    '-Throw': '%d' % a0,
                    '-Trunc': '%d' % (abs(aN - numbOfFrames + 1)),
                    '-PixSize': inputMovies.getSamplingRate(),
                    '-kV': inputMovies.getAcquisition().getVoltage(),
                    '-LogFile': logFileBase,
                    }
        argsDict['-InitDose'] = preExp
        argsDict['-OutStack'] = 1 if self.doSaveMovie else 0

        if self.defectFile.get():
            argsDict['-DefectFile'] = self.defectFile.get()

            if self.versionGE('1.0.1'):  # Patch overlap was introduced in 1.0.1
                patchOverlap = self.getAttributeValue('patchOverlap', None)
                if patchOverlap:  # 0 or None is False
                    argsDict['-Patch'] += " %d" % patchOverlap

        if self.doMagCor:
            if self.useEst:
                inputEst = self.inputEst.get().getOutputLog()
                input_params = parseMagEstOutput(inputEst)
                argsDict['-Mag'] = '%0.3f %0.3f %0.3f' % (
                    input_params[1],
                    input_params[2],
                    input_params[0])
            else:
                argsDict['-Mag'] = '%0.3f %0.3f %0.3f' % (self.scaleMaj,
                                                          self.scaleMin,
                                                          self.angDist)

        ext = pwutils.getExt(movie.getFileName()).lower()

        if ext in ['.mrc', '.mrcs']:
            args = ' -InMrc "%s" ' % movie.getBaseName()
        elif ext in ['.tif', '.tiff']:
            args = ' -InTiff "%s" ' % movie.getBaseName()
        else:
            raise Exception("Unsupported format: %s" % ext)

        args += ' '.join(['%s %s' % (k, v)
                          for k, v in argsDict.iteritems()])

        if inputMovies.getGain():
            args += ' -Gain "%s" ' % inputMovies.getGain()
            args += ' -RotGain %d' % self.gainRot.get()
            args += ' -FlipGain %d' % self.gainFlip.get()

        if inputMovies.getDark():
            args += ' -Dark "%s"' % inputMovies.getDark()

        args += ' -Gpu %(GPU)s'
        args += ' ' + self.extraParams2.get()

        try:
            self.runJob(program, args, cwd=movieFolder,
                        env=motioncorr.Plugin.getEnviron())
            self._fixMovie(movie)

            # Compute PSDs
            outMicFn = self._getExtraPath(self._getOutputMicName(movie))
            if not os.path.exists(outMicFn):
                # if only DW mic is saved
                outMicFn = self._getExtraPath(self._getOutputMicWtName(movie))

            def _extraWork():
                if self.doComputePSD:
                    # Compute uncorrected avg mic
                    roi = [self.cropOffsetX.get(), self.cropOffsetY.get(),
                           self.cropDimX.get(), self.cropDimY.get()]
                    fakeShiftsFn = self.writeZeroShifts(movie)
                    # FIXME: implement gain flip/rotation
                    self.averageMovie(movie, fakeShiftsFn, aveMicFn,
                                      binFactor=self.binFactor.get(),
                                      roi=roi, dark=inputMovies.getDark(),
                                      gain=inputMovies.getGain())

                    self.computePSDs(movie, aveMicFn, outMicFn,
                                     outputFnCorrected=self._getPsdJpeg(movie))

                self._saveAlignmentPlots(movie)

                if self._doComputeMicThumbnail():
                    self.computeThumbnail(outMicFn,
                                          outputFn=self._getOutputMicThumbnail(
                                              movie))
                # This protocols cleans up the temporary movie folder
                # which is required mainly when using a thread for this extra work
                self._cleanMovieFolder(movieFolder)

            if self._useWorkerThread():
                thread = Thread(target=_extraWork)
                thread.start()
            else:
                _extraWork()
        except:
            print("ERROR: Movie %s failed\n" % movie.getName())

    def _insertFinalSteps(self, deps):
        stepId = self._insertFunctionStep('waitForThreadStep', prerequisites=deps)
        return [stepId]

    def waitForThreadStep(self):
        # Quick and dirty (maybe desperate) way to wait
        # if the PSD and thumbnail were computed with a thread
        # If running in streaming this will not be necessary
        if self._useWorkerThread():
            time.sleep(60)  # wait 1 min to give some time the thread to finish

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        # Check base validation before the specific ones
        errors = ProtAlignMovies._validate(self)

        if self.doApplyDoseFilter and self.inputMovies.get():
            inputMovies = self.inputMovies.get()
            doseFrame = inputMovies.getAcquisition().getDosePerFrame()

            if doseFrame == 0.0 or doseFrame is None:
                errors.append('Dose per frame for input movies is 0 or not '
                              'set. You cannot apply dose filter.')

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getVersion(self):
        return motioncorr.Plugin.getActiveVersion()

    def _getProgram(self):
        return motioncorr.Plugin.getProgram()

    def _getMovieLogFile(self, movie):
        if self.patchX == 0 and self.patchY == 0:
            return 'micrograph_%06d_0-Full.log' % movie.getObjId()
        else:
            return 'micrograph_%06d_0-Patch-Full.log' % movie.getObjId()

    def _getAbsPath(self, baseName):
        return os.path.abspath(self._getExtraPath(baseName))

    def _getRelPath(self, baseName, refPath):
        return os.path.relpath(self._getExtraPath(baseName), refPath)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd_comparison', 'psd', extra=True)

    def _getPsdJpeg(self, movie):
        return self._getNameExt(movie, '_psd', 'jpeg', extra=True)

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts are in pixels irrespective of any binning.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movie))
        xShifts, yShifts = parseMovieAlignment2(logPath)

        return xShifts, yShifts

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = em.Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = em.Image(location=self._getPsdCorr(movie))
            mic.psdJpeg = em.Image(location=self._getPsdJpeg(movie))
        if self._doComputeMicThumbnail():
            mic.thumbnail = em.Image(
                location=self._getOutputMicThumbnail(movie))

    def _saveAlignmentPlots(self, movie):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _fixMovie(self, movie):
        if self.doSaveMovie:
            outputMicFn = self._getExtraPath(self._getOutputMicName(movie))
            outputMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
            movieFn = outputMicFn.replace('_aligned_mic.mrc',
                                          '_aligned_mic_Stk.mrc')
            pwutils.moveFile(movieFn, outputMovieFn)

        if not self.doSaveUnweightedMic:
            fnToDelete = self._getExtraPath(self._getOutputMicName(movie))
            pwutils.cleanPath(fnToDelete)

    def writeZeroShifts(self, movie):
        # TODO: find another way to do this
        shiftsMd = self._getTmpPath('zero_shifts.xmd')
        pwutils.cleanPath(shiftsMd)
        xshifts = [0] * movie.getNumberOfFrames()
        yshifts = xshifts
        alignment = MovieAlignment(first=1, last=movie.getNumberOfFrames(),
                                   xshifts=xshifts, yshifts=yshifts)
        roiList = [0, 0, 0, 0]
        alignment.setRoi(roiList)
        movie.setAlignment(alignment)
        writeShiftsMovieAlignment(movie, shiftsMd,
                                  1, movie.getNumberOfFrames())
        return shiftsMd

    def _getRange(self, movie, prefix):

        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, prefix)

        if iniFrame != indxFrame:
            first -= iniFrame
            last -= iniFrame
        else:
            first -= 1
            last -= 1

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
            else:
                return movie.getNumberOfFrames()
        else:
            return movie.getNumberOfFrames()

    def _createOutputMicrographs(self):
        createWeighted = self._createOutputWeightedMicrographs()
        # To create the unweighted average micrographs
        # we only consider the 'doSaveUnweightedMic' flag if the
        # weighted ones should be created.
        return (self.doSaveAveMic and
                (not createWeighted or self.doSaveUnweightedMic))

    def _createOutputWeightedMicrographs(self):
        return (self.doSaveAveMic and self.doApplyDoseFilter)

    def _doComputeMicThumbnail(self):
        return (self.doSaveAveMic and self.doComputeMicThumbnail)

    def _useWorkerThread(self):
        return '--use_worker_thread' in self.extraProtocolParams.get()

    def _doMovieFolderCleanUp(self):
        return False


def createGlobalAlignmentPlot(meanX, meanY, first):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Alignment based upon full frames')
    ax.set_xlabel('Shift x (pixels)')
    ax.set_ylabel('Shift y (pixels)')

    i = first
    # The output and log files list the shifts relative to the first frame.
    # or middle frame for motioncor2 1.0.2?

    # ROB unit seems to be pixels since sampling rate is only asked
    # by the program if dose filtering is required
    skipLabels = ceil(len(meanX)/10.0)
    labelTick = 1

    for x, y in izip(meanX, meanY):
        sumMeanX.append(x)
        sumMeanY.append(y)
        if labelTick == 1:
            ax.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    ax.plot(sumMeanX, sumMeanY, color='b')
    ax.plot(sumMeanX, sumMeanY, 'yo')

    plotter.tightLayout()

    return plotter
