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
from threading import Thread

import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_PARALLEL
from pwem.emlib.image import ImageHandler
from pwem.objects import Image, Float
from pwem.protocols import ProtAlignMovies

from .. import Plugin
from ..constants import *
from ..convert import *
from relion.convert.convert31 import OpticsGroups


class ProtMotionCorr(ProtAlignMovies):
    """ This protocol wraps motioncor2 movie alignment program developed at UCSF.

    Motioncor2 performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'movie alignment'
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.isEER = False

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif', '.eer'] else 'mrc'

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
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Generate odd and even sums using odd and even frames '
                               'respectively when this option is enabled.')

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save aligned movie?")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      label="Compute PSD?",
                      help="If Yes, the protocol will compute for each "
                           "aligned micrograph the PSD using EMAN2.")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      default=False,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail with EMAN2 and keep it with the '
                           'micrograph object for visualization purposes. ')

        form.addParam('extraProtocolParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional protocol parameters',
                      help="Here you can provide some extra parameters for the "
                           "protocol, not the underlying motioncor2 program."
                           "You can provide many options separated by space. "
                           "\n\n*Options:* \n\n"
                           "--dont_use_worker_thread \n"
                           " Now by default we use a separate thread to compute"
                           " PSD and thumbnail (if is required). This allows "
                           " more effective use of the GPU card, but requires "
                           " an extra CPU. Use this option (NOT RECOMMENDED) if "
                           " you want to prevent this behaviour")

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

        form.addParam('patchOverlap', params.IntParam, default=0,
                      label='Patches overlap (%)',
                      help='Specify the overlapping between patches. '
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

        form.addParam('doSaveUnweightedMic', params.BooleanParam, default=False,
                      condition='doApplyDoseFilter',
                      label="Save unweighted micrographs?",
                      help="Aligned but non-dose weighted images are sometimes "
                           "useful in CTF estimation, although there is no "
                           "difference in most cases.")

        form.addParam('extraParams2', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="Extra command line parameters. See MotionCor2 help.")

        form.addSection(label="Gain and defects")
        form.addParam('gainRot', params.EnumParam,
                      choices=['no rotation', '90 degrees',
                               '180 degrees', '270 degrees'],
                      label="Rotate gain reference:",
                      default=NO_ROTATION,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help="Rotate gain reference counter-clockwise.")

        form.addParam('gainFlip', params.EnumParam,
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

        form.addParam('defectMap', params.FileParam, allowsNull=True,
                      label='Camera defects map',
                      help='1. Defect map is a binary (0 or 1) map where defective '
                           ' pixels are assigned value of 1 and good pixels have '
                           'value of 0.\n2. The defective pixels are corrected '
                           'with a random pick of good pixels in its neighborhood. '
                           '\n3. This is map must have the same dimension and '
                           'orientation as the input movie frame.\n4. This map '
                           'can be provided as either MRC or TIFF file that has '
                           'MRC mode of 0 or 5 (unsigned 8 bit).')

        form.addSection("EER")
        form.addParam('EERtext', params.LabelParam,
                      label="These options are ignored for non-EER movies.")
        form.addParam('eerGroup', params.IntParam, default=32,
                      label='EER fractionation',
                      help="The number of hardware frames to group into one "
                           "fraction. This option is relevant only for Falcon4 "
                           "movies in the EER format. Falcon 4 operates at "
                           "248 frames/s.\nFractionate such that each fraction "
                           "has about 0.5 to 1.25 e/A2.")
        form.addParam('eerSampling', params.EnumParam, default=0,
                      choices=['1x', '2x', '4x'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='EER upsampling',
                      help="EER upsampling (1x = 4K, 2x = 8K, 3x=16K)")

        form.addSection(label="Mag. correction")
        form.addParam('doMagCor', params.BooleanParam, default=False,
                      label='Correct anisotropic magnification?',
                      help='Correct anisotropic magnification by '
                           'stretching image along the major axis, '
                           'the axis where the lower magnification is '
                           'detected.')
        form.addParam('scaleMaj', params.FloatParam, default=1.0,
                      condition='doMagCor',
                      label='Major scale factor')
        form.addParam('scaleMin', params.FloatParam, default=1.0,
                      condition='doMagCor',
                      label='Minor scale factor')
        form.addParam('angDist', params.FloatParam, default=0.0,
                      condition='doMagCor',
                      label='Distortion angle (deg)')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def _convertInputStep(self):
        inputMovies = self.getInputMovies()
        # parse EER gain file before its conversion to mrc
        if self.isEER and inputMovies.getGain():
            defects = parseEERDefects(inputMovies.getGain())
            if defects:
                with open(self._getExtraPath("defects_eer.txt"), "w") as f:
                    for d in defects:
                        f.write(" ".join(str(i) for i in d) + "\n")

        ProtAlignMovies._convertInputStep(self)
        if self.isEER:
            # write FmIntFile
            _, numbOfFrames, _ = inputMovies.getFramesRange()
            if self.doApplyDoseFilter:
                _, dose = self._getCorrectedDose(inputMovies)
            else:
                dose = 0.0
            with open(self._getExtraPath("FmIntFile.txt"), "w") as f:
                f.write("%d %d %f" % (numbOfFrames, self.eerGroup.get(), dose))

    def _processMovie(self, movie):
        inputMovies = self.getInputMovies()
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getOutputMicName(movie)
        program = Plugin.getProgram()
        frame0, frameN = self._getFrameRange()
        _, numbOfFrames, _ = inputMovies.getFramesRange()

        argsDict = self._getArgs()
        argsDict.update({'-OutMrc': '"%s"' % outputMicFn,
                         '-Throw': '%d' % 0 if self.isEER else (frame0 - 1),
                         '-Trunc': '%d' % 0 if self.isEER else (numbOfFrames - frameN)
                         })

        if Plugin.versionGE('1.4.7'):
            argsDict.update({'-LogDir': './'})
        else:
            logFileBase = self._getMovieRoot(movie) + "_"
            argsDict.update({'-LogFile': logFileBase})

        if self.splitEvenOdd:
            argsDict['-SplitSum'] = 1

        ext = pwutils.getExt(movie.getFileName()).lower()
        if ext in ['.mrc', '.mrcs']:
            args = ' -InMrc "%s" ' % movie.getBaseName()
        elif ext in ['.tif', '.tiff']:
            args = ' -InTiff "%s" ' % movie.getBaseName()
        elif ext == '.eer':
            args = ' -InEer "%s" ' % movie.getBaseName()
        else:
            raise Exception("Unsupported format: %s" % ext)

        args += ' '.join(['%s %s' % (k, v)
                          for k, v in argsDict.items()])
        args += ' ' + self.extraParams2.get()

        try:
            self.runJob(program, args, cwd=movieFolder,
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
                    self.error("ERROR: Extra work "
                               "(i.e plots, PSD, thumbnail) has failed for %s\n"
                               % movie.getFileName())

            if self._useWorkerThread():
                thread = Thread(target=_extraWork)
                thread.start()
            else:
                _extraWork()

        except Exception as e:
            self.error("ERROR: Motioncor2 has failed for %s. --> %s\n" % (movie.getFileName(), str(e)))
            import traceback
            traceback.print_exc()

    def _insertFinalSteps(self, deps):
        stepsId = []
        if self._useWorkerThread():
            stepsId.append(self._insertFunctionStep('waitForThreadStep',
                                                    prerequisites=deps))
        return stepsId

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
            summary.append('Aligned %d movies using motioncor2.'
                           % self.inputMovies.get().getSize())
            if self.splitEvenOdd and self._createOutputWeightedMicrographs():
                summary.append('Even/odd outputs are dose-weighted!')
        else:
            summary.append('Output is not ready')

        return summary

    def _methods(self):
        methods = []

        if self.doApplyDoseFilter:
            methods.append(' - Applied dose filtering')
        if self.patchX > 1 and self.patchY > 1:
            methods.append(' - Used patch-based alignment')
        if self.group > 1:
            methods.append(' - Grouped %d frames' % self.group)

        return methods

    def _validate(self):
        errors = []
        inputMovies = self.inputMovies.get()

        # check if the first movie exists
        firstMovie = self.inputMovies.get().getFirstItem()
        if not os.path.exists(firstMovie.getFileName()):
            errors.append("The input movie files do not exist!!! "
                          "Since usually input movie files are symbolic links, "
                          "please check that links are not broken if you "
                          "moved the project folder. ")

        # check frames range
        _, lastFrame, _ = inputMovies.getFramesRange()
        self.isEER = pwutils.getExt(firstMovie.getFileName()) == ".eer"
        if self.isEER:
            if not Plugin.versionGE('1.4.0'):
                errors.append("EER is only supported for motioncor2 v1.4.0 or newer.")
            if self.alignFrame0.get() != 1 or self.alignFrameN.get() not in [0, lastFrame]:
                errors.append("For EER data please set frame range from 1 to 0 (or 1 to %d)." % lastFrame)

        msg = "Frames range must be within %d - %d" % (1, lastFrame)
        if self.alignFrameN.get() == 0:
            self.alignFrameN.set(lastFrame)

        if not (1 <= self.alignFrame0 < lastFrame):
            errors.append(msg)
        elif not (self.alignFrameN <= lastFrame):
            errors.append(msg)
        elif not (self.alignFrame0 < self.alignFrameN):
            errors.append(msg)

        # check dose for DW
        acq = inputMovies.getAcquisition()
        if self.doApplyDoseFilter:
            dose = acq.getDosePerFrame()
            if dose is None or dose < 0.00001:
                errors.append("Input movies do not contain the dose per frame, "
                              "dose-weighting can not be performed.")

        # check eman2 plugin
        if self.doComputeMicThumbnail or self.doComputePSD:
            try:
                from pwem import Domain
                eman2 = Domain.importFromPlugin('eman2', doRaise=True)
            except:
                errors.append("EMAN2 plugin not found!\nComputing thumbnails "
                              "or PSD requires EMAN2 plugin and binaries installed.")

        # check gain dimensions and extension
        if inputMovies.getGain():
            ih = ImageHandler()
            gain = inputMovies.getGain()
            gainx, gainy, _, _ = ih.getDimensions(gain)
            movie = firstMovie.getFileName()
            imgx, imgy, _, _ = ih.getDimensions(movie)

            if sorted([gainx, gainy]) != sorted([imgx, imgy]):
                errors.append("Gain image dimensions (%d x %d) "
                              "do not match the movies (%d x %d)!" %
                              (gainx, gainy, imgx, imgy))

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getArgs(self):
        """ Prepare most arguments for mc2 run. """
        inputMovies = self.getInputMovies()
        # default values for motioncor2 are (1, 1)
        cropDimX = self.cropDimX.get() or 1
        cropDimY = self.cropDimY.get() or 1

        # reset values = 1 to 0 (motioncor2 does it automatically,
        # but we need to keep this for consistency)
        if self.patchX.get() == 1:
            self.patchX.set(0)
        if self.patchY.get() == 1:
            self.patchY.set(0)

        argsDict = {'-Patch': '%d %d' % (self.patchX, self.patchY),
                    '-MaskCent': '%d %d' % (self.cropOffsetX,
                                            self.cropOffsetY),
                    '-MaskSize': '%d %d' % (cropDimX, cropDimY),
                    '-FtBin': self.binFactor.get(),
                    '-Tol': self.tol.get(),
                    '-Group': self.group.get(),
                    '-PixSize': inputMovies.getSamplingRate(),
                    '-kV': inputMovies.getAcquisition().getVoltage(),
                    '-OutStack': 1 if self.doSaveMovie else 0,
                    '-Gpu': '%(GPU)s',
                    '-SumRange': "0.0 0.0",  # switch off writing out DWS,
                    # '-FmRef': 0
                    }

        if self.isEER:
            argsDict['-EerSampling'] = self.eerSampling.get() + 1
            argsDict['-FmIntFile'] = "../../extra/FmIntFile.txt"
        elif self.doApplyDoseFilter:
            preExp, dose = self._getCorrectedDose(inputMovies)
            argsDict.update({'-FmDose': dose,
                             '-InitDose': preExp})

        if self.defectFile.get():
            argsDict['-DefectFile'] = "%s" % self.defectFile.get()
        elif self.defectMap.get():
            argsDict['-DefectMap'] = "%s" % self.defectMap.get()
        elif os.path.exists(self._getExtraPath("defects_eer.txt")):
            argsDict['-DefectFile'] = "../../extra/defects_eer.txt"

        patchOverlap = self.getAttributeValue('patchOverlap')
        if patchOverlap:  # 0 or None is False
            argsDict['-Patch'] += " %d" % patchOverlap

        if self.doMagCor:
            argsDict['-Mag'] = '%0.3f %0.3f %0.3f' % (self.scaleMaj,
                                                      self.scaleMin,
                                                      self.angDist)

        if inputMovies.getGain():
            argsDict.update({'-Gain': '"%s"' % inputMovies.getGain(),
                             '-RotGain': self.gainRot.get(),
                             '-FlipGain': self.gainFlip.get()})

        if inputMovies.getDark():
            argsDict['-Dark'] = "%s" % inputMovies.getDark()

        return argsDict

    def _getBinFactor(self):
        if not self.isEER:
            return self.binFactor.get()
        else:
            return self.binFactor.get() / (self.eerSampling.get() + 1)

    def _getCwdPath(self, movie, path):
        return os.path.join(self._getOutputMovieFolder(movie), path)

    def _getMovieLogFile(self, movie):
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s%s%s-Full.log' % (self._getMovieRoot(movie),
                                    self._getLogSuffix(),
                                    '-Patch' if usePatches else '')

    def _getLogSuffix(self):
        return '_aligned_mic' if Plugin.versionGE('1.4.7') else '_0'

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
        args = "%s %s " % (inputFn, outputFn)
        args += "--process=math.realtofft --meanshrink %s " % scaleFactor
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

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts are in pixels irrespective of any binning.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movie))
        xShifts, yShifts = parseMovieAlignment2(logPath)

        return xShifts, yShifts

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

    def _saveAlignmentPlots(self, movie, pixSize):
        """ Compute alignment shift plots and save to file as png images. """
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange()
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
        a0, aN = self._getFrameRange()
        nframes = aN - a0 + 1
        preExp, dose = self._getCorrectedDose(self.getInputMovies())
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
            self.error("Expected %d frames, found less. Check movie %s !" % (
                nframes, movie.getFileName()))

    def _getFrameRange(self, n=None, prefix=None):
        # Reimplement this method
        if self.isEER:
            return self.alignFrame0.get(), self.alignFrameN.get() // self.eerGroup.get()
        else:
            return self.alignFrame0.get(), self.alignFrameN.get()


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
