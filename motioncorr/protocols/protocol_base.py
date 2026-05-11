# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
from os.path import exists, abspath, basename

import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
from pwem.emlib.image.image_readers import Dm4ImageReader
from pyworkflow.protocol import STEPS_PARALLEL
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.objects import Movie
from pyworkflow.utils import cyanStr
from .. import Plugin, DEFECTS_FILE_EER

from ..constants import NO_FLIP, NO_ROTATION
from ..convert import parseMovieAlignment2, parseEERDefects


logger = logging.getLogger(__name__)


class ProtMotionCorrBase(EMProtocol):
    """
    Base protocol for movie motion correction using MotionCor-based alignment.

    AI Generated:

    Motion Correction Base (ProtMotionCorrBase) — User Manual
        Overview

        The ProtMotionCorrBase protocol provides the common framework for
        correcting beam-induced motion in cryo-EM movie stacks. During electron
        exposure, particles move slightly due to specimen drift and beam-induced
        deformation. If uncorrected, this motion reduces high-resolution signal
        and degrades downstream processing.

        This protocol aligns movie frames, optionally applies dose-weighting,
        supports local patch-based correction, handles gain and dark reference
        corrections, accepts defect information, and supports modern EER
        (Electron Event Representation) movies.

        It acts as a reusable base protocol: subclasses define how input movies
        are accessed and how per-movie outputs are organized, while this class
        provides the shared logic for validation, preprocessing, parameter
        handling, and MotionCor argument generation.

        Inputs

        The protocol expects an input set of movies. Each movie should contain:

            - A valid stack of frames.
            - Acquisition metadata such as voltage and dose per frame.
            - Optionally, gain and dark reference images.

        The input movies must exist on disk. Since cryo-EM projects often rely
        on symbolic links, broken links are explicitly checked during validation.

        Main Workflow

        The protocol follows these general stages:

            1. Input preparation
               Gain, dark, and EER-related files are prepared and converted
               if necessary.

            2. Validation
               Frame ranges, gain dimensions, dose metadata, and EER-specific
               restrictions are checked before execution.

            3. MotionCor argument generation
               All alignment parameters are translated into command-line options
               for the MotionCor binary.

            4. Execution
               Subclasses launch the external program.

            5. Output parsing
               Frame shifts and alignment metadata can later be extracted.

        Frame Alignment

        Motion correction aligns frames over a user-defined frame range.

        For conventional movie stacks:

            - Alignment uses the selected start and end frame directly.

        For EER movies:

            - Alignment is performed on grouped hardware frames.
            - Frame range must begin at frame 1.
            - End frame must be either 0 (all frames) or the total number of
              available frames.

        This restriction is enforced because EER fractionation is internally
        handled differently from standard frame stacks.

        Patch-Based Local Alignment

        Local alignment can be performed by dividing each frame into patches.

        Parameters include:

            - patchX, patchY
                Number of patches along X and Y.

            - patchOverlap
                Percentage of overlap between neighboring patches.

        Practical interpretation:

            - 0 0 patches:
                Only global motion correction is performed.

            - Multiple patches:
                Local deformation is estimated and corrected.

        For biological datasets:

            - Smaller particles or thin ice often benefit from local alignment.
            - Large particles with strong beam-induced local motion usually
              benefit substantially from patch correction.

        Frame Grouping

        Frames can be grouped before alignment.

            - group
                Grouping for global alignment.

            - groupLocal
                Grouping for local alignment.

        Grouping reduces noise by summing neighboring frames prior to alignment.
        This is often useful when individual frames have very low signal.

        Dose Weighting

        If enabled, the protocol applies dose-dependent filtering before summing
        aligned frames.

        This requires valid dose metadata:

            - Initial dose
            - Dose per frame

        Biological rationale:

            Early frames usually preserve higher-resolution information, whereas
            later frames accumulate radiation damage. Dose-weighting downweights
            damaged high-frequency signal and generally improves reconstruction
            quality.

        For tomography-specific acquisition orders, the protocol also supports
        corrected pre-exposure calculation.

        Gain, Dark, and Detector Defects

        Gain and dark correction references are supported.

        Special behavior:

            - DM4 gain references are automatically converted to MRC.
            - Converted files are cached to avoid repeated conversion.

        Defect correction can be provided in three ways:

            - Defect file
                Rectangular bad-pixel regions.

            - Defect map
                Binary defect mask.

            - EER-derived defect file
                Automatically parsed from EER gain files.

        These corrections are important because detector defects can introduce
        artifacts that bias alignment and downstream CTF estimation.

        EER Support

        EER movies are handled explicitly.

        Additional EER-specific behavior includes:

            - Parsing detector defects from the gain file.
            - Writing an FmIntFile describing:
                * total number of frames
                * EER grouping factor
                * dose information

        Important EER parameters:

            - eerGroup
                Hardware frame fractionation.

            - eerSampling
                Upsampling level.

        Since Falcon EER data can have very high temporal sampling, correct
        grouping is important for balancing temporal resolution and signal.

        Magnification Correction

        The protocol optionally supports anisotropic magnification correction.

        Parameters include:

            - scaleMaj
                Major scale factor.

            - scaleMin
                Minor scale factor.

            - angDist
                Distortion angle.

        This correction compensates for directional magnification distortions,
        which can otherwise slightly bias high-resolution reconstructions.

        Validation Rules

        Before execution, the protocol verifies:

            - Input movies exist.
            - Frame range is valid.
            - EER frame constraints are respected.
            - Dose metadata is available if dose-weighting is enabled.
            - Gain dimensions match movie dimensions.

        These checks prevent common user errors that would otherwise produce
        invalid motion correction results.

        MotionCor Argument Construction

        The protocol prepares command-line arguments for MotionCor.

        Internally generated arguments include:

            - frame throwing and truncation
            - patch geometry
            - crop center and crop size
            - binning factor
            - tolerance
            - pixel size
            - accelerating voltage
            - dose-weighting parameters
            - gain and dark references
            - defect handling
            - EER parameters
            - anisotropic magnification correction

        This centralized argument generation is one of the main roles of the
        base class.

        Output Utilities

        After execution, subclasses can use shared utilities to retrieve:

            - movie log file paths
            - per-frame X/Y shifts
            - corrected frame ranges
            - corrected dose parameters
            - effective binning factors

        The extracted shifts are always reported in pixels, independent of any
        binning applied during alignment.

        Practical Recommendations

        For routine cryo-EM single-particle processing:

            - Use dose-weighting whenever dose metadata is available.
            - Use local patch alignment for most modern datasets.
            - Keep frame grouping modest unless frames are extremely noisy.

        For EER datasets:

            - Carefully choose EER grouping so each fraction contains a useful
              electron dose.
            - Ensure the EER frame range follows the required convention.

        For high-resolution work:

            - Check gain orientation carefully.
            - Verify gain dimensions.
            - Use anisotropic magnification correction when available.

        Final Perspective

        Motion correction is one of the earliest but most important stages of
        cryo-EM image processing.

        The ProtMotionCorrBase protocol provides the shared infrastructure
        needed to prepare movie stacks for accurate downstream analysis.
        Proper frame range selection, correct detector references, reliable
        dose metadata, and sensible local alignment settings all directly
        influence the quality of particle alignment, CTF estimation, and final
        reconstruction resolution.
    """
    _label = None
    stepsExecutionMode = STEPS_PARALLEL
    program = Plugin.getProgram()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isEER = False

    # -------------------------- DEFINE param functions -----------------------
    @staticmethod
    def _defineCommonParams(form, allowDW=True):
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=cons.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " Motioncor can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addSection(label="Motioncor params")
        if allowDW:
            form.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                          label='Apply dose filter',
                          help='Apply a dose-dependent filter to frames before '
                               'summing them. Pre-exposure and dose per frame '
                               'should be specified during movies import.')

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

        line = form.addLine('Group N frames',
                            help='Group every specified number of frames by adding '
                                 'them together. The alignment is then performed on '
                                 'the summed frames. By default, no grouping is '
                                 'performed.')
        line.addParam('group', params.IntParam, default='1',
                      label='global align')
        line.addParam('groupLocal', params.IntParam, default='4',
                      label='local align')

        form.addParam('tol', params.FloatParam, default='0.2',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Tolerance (px)',
                      help='Tolerance for iterative alignment, default *0.2px*.')
        if allowDW:
            form.addParam('doSaveUnweightedMic', params.BooleanParam, default=False,
                          condition='doApplyDoseFilter',
                          label="Save unweighted images?",
                          help="Aligned but non-dose weighted images are sometimes "
                               "useful in CTF estimation, although there is no "
                               "difference in most cases.")

        form.addParam('extraParams2', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="Extra command line parameters. See MotionCor help.")

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
                           'coordinates of the lower left corner, width, and '
                           'height, respectively.')

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
                           "fraction. This option is relevant only for Falcon 4 "
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
        logger.info(cyanStr('Converting the inputs (gain, dark and /or EER)...'))
        inputMovies = self.getInputMovies()
        self._prepareEERFiles()
        pwutils.makePath(self._getExtraPath('DONE'))

        # Convert gain
        gain = inputMovies.getGain()
        inputMovies.setGain(self.__convertCorrectionImage(gain))

        # Convert dark
        dark = inputMovies.getDark()
        inputMovies.setDark(self.__convertCorrectionImage(dark))

    def _prepareEERFiles(self):
        """ Parse .gain file for defects and create dose distribution file.
        EER gain must be parsed before conversion to mrc. """
        gainFile = self.getInputMovies().getGain()
        if self.isEER and gainFile:
            defects = parseEERDefects(gainFile)
            if defects:
                with open(self._getExtraPath(DEFECTS_FILE_EER), "w") as f:
                    for d in defects:
                        f.write(" ".join(str(i) for i in d) + "\n")
        if self.isEER:
            # write FmIntFile
            numbOfFrames = self._getNumberOfFrames()
            if self.doApplyDoseFilter:
                if self.getClassName() == "ProtTsMotionCorr":
                    acqOrder = 1
                else:
                    acqOrder = None
                _, dose = self._getCorrectedDose(acqOrder)
            else:
                dose = 0.0
            with open(self._getExtraPath("FmIntFile.txt"), "w") as f:
                f.write(f"{numbOfFrames} {self.eerGroup.get()} {dose}")

    def allowsDelete(self, obj):
        return True

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        inputMovies = self.getInputMovies()

        # check if the first movie exists
        firstMovie = inputMovies.getFirstItem()
        if not isinstance(firstMovie, Movie):
            firstMovie = firstMovie.getFirstItem()

        if not exists(firstMovie.getFileName()):
            errors.append("The input movie files do not exist!!! "
                          "Since usually input movie files are symbolic links, "
                          "please check that links are not broken if you "
                          "moved the project folder. ")

        # check frames range
        lastFrame = self._getNumberOfFrames()
        self.isEER = pwutils.getExt(firstMovie.getFileName()) == ".eer"
        if self.isEER:
            if self.alignFrame0.get() != 1 or self.alignFrameN.get() not in [0, lastFrame]:
                errors.append(f"For EER data please set frame range "
                              f"from 1 to 0 (or 1 to {lastFrame}).")

        msg = f"Frames range must be within 1 - {lastFrame}"
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
                errors.append("Input movies do not contain the dose information, "
                              "dose-weighting can not be performed.")

        # check gain dimensions and extension
        if inputMovies.getGain() and exists(inputMovies.getGain()):
            ih = ImageHandler()
            gain = inputMovies.getGain()
            gainx, gainy, _, _ = ih.getDimensions(gain)
            movie = firstMovie.getFileName()
            imgx, imgy, _, _ = ih.getDimensions(movie)

            if sorted([gainx, gainy]) != sorted([imgx, imgy]):
                errors.append(f"Gain image dimensions ({gainx} x {gainy}) "
                              f"do not match the movies ({imgx} x {imgy})!")

        return errors

    def _methods(self):
        methods = []

        if self.doApplyDoseFilter:
            methods.append(' - Applied dose filtering')
        if self.patchX > 1 and self.patchY > 1:
            methods.append(' - Used patch-based alignment')
        if self.group > 1:
            methods.append(f' - Grouped {self.group} frames')

        return methods

    # --------------------------- UTILS functions -----------------------------
    def getInputMovies(self):
        """ Should be implemented in subclasses. """
        raise NotImplementedError

    def _getMcArgs(self, acqOrder=None):
        """ Prepare most arguments for the binary. """
        inputMovies = self.getInputMovies()

        # default values for motioncor are (1, 1)
        cropDimX = self.cropDimX.get() or 1
        cropDimY = self.cropDimY.get() or 1

        frame0, frameN = self._getFramesRange()
        numbOfFrames = self._getNumberOfFrames()

        # reset values = 1 to 0 (motioncor does it automatically,
        # but we need to keep this for consistency)
        if self.patchX.get() == 1:
            self.patchX.set(0)
        if self.patchY.get() == 1:
            self.patchY.set(0)

        argsDict = {
            '-Throw': 0 if self.isEER else (frame0 - 1),
            '-Trunc': 0 if self.isEER else (numbOfFrames - frameN),
            '-Patch': f"{self.patchX} {self.patchY}",
            '-MaskCent': f"{self.cropOffsetX} {self.cropOffsetY}",
            '-MaskSize': f"{cropDimX} {cropDimY}",
            '-FtBin': self.binFactor.get(),
            '-Tol': self.tol.get(),
            '-PixSize': inputMovies.getSamplingRate(),
            '-kV': inputMovies.getAcquisition().getVoltage(),
            '-Cs': 0,
            '-OutStack': 1 if self.doSaveMovie else 0,
            '-Gpu': '%(GPU)s',
            '-SumRange': "0.0 0.0",  # switch off writing out DWS,
            # '-LogDir': './'
            # '-FmRef': 0
        }
        if self.isEER:
            argsDict.update({'-EerSampling': self.eerSampling.get() + 1,
                             '-FmIntFile': self._getExtraPath("FmIntFile.txt")})

        if self.doApplyDoseFilter:
            preExp, dose = self._getCorrectedDose(acqOrder)
            argsDict['-InitDose'] = preExp if preExp > 0.001 else 0
            if not self.isEER:
                argsDict['-FmDose'] = dose

        argsDict['-Group'] = f'{self.group.get()} {self.groupLocal.get()}'

        if self.splitEvenOdd:
            argsDict['-SplitSum'] = 1
        if self.defectFile.get():
            argsDict['-DefectFile'] = self.defectFile.get()
        elif self.defectMap.get():
            argsDict['-DefectMap'] = self.defectMap.get()
        elif exists(self._getExtraPath(DEFECTS_FILE_EER)):
            argsDict['-DefectFile'] = self._getExtraPath(DEFECTS_FILE_EER)

        if inputMovies.getGain():
            argsDict.update({'-Gain': f'"{inputMovies.getGain()}"',
                             '-RotGain': self.gainRot.get(),
                             '-FlipGain': self.gainFlip.get()})

        if inputMovies.getDark():
            argsDict['-Dark'] = inputMovies.getDark()

        patchOverlap = self.getAttributeValue('patchOverlap')
        if patchOverlap:  # 0 or None is False
            argsDict['-Patch'] += f" {patchOverlap}"

        if self.doMagCor:
            argsDict['-Mag'] = f"{self.scaleMaj} {self.scaleMin} {self.angDist}"

        return argsDict

    @staticmethod
    def _getInputFormat(inputFn, absPath=False):
        if absPath:
            inputFn = abspath(inputFn)
        else:
            inputFn = basename(inputFn)
        ext = pwutils.getExt(inputFn).lower()
        if ext in ['.mrc', '.mrcs']:
            args = f' -InMrc "{inputFn}" '
        elif ext in ['.tif', '.tiff']:
            args = f' -InTiff "{inputFn}" '
        elif ext == '.eer':
            args = f' -InEer "{inputFn}" '
        else:
            raise ValueError(f"Unsupported format: {ext}")

        return args

    def _getFramesRange(self):
        """ Returns frames range for alignment. """
        if self.isEER:
            return self.alignFrame0.get(), self.alignFrameN.get() // self.eerGroup.get()
        else:
            return self.alignFrame0.get(), self.alignFrameN.get()

    def _getBinFactor(self):
        # Reimplement this method
        if not self.isEER:
            return self.binFactor.get()
        else:
            return self.binFactor.get() / (self.eerSampling.get() + 1)

    def _getMovieLogFile(self, movie):
        """ Should be implemented in subclasses. """
        raise NotImplementedError

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts are in pixels irrespective of any binning.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movie))
        xShifts, yShifts = parseMovieAlignment2(logPath)

        return xShifts, yShifts

    def _getNumberOfFrames(self):
        """ Dirty hack because of https://github.com/scipion-em/scipion-em-tomo/issues/334 """
        _, frames, _ = self.__getFramesRange()
        if not frames:
            frames = self.getInputMovies().getFirstItem().getDim()[2]

        return frames

    def _getDoseParams(self):
        """ Precalculate params in advance. """
        acq = self.getInputMovies().getAcquisition()
        preExp = acq.getDoseInitial()
        dose = acq.getDosePerFrame()

        return acq, preExp, dose

    def __getFramesRange(self):
        """ Returns frames range for input movies. """
        return self.getInputMovies().getFramesRange()

    def _getCorrectedDose(self, acqOrder=None):
        """ Reimplement this because of a special tomo case. """
        acq, preExp, dose = self._getDoseParams()

        if acqOrder is not None:
            # in tomo case dose = dosePerTilt
            preExp += (acqOrder - 1) * dose
            dosePerFrame = dose / self._getNumberOfFrames() if dose else 0.0
            return preExp, dosePerFrame
        else:
            firstFrame, _, _ = self.__getFramesRange()
            preExp += dose * (firstFrame - 1)
            return preExp, dose

    def __convertCorrectionImage(self, image):
        """ Reimplement to convert dm4 gain only. """
        if image is None:
            return None
        elif image.endswith(".dm4"):
            # Get final correction image file
            finalName = self._getExtraPath(pwutils.replaceBaseExt(image, "mrc"))

            if not exists(finalName):
                Dm4ImageReader.dmToMrc(image, finalName)
                logger.info(cyanStr(f"Converting {image} to {finalName}"))

            # return final name
            return abspath(finalName)
        else:
            return image
