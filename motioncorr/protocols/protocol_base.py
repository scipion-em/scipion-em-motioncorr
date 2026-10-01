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
        form.addParam('eerGroup', params.IntParam, default=8,
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

        # Convert correction images once and keep the prepared paths on
        # the protocol. The input Set is periodically refreshed during
        # streaming and may legitimately restore its persisted original
        # gain/dark metadata; processing must keep using the correction
        # images prepared by this shared step.
        gain = inputMovies.getGain()
        self._preparedGain = self.__convertCorrectionImage(gain)
        inputMovies.setGain(self._preparedGain)

        dark = inputMovies.getDark()
        self._preparedDark = self.__convertCorrectionImage(dark)
        inputMovies.setDark(self._preparedDark)

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

        # Streaming inputs may legitimately still be empty when
        # validation runs. Defer movie-dependent checks until the first
        # item exists instead of treating missing-yet metadata as invalid.
        firstMovie = inputMovies.getFirstItem()
        if firstMovie is None:
            isStreamOpen = getattr(inputMovies, 'isStreamOpen', None)
            if callable(isStreamOpen) and isStreamOpen():
                return errors
            errors.append("Input movies are empty.")
            return errors

        if not isinstance(firstMovie, Movie):
            getFirstItem = getattr(firstMovie, 'getFirstItem', None)
            firstMovie = getFirstItem() if callable(getFirstItem) else firstMovie

        if firstMovie is None:
            errors.append("Input movies are empty.")
            return errors

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
        if self.alignFrameN.get() > 0:
            if not (1 <= self.alignFrame0 < lastFrame):
                errors.append(msg)
            if self.alignFrameN > lastFrame:
                errors.append(msg)
            if self.alignFrame0 >= self.alignFrameN:
                errors.append(msg)

        # check dose for DW
        if self.doApplyDoseFilter:
            _, _, dose = self._resolveAcquisitionValues(
                inputMovies=inputMovies,
                firstMovie=firstMovie,
            )
            if dose is None or dose < 0.00001:
                errors.append(
                    "Input movies do not contain the dose information, "
                    "dose-weighting can not be performed. Either provide "
                    "the dose per frame in the input movies' acquisition, "
                    "or turn off 'Apply dose filtering'."
                )

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
        # Trunc value
        truncVal = numbOfFrames - frameN
        if self.isEER or frameN == 0:
            truncVal = 0
        argsDict = {
            '-Throw': 0 if self.isEER else (frame0 - 1),
            '-Trunc': truncVal,
            '-Patch': f"{self.patchX} {self.patchY}",
            '-MaskCent': f"{self.cropOffsetX} {self.cropOffsetY}",
            '-MaskSize': f"{cropDimX} {cropDimY}",
            '-FtBin': self.binFactor.get(),
            '-Tol': self.tol.get(),
            '-PixSize': inputMovies.getSamplingRate(),
            '-kV': self._getCachedAcquisitionValues()[0],
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

        if self.doApplyDoseFilter and self._hasValidDose():
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

        gain = getattr(
            self,
            '_preparedGain',
            inputMovies.getGain(),
        )
        dark = getattr(
            self,
            '_preparedDark',
            inputMovies.getDark(),
        )

        if gain:
            argsDict.update({'-Gain': f'"{gain}"',
                             '-RotGain': self.gainRot.get(),
                             '-FlipGain': self.gainFlip.get()})

        if dark:
            argsDict['-Dark'] = dark

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
        voltage, preExp, dose = self._getCachedAcquisitionValues()

        return voltage, preExp, dose

    def __getFramesRange(self):
        """ Returns frames range for input movies. """
        return self.getInputMovies().getFramesRange()

    @staticmethod
    def _readAcquisitionValues(acquisition):
        if acquisition is None:
            return None, None, None

        def _get(name):
            getter = getattr(acquisition, name, None)
            return getter() if callable(getter) else None

        return (
            _get('getVoltage'),
            _get('getDoseInitial'),
            _get('getDosePerFrame'),
        )

    def _resolveAcquisitionValues(self, inputMovies=None, firstMovie=None):
        # Resolve Set-level Acquisition metadata first. When a streaming
        # runtime Set is visible before all nested properties have arrived,
        # use the first input item only to fill still-missing values.
        inputMovies = inputMovies or self.getInputMovies()

        getAcquisition = getattr(inputMovies, 'getAcquisition', None)
        setAcquisition = (
            getAcquisition()
            if callable(getAcquisition)
            else None
        )
        values = list(
            self._readAcquisitionValues(setAcquisition)
        )

        if None not in values:
            return tuple(values)

        if firstMovie is None:
            getFirstItem = getattr(inputMovies, 'getFirstItem', None)
            firstMovie = (
                getFirstItem()
                if callable(getFirstItem)
                else None
            )

        candidate = firstMovie
        getCandidateAcquisition = getattr(
            candidate,
            'getAcquisition',
            None,
        )

        if not callable(getCandidateAcquisition):
            getFirstItem = getattr(candidate, 'getFirstItem', None)
            candidate = (
                getFirstItem()
                if callable(getFirstItem)
                else candidate
            )
            getCandidateAcquisition = getattr(
                candidate,
                'getAcquisition',
                None,
            )

        itemAcquisition = (
            getCandidateAcquisition()
            if callable(getCandidateAcquisition)
            else None
        )
        itemValues = self._readAcquisitionValues(
            itemAcquisition
        )

        return tuple(
            itemValue if value is None else value
            for value, itemValue in zip(values, itemValues)
        )

    def _getCachedAcquisitionValues(self):
        # Cache plain scalars so later loadAllProperties() refreshes cannot
        # mutate the Acquisition object used by an already-running protocol.
        if not hasattr(self, '_cachedAcquisitionValues'):
            with self._lock:
                if not hasattr(self, '_cachedAcquisitionValues'):
                    self._cachedAcquisitionValues = (
                        self._resolveAcquisitionValues()
                    )

        return self._cachedAcquisitionValues

    def _useDoseWeightedOutput(self):
        doApplyDoseFilter = self.doApplyDoseFilter
        getter = getattr(doApplyDoseFilter, 'get', None)
        enabled = (
            getter()
            if callable(getter)
            else bool(doApplyDoseFilter)
        )
        return bool(enabled) and self._hasValidDose()

    def _hasValidDose(self):
        """ Whether the input movies' acquisition carries a usable dose
        per frame, so dose-weighting can actually be requested from the
        underlying tool. _validate() already blocks a direct launch
        with doApplyDoseFilter on and no dose, but that check is not
        necessarily re-enforced on every launch path (e.g. a protocol
        started as part of a resumed/chained workflow) - asking
        MotionCor2 to dose-weight with an unusable dose makes it
        silently skip writing the dose-weighted output altogether, so
        this must also be checked at the point the args are built.

        Computed once and cached: dose availability is a property of
        the whole input Set, constant for the entire run (same as
        self.isEER), not something that can legitimately change
        between movies. """
        if not hasattr(self, '_hasValidDoseCache'):
            _, _, dose = self._getCachedAcquisitionValues()
            self._hasValidDoseCache = dose is not None and dose >= 0.00001

        return self._hasValidDoseCache

    def _getCorrectedDose(self, acqOrder=None):
        """ Reimplement this because of a special tomo case. """
        _, preExp, dose = self._getDoseParams()
        # Acquisition metadata may not always carry dose values (e.g.
        # an import that didn't set them, or a movie whose acquisition
        # dose is not visible yet at processing time) - default to 0
        # instead of raising a TypeError on the arithmetic below. The
        # tomo branch already guarded against this; the non-tomo one
        # (the else branch) did not.
        preExp = preExp or 0.0
        dose = dose or 0.0

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
