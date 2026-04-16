# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import traceback
from enum import Enum
from os.path import basename, join, abspath
from pathlib import Path
from typing import List
import mrcfile
import numpy as np
from pyworkflow import BETA
from pyworkflow.object import Set, Pointer
from pyworkflow.protocol import PointerParam, IntParam, BooleanParam, FloatParam, LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.utils import cyanStr, makePath, Message, cleanPath, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltImageM, TiltImage, TiltSeries, SetOfTiltSeries, TiltSeriesM, SetOfTiltSeriesM
from .. import Plugin
from .protocol_base import ProtMotionCorrBase


logger = logging.getLogger(__name__)
MRCS_EXT = '.mrcs'
MRC_EXT = '.mrc'
EVEN_SUFFIX = '_EVN'
ODD_SUFFIX = '_ODD'
OUTPUT_TSM_FAILED_NAME = "FailedTiltSeriesMovies"


class TSMcorrOutputs(Enum):
    tiltSeries = SetOfTiltSeries

class ProtTsMotionCorr(ProtMotionCorrBase):
    """ This protocol wraps motioncor movie alignment program developed at UCSF.

    Motioncor performs anisotropic drift correction
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'align tilt-series movies'
    _possibleOutputs = TSMcorrOutputs
    _devStatus = BETA
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.doSaveMovie = False
        self.doApplyDoseFilter = False
        self.tsDict = None
        self.patchStr = None
        self.sRate = -1
        self.failedItems = []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
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

        form.addParam('splitEvenOdd', BooleanParam,
                      default=False,
                      label='Split & sum odd/even frames?',
                      help='(Used for denoising data preparation). If set to Yes, 2 additional movies/tilt '
                           'series will be generated, one generated from the even frames and the other from the '
                           'odd ones using the same alignment for the whole stack of frames.')

        form.addParam('removeIndivImgs', BooleanParam,
                      label='Remove the generated unstacked images?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to True, the individual unstacked images will be removed before the '
                           'protocol execution ends. Only the staked tilt-series will remain.')

        self._defineCommonParams(form, allowDW=False)
        # Patch alignment is not recommended
        form.getParam('patchX').setDefault(0)
        form.getParam('patchY').setDefault(0)
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        cInPId = self._insertFunctionStep(self._convertInputStep,
                                          prerequisites=[],
                                          needsGPU=False)
        for tsId, tsM in self.tsMDict.items():
            makePath(self._getTsResultsPath(tsId))
            pasPId = None
            for counter, tiM in enumerate(tsM.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD)):
                tiM = tiM.clone()
                pasPId = self._insertFunctionStep(self.processAngularStackStep,
                                         tsId,
                                         tiM,
                                         counter + 1,
                                         prerequisites=cInPId,
                                         needsGPU=True)
            if not pasPId:
                logger.error(redStr(f'tsId = {tsId}: no tilt-image movies were found.'))
                continue
            cOutPId = self._insertFunctionStep(self.createOutputStep,
                                               tsId,
                                               tsM,
                                               prerequisites=pasPId,
                                               needsGPU=False)
            closeSetStepDeps.append(cOutPId)
        self._insertFunctionStep(self._closeOutputSet,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.tsMDict = {tsM.getTsId(): tsM.clone() for tsM in self.getInputMovies().iterItems()}
        usePatches = self.patchX != 0 or self.patchY != 0
        self.patchStr = '-Patch' if usePatches else ''
        self.sRate = self.getInputMovies().getSamplingRate()

    def processAngularStackStep(self, tsId: str, tiM: TiltImageM, counter: int):
        tiMFileName = tiM.getFileName()
        logger.info(cyanStr(f'tsId = {tsId} - processing angular stack {basename(tiMFileName)}'))
        try:
            argsDict = self._getMcArgs(acqOrder=tiM.getAcquisitionOrder())
            argsDict['-OutMrc'] = f'"{abspath(self._getOutTsFName(tsId, counter))}"'
            argsDict['-LogDir'] = f'"{abspath(self._getTsResultsPath(tsId))}"'
            params = self._getInputFormat(tiMFileName, absPath=True)
            params += ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
            params += ' ' + self.extraParams2.get()
            self.runJob(self.program, params, env=Plugin.getEnviron())
        except Exception as e:
            self.error(f"ERROR: Motioncor has failed for {tiMFileName} --> {str(e)}\n")
            traceback.print_exc()

    def createOutputStep(self, tsId: str, tsM: TiltSeriesM):
        if tsId in self.failedItems:
            self.createOutputFailedSet(tsM)
            return
        logger.info(cyanStr(f'===> tsId = {tsId}: Creating the resulting tilt-series...'))
        outStackFn = self._mountFinalStack(tsId)
        outStackFnOdd, outStackFnEven = '', ''
        if self.splitEvenOdd.get():
            outStackFnOdd = self._mountFinalStack(tsId, suffix=ODD_SUFFIX)
            outStackFnEven = self._mountFinalStack(tsId, suffix=EVEN_SUFFIX)
        self._registerOutput(tsM, outStackFn, outStackFnEven, outStackFnOdd)

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
    def getInputMovies(self, asPointer=False):
        return self.inputTiltSeriesM if asPointer else self.inputTiltSeriesM.get()

    def _getTsResultsPath(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def _getOutTsFName(self, tsId: str, counter: int) -> str:
        return join(self._getTsResultsPath(tsId), f'{tsId}_{counter:03d}.mrc')

    def _getOutStackFName(self, tsId: str,  suffix: str = ''):
        return join(self._getTsResultsPath(tsId), f'{tsId}{suffix}{MRCS_EXT}')

    def _mountFinalStack(self, tsId: str, suffix: str = '') -> str:
        logger.info(cyanStr(f'===> tsId = {tsId}{suffix}: mounting the stack file...'))
        outStackFile = self._getOutStackFName(tsId, suffix=suffix)
        resultImgs = self._getResultImgs(tsId, suffix=suffix)
        # Read the first image to get the dimensions
        with mrcfile.mmap(resultImgs[0], mode='r+') as mrc:
            img = mrc.data
            nx, ny = img.shape

        # Create an empty array in which the stack of images will be stored
        shape = (len(resultImgs), nx, ny)
        stackArray = np.empty(shape, dtype=img.dtype)

        # Fill it with the images sorted by angle
        for i, img in enumerate(resultImgs):
            with mrcfile.mmap(img) as mrc:
                logger.info(f'Inserting image - index [{i}], {img}')
                stackArray[i] = mrc.data

        # Save the stack in a new mrc file
        with mrcfile.new_mmap(outStackFile, shape, overwrite=True) as mrc:
            mrc.set_data(stackArray)
            mrc.update_header_from_data()
            mrc.update_header_stats()
            mrc.voxel_size = self.sRate

        # Remove the individual unstacked images if requested
        if self.removeIndivImgs.get():
            cleanPath(*resultImgs)
        return outStackFile

    def _getResultImgs(self, tsId: str, suffix: str = '') -> List[str]:
        imagesDir = Path(self._getTsResultsPath(tsId))
        pattern = f'*{suffix}{MRC_EXT}'
        exclusionWords = [EVEN_SUFFIX, ODD_SUFFIX] if suffix == '' else []
        if exclusionWords:
            finalList = [str(p) for p in imagesDir.glob(pattern) if
                         not any(exclusionWord in p.name for exclusionWord in exclusionWords)]
        else:
            finalList = [str(p) for p in imagesDir.glob(pattern)]
        return sorted(finalList)

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self,
                        inTsM: TiltSeriesM,
                        tsFName: str,
                        tsFnameEven: str,
                        tsFnameOdd: str):
        with self._lock:
            # Mount the resulting tilt-series
            outTsSet = self._getOutputTsSet()
            binFactor = self.binFactor.get()
            newTs = TiltSeries()
            newTs.copyInfo(inTsM)
            newTs.setSamplingRate(inTsM.getSamplingRate() * binFactor)
            outTsSet.append(newTs)

            for index, inTi in enumerate(inTsM.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD)):
                newTi = TiltImage()
                newTi.copyInfo(inTi)
                newTi.setSamplingRate(inTi.getSamplingRate()*binFactor)
                newTi.setAcquisition(inTi.getAcquisition())
                newTi.setFileName(tsFName)
                newTi.setIndex(index + 1)
                if self.splitEvenOdd.get():
                    newTi.setOddEven([tsFnameOdd, tsFnameEven])
                newTs.append(newTi)

            newTs.write()
            outTsSet.update(newTs)
            outTsSet.write()
            self._store(outTsSet)
            for outputName in self._possibleOutputs:
                output = getattr(self, outputName.name, None)
                if output:
                    output.close()

    def _getOutputTsSet(self) -> SetOfTiltSeries:
        outSetSetAttrib = self._possibleOutputs.tiltSeries.name
        outTsSet = getattr(self, outSetSetAttrib, None)
        if outTsSet:
            outTsSet.enableAppend()
        else:
            outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
            inM = self.getInputMovies()
            outTsSet.copyInfo(inM)
            outTsSet.setSamplingRate(inM.getSamplingRate()*self.binFactor.get())
            outTsSet.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outSetSetAttrib: outTsSet})
            self._defineSourceRelation(self.getInputMovies(asPointer=True), outTsSet)
        return outTsSet

    @retry_on_sqlite_lock(log=logger)
    def createOutputFailedSet(self, item: TiltSeriesM):
        """ Just copy input item to the failed output set. """
        with self._lock:
            logger.info(f'Failed TS ---> {item.getTsId()}')
            inputSetPointer = self.getInputMovies(asPointer=True)
            output = self._getOutputFailedSet(inputSetPointer)
            newItem = item.clone()
            newItem.copyInfo(item)
            output.append(newItem)

            if isinstance(item, TiltSeries):
                newItem.copyItems(item)
                newItem.write(properties=False)

            output.update(newItem)
            output.write()
            self._store(output)

            # Close explicitly the outputs (for streaming)
            output.close()

    def _getOutputFailedSet(self, inputPtr: Pointer):
        """ Create output set for failed TSM. """
        inputSet = inputPtr.get()
        failedTs = getattr(self, OUTPUT_TSM_FAILED_NAME, None)
        if failedTs:
            failedTs.enableAppend()
        else:
            logger.info(cyanStr('Create the set of failed TSM'))
            failedTs = SetOfTiltSeriesM.create(self._getPath(), template='tiltseriesM', suffix='Failed')
            failedTs.copyInfo(inputSet)
            failedTs.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{OUTPUT_TSM_FAILED_NAME: failedTs})
            self._defineSourceRelation(inputPtr, failedTs)

        return failedTs


