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

import os

import pyworkflow.utils as pwutils
from pyworkflow import BETA
import pyworkflow.object as pwobj

from tomo.protocols import ProtTsCorrectMotion

from .. import Plugin
from .protocol_base import ProtMotionCorrBase


class ProtTsMotionCorr(ProtMotionCorrBase, ProtTsCorrectMotion):
    """ This protocol wraps motioncor2 movie alignment program developed at UCSF.

    Motioncor2 performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'align tilt-series movies'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtTsCorrectMotion.__init__(self, **kwargs)
        ProtMotionCorrBase.__init__(self, **kwargs)
        self.doSaveMovie = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtTsCorrectMotion._defineParams(self, form)
        self._defineCommonParams(form)
        form.getParam('doApplyDoseFilter').default = pwobj.Boolean(False)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, inputId):
        self._prepareEERFiles()
        ProtTsCorrectMotion.convertInputStep(self, inputId)

    def _processTiltImageM(self, workingFolder, tiltImageM,
                           initialDose, dosePerFrame, gain, dark, *args):
        outputFn, outputFnDW = self._getOutputTiltImagePaths(tiltImageM)

        def _getPath(path):
            """ shortcut to get relative path from workingFolder. """
            return os.path.relpath(path, workingFolder)

        self.info(f"workingFolder: {workingFolder}")
        self.info(f"outputFn: {outputFn}")

        argsDict = self._getMcArgs()
        argsDict['-OutMrc'] = f'"{_getPath(outputFn)}"'

        numbOfFrames = self._getNumberOfFrames()
        dose = dosePerFrame / numbOfFrames if dosePerFrame else 0.0

        if self.isEER:
            argsDict['-EerSampling'] = self.eerSampling.get() + 1
            argsDict['-FmIntFile'] = "FmIntFile.txt"

            with open("FmIntFile.txt", "w") as f:
                f.write(f"{numbOfFrames} {self.eerGroup.get()} {dose}")

        elif self.doApplyDoseFilter:
            order = tiltImageM.getAcquisitionOrder()
            initDose = initialDose + (order - 1) * dosePerFrame
            argsDict.update({'-FmDose': dose,
                             '-InitDose': initDose if initDose > 0.001 else 0})

        tiFn = tiltImageM.getFileName()
        inputFn = os.path.abspath(tiFn)

        self.info(f"inputFn: {tiFn}")

        args = self._getInputFormat(inputFn, absPath=True)
        args += ' '.join(['%s %s' % (k, v)
                          for k, v in argsDict.items()])

        if gain:
            args += f' -Gain "{gain}" '
            args += f' -RotGain {self.gainRot.get()}'
            args += f' -FlipGain {self.gainFlip.get()}'

        if dark:
            args += f' -Dark "{dark}"'

        args += ' ' + self.extraParams2.get()

        try:
            self.runJob(Plugin.getProgram(), args,
                        cwd=workingFolder,
                        env=Plugin.getEnviron())

            # Move output log to extra dir
            logFn = os.path.join(workingFolder, self._getMovieLogFile(tiltImageM))
            logFnExtra = self._getExtraPath(self._getMovieLogFile(tiltImageM))
            pwutils.moveFile(logFn, logFnExtra)

            if self.doApplyDoseFilter and not self.doSaveUnweightedMic:
                pwutils.cleanPath(outputFn)

        except Exception as e:
            self.error(f"ERROR: Motioncor2 has failed for {tiFn} --> {str(e)}\n")
            import traceback
            traceback.print_exc()

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        tiltImageM = self._tsDict.getTi(tsId, tiltImageId)
        workingFolder = self._getTmpPath(self._getTiltImageMRoot(tiltImageM))
        pwutils.makePath(workingFolder)
        self._processTiltImageM(workingFolder, tiltImageM, *args)

        if self._doSplitEvenOdd():
            baseName = self._getTiltImageMRoot(tiltImageM)
            evenName = os.path.abspath(self._getExtraPath(baseName + '_EVN.mrc'))
            oddName = os.path.abspath(self._getExtraPath(baseName + '_ODD.mrc'))

            # Store the corresponding tsImM to use its data later in the even/odd TS
            self.tsMList.append(tiltImageM)

            # Update even and odd average lists
            self.evenAvgFrameList.append(evenName)
            self.oddAvgFrameList.append(oddName)

        tiFn, tiFnDW = self._getOutputTiltImagePaths(tiltImageM)
        if not os.path.exists(tiFn) and not os.path.exists(tiFnDW):
            raise Exception(f"Expected output file(s) '{tiFn}' not produced!")

        if not pwutils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pwutils.cleanPath(workingFolder)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'outputTiltSeries'):
            summary.append('Aligned %d tilt series movies using motioncor2.'
                           % self.getInputMovies().getSize())
            if self.splitEvenOdd and self._createOutputWeightedTS():
                summary.append('Even/odd outputs are dose-weighted!')
        else:
            summary.append('Output is not ready')

        return summary

    # --------------------------- UTILS functions -----------------------------
    def getInputMovies(self):
        return self.inputTiltSeriesM.get()

    def _getMovieLogFile(self, tiltImageM):
        usePatches = self.patchX != 0 or self.patchY != 0

        if Plugin.versionGE('1.6.2'):
            return '%s%s-Full.log' % (pwutils.removeBaseExt(tiltImageM.getFileName()),
                                      '-Patch' if usePatches else '')
        else:
            return '%s%s-Full.log' % (self._getTiltImageMRoot(tiltImageM),
                                      '-Patch' if usePatches else '')

    def _createOutputWeightedTS(self):
        return self.doApplyDoseFilter
