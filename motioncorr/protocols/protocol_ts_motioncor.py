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

from tomo.protocols import ProtTsCorrectMotion

from .. import Plugin
from .protocol_base import ProtMotionCorrBase


class ProtTsMotionCorr(ProtMotionCorrBase, ProtTsCorrectMotion):
    """ This protocol wraps motioncor3 movie alignment program developed at UCSF.

    Motioncor3 performs anisotropic drift correction
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'align tilt-series movies'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtTsCorrectMotion.__init__(self, **kwargs)
        ProtMotionCorrBase.__init__(self, **kwargs)
        self.doSaveMovie = False
        self.doApplyDoseFilter = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, addEvenOddParam=True):
        ProtTsCorrectMotion._defineParams(self, form)
        self._defineCommonParams(form, allowDW=False)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self, inputId):
        if self.isEER:
            self._prepareEERFiles()
        else:
            ProtTsCorrectMotion.convertInputStep(self, inputId)

    def _processTiltImageM(self, workingFolder, tiltImageM, *args):
        outputFn, _ = self._getOutputTiltImagePaths(tiltImageM)

        def _getPath(path):
            """ shortcut to get relative path from workingFolder. """
            return os.path.relpath(path, workingFolder)

        self.info(f"workingFolder: {workingFolder}")
        self.info(f"outputFn: {outputFn}")

        argsDict = self._getMcArgs(acqOrder=tiltImageM.getAcquisitionOrder())
        argsDict['-OutMrc'] = f'"{_getPath(outputFn)}"'

        tiFn = tiltImageM.getFileName()
        inputFn = os.path.abspath(tiFn)

        self.info(f"inputFn: {tiFn}")

        params = self._getInputFormat(inputFn, absPath=True)
        params += ' '.join(['%s %s' % (k, v)
                            for k, v in argsDict.items()])

        params += ' ' + self.extraParams2.get()

        try:
            self.runJob(Plugin.getProgram(), params,
                        cwd=workingFolder,
                        env=Plugin.getEnviron())

            # Move output log to extra dir
            logFn = os.path.join(workingFolder, self._getMovieLogFile(tiltImageM))
            logFnExtra = self._getExtraPath(self._getMovieLogFile(tiltImageM))
            pwutils.moveFile(logFn, logFnExtra)

        except Exception as e:
            self.error(f"ERROR: Motioncor3 has failed for {tiFn} --> {str(e)}\n")
            import traceback
            traceback.print_exc()

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        tiltImageM = self._tsDict.getTi(tsId, tiltImageId)
        workingFolder = self._getTmpPath(self._getTiltImageMRoot(tiltImageM))
        pwutils.makePath(workingFolder)
        self._processTiltImageM(workingFolder, tiltImageM, *args)

        if self._doSplitEvenOdd():
            baseName = self._getTiltImageMRoot(tiltImageM)

            evenName = (os.path.abspath(self._getExtraPath(baseName + '_EVN.mrc')))
            oddName = (os.path.abspath(self._getExtraPath(baseName + '_ODD.mrc')))

            # Store the corresponding tsImM to use its data later in the even/odd TS
            tiltImageM.setOddEven([oddName, evenName])

        tiFn, _ = self._getOutputTiltImagePaths(tiltImageM)
        if not os.path.exists(tiFn):
            raise FileNotFoundError(f"Expected output file '{tiFn}' not produced!")

        if not pwutils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pwutils.cleanPath(workingFolder)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, 'outputTiltSeries'):
            summary.append('Aligned %d tilt series movies using Motioncor3.'
                           % self.getInputMovies().getSize())
        else:
            summary.append('Output is not ready')

        return summary

    def _validate(self):
        errors = ProtMotionCorrBase._validate(self)

        return errors

    # --------------------------- UTILS functions -----------------------------
    def getInputMovies(self):
        return self.inputTiltSeriesM.get()

    def _getMovieLogFile(self, tiltImageM):
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s%s-Full.log' % (pwutils.removeBaseExt(tiltImageM.getFileName()),
                                  '-Patch' if usePatches else '')

    def _createOutputWeightedTS(self):
        return False
