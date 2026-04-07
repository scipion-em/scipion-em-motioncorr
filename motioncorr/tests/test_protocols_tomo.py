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
from pyworkflow.tests import setupTestProject
from pyworkflow.utils import magentaStr, cyanStr
from tomo.objects import SetOfTiltSeries, TomoAcquisition
from tomo.tests import DataSet
from tomo.protocols import ProtImportTsMovies
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from ..protocols import ProtTsMotionCorr


class TestMotioncorTiltSeriesM(TestBaseCentralizedLayer):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.getFileM = cls.dataset.getFile('empiar')
        cls.samplingRate = 0.675
        cls.testAcq = TomoAcquisition(
            angleMin=-60,
            angleMax=60,
            step=60,  # Only 3 angular stacks (-60, 0, 60) deg.
            tiltAxisAngle=84.1,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            doseInitial=0,
            dosePerFrame=3.0,
            accumDose=120.0
        )
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        cls._runPreviousProtocols()
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTsM = cls._runImportTiltSeriesM()

    @classmethod
    def _runImportTiltSeriesM(cls, filesPattern='{TS}_{TO}_{TA}.mrc'):
        print(magentaStr("\n==> Importing data - tilt-series movies:"))
        protImport = cls.newProtocol(ProtImportTsMovies,
                                     filesPath=cls.getFileM,
                                     filesPattern=filesPattern,
                                     voltage=cls.testAcq.getVoltage(),
                                     magnification=cls.testAcq.getMagnification(),
                                     sphericalAberration=cls.testAcq.getSphericalAberration(),
                                     amplitudeContrast=cls.testAcq.getAmplitudeContrast(),
                                     samplingRate=cls.samplingRate,
                                     doseInitial=cls.testAcq.getDoseInitial(),
                                     dosePerFrame=cls.testAcq.getDosePerFrame(),
                                     tiltAxisAngle=cls.testAcq.getTiltAxisAngle())
        cls.launchProtocol(protImport)
        return getattr(protImport, protImport.OUTPUT_NAME, None)

    def _runTiltSeriesMotionCorr(self, splitEvenOdd: bool = False):
        print(magentaStr("\n==> Running the tilt-series motion correction:"
                         f"\n\t- Split odd/even = {splitEvenOdd}"))
        protMc = self.newProtocol(ProtTsMotionCorr,
                                  inputTiltSeriesM=self.importedTsM,
                                  splitEvenOdd=splitEvenOdd)
        label = 'Motion Correction'
        label += ' OE' if splitEvenOdd else ''
        protMc.setObjLabel(label)
        self.launchProtocol(protMc)
        outTsSet = getattr(protMc, protMc._possibleOutputs.tiltSeries.name, None)
        self._checkTsSet(outTsSet, splitEvenOdd=splitEvenOdd)

    def _checkTsSet(self, outTsSet: SetOfTiltSeries, splitEvenOdd: bool = False):
        self.checkTiltSeries(
            outTsSet,
            expectedSetSize=2,
            expectedSRate=self.samplingRate,
            expectedDimensions=[7420, 7676, 3],
            testAcqObj=self.testAcq,
            imported=True,
            hasOddEven=splitEvenOdd,
            anglesCount=3,
            presentTsIds=['TS_01', 'TS_02']
        )

    def testTsMotionCorr_01(self):
        self._runTiltSeriesMotionCorr()

    def testTsMotionCorr_02(self):
        self._runTiltSeriesMotionCorr(splitEvenOdd=True)
