# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
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
# **************************************************************************

from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr
from tomo.tests import DataSet  # initialization of tomo data set definition
from tomo.protocols import ProtImportTsMovies
from ..protocols import ProtTsMotionCorr


class TestMotioncor2TiltSeriesAlignMovies(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.getFileM = cls.dataset.getFile('empiar')

    def _runImportTiltSeriesM(self, filesPattern='{TS}_{TO}_{TA}.mrc'):
        protImport = self.newProtocol(ProtImportTsMovies,
                                      filesPath=self.getFileM,
                                      filesPattern=filesPattern,
                                      voltage=300,
                                      magnification=105000,
                                      sphericalAberration=2.7,
                                      amplitudeContrast=0.1,
                                      samplingRate=0.675,
                                      doseInitial=0,
                                      dosePerFrame=3.0,
                                      tiltAxisAngle=84.1)
        self.launchProtocol(protImport)
        return protImport

    def test_tiltseries_motioncor2(self):
        print(magentaStr("\n==> Importing data - TiltSeries:"))
        protImport = self._runImportTiltSeriesM()
        print(magentaStr("\n==> Testing motioncor2 - patch-based:"))
        protMc = self.newProtocol(ProtTsMotionCorr)
        protMc.inputTiltSeriesM.set(protImport.outputTiltSeriesM)
        self.launchProtocol(protMc)
        self.checkTSSet(protMc.outputTiltSeries, 2, 3, checkIds=True)

    def checkTSSet(self, set, size, anglesCount, checkIds=False):
        """
        Check basic attributes of a TS set
        :param set: TiltSeries set (Movies or Images)
        :param size: Expected size
        :param anglesCount: Expected number of tilts
        :param checkIds: check if ids start with 1 and increments by one
        :return: None
        """
        self.assertSetSize(set, size)
        for ts in set:
            self.assertEquals(ts.getSize(), anglesCount,
                              "Size of tilt images is wrong.")
            for i, ti in enumerate(ts):
                if checkIds:
                    self.assertEqual(i + 1, ti.getObjId(),
                                     "Tilt image Movie objId is incorrect")

                self.assertEqual(ts.getTsId(), ti.getTsId(),
                                 "Tilt image id is incorrect")
                self.assertEqual(ti.getSamplingRate(), ts.getSamplingRate(),
                                 'Tilt image sampling rate must be equal to '
                                 'tilt serie')
