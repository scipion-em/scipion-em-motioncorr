# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es) [1]
# *             Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [2]
# *             Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [3]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Department of Anatomy and Cell Biology, McGill University
# * [3] MRC Laboratory of Molecular Biology (MRC-LMB)
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
import os.path

from pwem.protocols import ProtImportMovies
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr

from ..protocols import ProtMotionCorr, ProtMotionCorrTasks


class TestMotioncorAlignMovies(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds1 = DataSet.getDataSet('movies')
        cls.ds2 = DataSet.getDataSet('relion30_tutorial')

    @classmethod
    def runImportMovies(cls, pattern, label, **kwargs):
        """ Run an Import micrograph protocol. """
        protImport = cls.newProtocol(ProtImportMovies, filesPattern=pattern,
                                     **kwargs)
        protImport.setObjLabel(f"import movies - {label}")
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def newProtocolMc(cls, *args, **kwargs):
        if int(os.environ.get('TEST_MC_TASKS', 0)):
            McProtClass = ProtMotionCorrTasks
        else:
            McProtClass = ProtMotionCorr
        return cls.newProtocol(McProtClass, *args, **kwargs)

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        print(magentaStr("\n==> Importing data - movies:"))
        cls.protImport1 = cls.runImportMovies(
            cls.ds1.getFile('c3-adp-se-xyz-0228_200.tif'),
            "tif + dm4 gain",
            samplingRate=0.554,
            voltage=300,
            sphericalAberration=2.7,
            dosePerFrame=1.3,
            gainFile=cls.ds1.getFile('SuperRef_c3-adp-se-xyz-0228_001.dm4')
        )

        cls.protImport2 = cls.runImportMovies(
            cls.ds1.getFile('Falcon*.mrcs'),
            "mrcs",
            samplingRate=1.1,
            voltage=300,
            sphericalAberration=2.7,
            dosePerFrame=1.2
        )

        cls.protImport3 = cls.runImportMovies(
            cls.ds2.getFile('Movies/20170629_00030_frameImage.tiff'),
            "tif + mrc gain",
            samplingRate=0.885,
            voltage=200,
            sphericalAberration=1.4,
            dosePerFrame=1.277,
            gainFile=cls.ds2.getFile('Movies/gain.mrc')
        )

    def _checkOutput(self, protocol):
        output = "outputMicrographsDoseWeighted"
        self.assertIsNotNone(getattr(protocol, output),
                             "Output SetOfMicrographs was not created.")

    def _checkGainFile(self, protocol):
        gainFile = protocol.inputMovies.get().getGain()
        self.assertTrue(os.path.exists(gainFile),
                        f"Gain file {gainFile} was not found.")

    def _checkAlignment(self, movie, goldRange, goldRoi):
        alignment = movie.getAlignment()
        rangeFrames = alignment.getRange()
        aliFrames = rangeFrames[1] - rangeFrames[0] + 1
        msgRange = "Alignment range must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRange, rangeFrames, msgRange % (goldRange,
                                                             type(goldRange),
                                                             rangeFrames,
                                                             type(rangeFrames)))
        roi = alignment.getRoi()
        shifts = alignment.getShifts()
        zeroShifts = (aliFrames * [0], aliFrames * [0])
        nrShifts = len(shifts[0])
        msgRoi = "Alignment ROI must be %s (%s) and it is %s (%s)"
        msgShifts = "Alignment SHIFTS must be non-zero!"
        self.assertEqual(goldRoi, roi, msgRoi % (goldRoi, type(goldRoi),
                                                 roi, type(roi)))
        self.assertNotEqual(zeroShifts, shifts, msgShifts)
        self.assertEqual(nrShifts, aliFrames, "Number of shifts is not equal"
                                              "to number of aligned frames.")

    def test_tif(self):
        print(magentaStr("\n==> Testing motioncor - tif movies:"))
        prot = self.newProtocolMc(
                                objLabel='tif - motioncor',
                                patchX=0, patchY=0, binFactor=2,
                                gainFlip=1) # flip upside down because of dm4
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkGainFile(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 38), [0, 0, 0, 0])

    def test_tif2(self):
        print(magentaStr("\n==> Testing motioncor - tif movies (2):"))
        prot = self.newProtocolMc(
                                objLabel='tif - motioncor (2)',
                                patchX=0, patchY=0)
        prot.inputMovies.set(self.protImport3.outputMovies)
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkGainFile(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 24), [0, 0, 0, 0])

    def test_mrcs(self):
        print(magentaStr("\n==> Testing motioncor - mrcs movies:"))
        prot = self.newProtocolMc(
                                objLabel='mrcs - motioncor',
                                patchX=0, patchY=0)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 16), [0, 0, 0, 0])

    def test_eer(self):
        print(magentaStr("\n==> Testing motioncor - eer movies:"))
        protImport = self.runImportMovies(
            self.ds1.getFile('FoilHole*.eer'),
            "eer + gain",
            samplingRate=1.2,
            voltage=300,
            sphericalAberration=2.7,
            dosePerFrame=0.07,
            gainFile=self.ds1.getFile('eer.gain')
        )
        prot = self.newProtocolMc(
                                objLabel='eer - motioncor',
                                patchX=0, patchY=0, eerGroup=14)
        prot.inputMovies.set(protImport.outputMovies)
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkGainFile(prot)
