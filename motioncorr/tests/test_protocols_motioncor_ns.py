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
import tempfile

from pwem.protocols import ProtImportMovies
from pwem.objects import CTFModel
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr

from ..protocols import ProtMotionCorrTasks, ProtMotionCorrNewStreaming


class TestMotioncorNSAlignMovies(BaseTest):
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
            McProtClass = ProtMotionCorrNewStreaming
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
        output = protocol._possibleOutputs.micrographsDW.name
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
        prot.inputMovies.set(getattr(self.protImport1, self.protImport1._possibleOutputs.outputMovies.name, None))
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkGainFile(prot)
        self._checkAlignment(getattr(prot, prot._possibleOutputs.movies.name)[1],
                             (1, 38), [0, 0, 0, 0])

    def test_tif2(self):
        print(magentaStr("\n==> Testing motioncor - tif movies (2):"))
        prot = self.newProtocolMc(
                                objLabel='tif - motioncor (2)',
                                patchX=0, patchY=0)
        prot.inputMovies.set(getattr(self.protImport3, self.protImport3._possibleOutputs.outputMovies.name, None))
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkGainFile(prot)
        self._checkAlignment(getattr(prot, prot._possibleOutputs.movies.name)[1],
                             (1, 24), [0, 0, 0, 0])

    def test_mrcs(self):
        print(magentaStr("\n==> Testing motioncor - mrcs movies:"))
        prot = self.newProtocolMc(
                                objLabel='mrcs - motioncor',
                                patchX=0, patchY=0)
        prot.inputMovies.set(getattr(self.protImport2, self.protImport2._possibleOutputs.outputMovies.name, None))
        self.launchProtocol(prot)

        self._checkOutput(prot)
        self._checkAlignment(getattr(prot, prot._possibleOutputs.movies.name)[1],
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


class TestMotioncorNSReadCtfModel(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _assertStandardizedAngle(self, ctfAngle, rawAngle, places=2):
        """setStandardDefocus normalizes astigmatism angle to [0, 180)."""
        expected = rawAngle % 180.0
        self.assertAlmostEqual(ctfAngle, expected, places=places)

    def test_readCtfModel_motioncor_txt(self):
        prot = self.newProtocol(ProtMotionCorrNewStreaming)

        ctf_content = (
            "# Columns: #1 micrograph number; #2 - defocus 1 [A]; #3 - defocus 2; #4 - azimuth\n"
            "# of astigmatism; #5 - additional phase shift [radian]; #6 - cross correlation;\n"
            "#7 - spacing (in Angstroms) up to which CTF rings were fit successfully\n"
            "  11744.11   11479.96   -50.57     0.10    0.06715   17.2973\n"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            ctf_fn = os.path.join(tmpdir, "movie_aligned_mic_Ctf.txt")
            psd_fn = os.path.join(tmpdir, "movie_aligned_mic_Ctf.mrc")

            with open(ctf_fn, "w") as fh:
                fh.write(ctf_content)

            with open(psd_fn, "wb") as fh:
                fh.write(b"\0")

            ctf = prot._readCtfModel(CTFModel(), ctf_fn, psd_fn)

            self.assertAlmostEqual(ctf.getDefocusU(), 11744.11, places=2)
            self.assertAlmostEqual(ctf.getDefocusV(), 11479.96, places=2)
            self._assertStandardizedAngle(ctf.getDefocusAngle(), -50.57)
            self.assertAlmostEqual(ctf.getFitQuality(), 0.06715, places=5)
            self.assertAlmostEqual(ctf.getResolution(), 17.2973, places=4)
            self.assertAlmostEqual(ctf.getPhaseShift(), 5.72957795, places=5)
            self.assertEqual(ctf.getPsdFile(), psd_fn)

    def test_readCtfModel_with_leading_index(self):
        prot = self.newProtocol(ProtMotionCorrNewStreaming)

        # Variant with leading micrograph index:
        # [idx, defocusU, defocusV, defocusAngle, phaseShiftRad, fit, resolution]
        ctf_content = "1 12153.57 11175.66 -52.43 0.00 0.09779 13.0612\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            ctf_fn = os.path.join(tmpdir, "movie_aligned_mic_Ctf.txt")
            with open(ctf_fn, "w") as fh:
                fh.write(ctf_content)

            ctf = prot._readCtfModel(CTFModel(), ctf_fn)

            self.assertAlmostEqual(ctf.getDefocusU(), 12153.57, places=2)
            self.assertAlmostEqual(ctf.getDefocusV(), 11175.66, places=2)
            self._assertStandardizedAngle(ctf.getDefocusAngle(), -52.43)
            self.assertAlmostEqual(ctf.getFitQuality(), 0.09779, places=5)
            self.assertAlmostEqual(ctf.getResolution(), 13.0612, places=4)
