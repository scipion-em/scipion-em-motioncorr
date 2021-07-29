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

from pwem.protocols import ProtImportMovies
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr

from ..protocols import ProtMotionCorr


class TestMotioncor2AlignMovies(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern,
                  'dosePerFrame': 1.3
                  }
        if 'samplingRate' not in kwargs:
            del params['samplingRate']
            params['samplingRateMode'] = 0
        else:
            params['samplingRateMode'] = 1

        params.update(kwargs)

        protImport = cls.newProtocol(ProtImportMovies, **params)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        print(magentaStr("\n==> Importing data - movies:"))
        cls.protImport1 = cls.runImportMovies(cls.ds.getFile('qbeta/qbeta.mrc'),
                                              magnification=50000)
        cls.protImport2 = cls.runImportMovies(cls.ds.getFile('cct/cct_1.em'),
                                              magnification=61000)

    def _checkMicrographs(self, protocol):
        self.assertIsNotNone(getattr(protocol, 'outputMicrographsDoseWeighted'),
                             "Output SetOfMicrographs were not created.")

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
                                              " number of aligned frames.")

    def test_cct_motioncor2_patch(self):
        print(magentaStr("\n==> Testing motioncor2 - patch-based:"))
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='cct - motioncor2 test1 (patch-based)',
                                patchX=2, patchY=2)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0])

    def test_qbeta_motioncor2_patch(self):
        print(magentaStr("\n==> Testing motioncor2 - patch-based with grouping:"))
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncor2 test2 (patch-based)',
                                patchX=2, patchY=2,
                                group=2)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0])

    def test_qbeta_motioncor2_sel(self):
        print(magentaStr("\n==> Testing motioncor2 - patch-based with frame range:"))
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncor2 test3 (frame range)',
                                patchX=2, patchY=2,
                                alignFrame0=2,
                                alignFrameN=6)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (2, 6), [0, 0, 0, 0])
