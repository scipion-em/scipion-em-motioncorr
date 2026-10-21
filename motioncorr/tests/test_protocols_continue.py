# **************************************************************************
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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""Continue aligns again what failed, on real runs.

A movie that cannot be read fails alone and the run still ends well. Once
it is readable again, Continue aligns it and leaves the others alone.
"""

import os

from pwem.protocols import ProtImportMovies
from pyworkflow.protocol.constants import MODE_RESUME
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr

from ..protocols import ProtMotionCorr

MOVIE_NAMES = ('20170629_00021_frameImage.tiff',
               '20170629_00022_frameImage.tiff',
               '20170629_00023_frameImage.tiff')


class TestMotioncorContinueAfterFailedMovies(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        ds = DataSet.getDataSet('relion30_tutorial')
        cls.moviesSource = ds.getFile('Movies/')
        cls.moviesDir = os.path.abspath(cls.proj.getTmpPath('movies'))
        os.makedirs(cls.moviesDir)

        # Linked, not copied: they are large, and a broken one is replaced
        # by a file of its own.
        for name in MOVIE_NAMES:
            os.symlink(os.path.join(cls.moviesSource, name),
                       os.path.join(cls.moviesDir, name))

        cls.protImport = cls.newProtocol(ProtImportMovies,
                                         filesPath=cls.moviesDir,
                                         filesPattern='*.tiff',
                                         samplingRate=0.885, voltage=200,
                                         sphericalAberration=1.4,
                                         dosePerFrame=1.277,
                                         gainFile=ds.getFile('Movies/gain.mrc'))
        cls.launchProtocol(cls.protImport)

    def _sizes(self, prot):
        return (prot.outputMovies.getSize(),
                prot.outputMicrographsDoseWeighted.getSize())

    def testAFailedMovieIsAlignedOnContinue(self):
        broken = os.path.join(self.moviesDir, MOVIE_NAMES[1])
        os.remove(broken)

        with open(broken, 'wb') as handle:
            handle.write(b'not a movie')

        print(magentaStr("\n==> Testing motioncor with a broken movie:"))
        prot = self.newProtocol(ProtMotionCorr, objLabel='mc broken movie',
                                patchX=0, patchY=0)
        prot.inputMovies.set(self.protImport.outputMovies)
        self.launchProtocol(prot)

        self.assertEqual((2, 2), self._sizes(prot),
                         "Only the broken movie is left out, of every output.")

        os.remove(broken)
        os.symlink(os.path.join(self.moviesSource, MOVIE_NAMES[1]), broken)
        print(magentaStr("\n==> Continue once the movie is readable:"))
        prot.runMode.set(MODE_RESUME)
        self.launchProtocol(prot)

        self.assertEqual((3, 3), self._sizes(prot),
                         "Continue aligns the movie that failed, once.")
