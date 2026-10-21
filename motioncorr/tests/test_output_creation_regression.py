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
"""Creating the streaming outputs must work when the Set kept by _defineOutputs is the one just created.

A freshly created, still empty Set has stored no properties yet, so asking it to reload them fails.
"""
import tempfile
import unittest

from pwem.objects import SetOfMovies, Acquisition

from motioncorr.protocols.protocol_motioncorr_ns import ProtMotionCorrNewStreaming, MotionCorrOutputs


class _ProtocolStub:
    """Keeps the very Set it is given as output, as a run whose outputs are the created Sets."""
    _possibleOutputs = MotionCorrOutputs

    def __init__(self, path):
        self.sRate = 1.5
        self._path = path
        self.inputMovies = SetOfMovies(filename=path + '/input_movies.db')
        self.inputMovies.setSamplingRate(1.0)
        self.inputMovies.setAcquisition(Acquisition(voltage=300, sphericalAberration=2.7, amplitudeContrast=0.1))

    def getInputMovies(self, asPointer=False):
        return self.inputMovies

    def _getPath(self, *paths):
        return self._path

    def _defineOutputs(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self, name, value)

    def _defineSourceRelation(self, source, output):
        pass


class TestNewEmptyOutputsCanBeCreated(unittest.TestCase):

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.protocol = _ProtocolStub(tmp.name)

    def testNewMoviesOutput(self):
        outputMovies = ProtMotionCorrNewStreaming._getOutputMovies(self.protocol)

        self.assertEqual(outputMovies.getSamplingRate(), 1.5)
        self.assertEqual(outputMovies.getAcquisition().getVoltage(), 300)

    def testNewMicrographsOutput(self):
        name = MotionCorrOutputs.micrographs.name
        outputMics = ProtMotionCorrNewStreaming._getOutputMics(self.protocol, name)

        self.assertIs(outputMics, getattr(self.protocol, name))
        self.assertEqual(outputMics.getSamplingRate(), 1.5)


if __name__ == '__main__':
    unittest.main()
