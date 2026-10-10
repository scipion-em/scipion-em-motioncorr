"""Regression tests for MotionCor3 classic protocol output paths.

Run via:
    scipion3 test motioncorr.tests.test_motioncorr_output_paths.TestMotioncorrOutputPaths
"""

import os
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from motioncorr import Plugin
from motioncorr.protocols.protocol_motioncorr import ProtMotionCorr


class _FakeMovie:
    def __init__(self, filename):
        self._filename = filename

    def getFileName(self):
        return self._filename

    def getAttributeValue(self, name, default=None):
        return default


class _NamingProtocol:
    patchX = 5
    patchY = 5
    _getMovieRoot = ProtMotionCorr._getMovieRoot
    _getOutputMicName = ProtMotionCorr._getOutputMicName
    _getMovieLogFile = ProtMotionCorr._getMovieLogFile


class TestMotioncorrOutputPaths(unittest.TestCase):
    def test_ClassicOutputNamesKeepTimestampAndMovieSuffix(self):
        protocol = _NamingProtocol()
        first = _FakeMovie("Runs/input/ja_115-1_0000_May14_21.14.32_1.tif")
        second = _FakeMovie("Runs/input/ja_115-1_0000_May14_21.14.32_2.tif")

        self.assertEqual(
            protocol._getMovieRoot(first), "ja_115-1_0000_May14_21.14.32_1"
        )
        self.assertEqual(
            protocol._getMovieRoot(second), "ja_115-1_0000_May14_21.14.32_2"
        )
        self.assertNotEqual(
            protocol._getOutputMicName(first),
            protocol._getOutputMicName(second),
        )

    def test_ClassicLogNameMatchesMotionCor3InputStem(self):
        protocol = _NamingProtocol()
        movie = _FakeMovie("Runs/input/ja_115-1_0000_May14_21.14.32_1.tif")
        self.assertEqual(
            protocol._getMovieLogFile(movie),
            "ja_115-1_0000_May14_21.14.32_1-Patch-Full.log",
        )

    def test_ClassicRunJobReceivesAbsoluteGainPath(self):
        relative_gain = "Runs/000027_XmippProtMovieGain/extra/orientedGain.mrc"

        class FakeProtocol:
            isEER = False
            doComputePSD = False
            extraParams2 = SimpleNamespace(get=lambda: "")
            defectFile = SimpleNamespace(get=lambda: "")
            defectMap = SimpleNamespace(get=lambda: "")

            def __init__(self):
                self.command_args = None

            def getInputMovies(self):
                return SimpleNamespace(getSamplingRate=lambda: 0.79)

            def _getOutputMovieFolder(self, movie):
                return "Runs/000025_ProtMotionCorr/tmp/movie_000002"

            def _getOutputMicName(self, movie):
                return "corrected.mrc"

            def _getMcArgs(self):
                return {"-Gain": '"' + relative_gain + '"'}

            def _getInputFormat(self, file_name):
                return ' -InTiff "movie.tif" '

            def runJob(self, program, args, **kwargs):
                self.command_args = args

            def _moveOutput(self, movie):
                pass

            def _useWorkerThread(self):
                return False

            def _saveAlignmentPlots(self, movie, sampling_rate):
                pass

            def _getExtraPath(self, path):
                return "Runs/000025_ProtMotionCorr/extra/" + path

            def _getMicFn(self, movie):
                return "corrected.mrc"

            def _doComputeMicThumbnail(self):
                return False

            def error(self, message):
                raise AssertionError(message)

        fake_protocol = FakeProtocol()
        movie = _FakeMovie("Runs/input/movie.tif")
        with patch.object(Plugin, "getProgram", return_value="MotionCor3"), \
                patch.object(Plugin, "getEnviron", return_value={}):
            ProtMotionCorr._processMovie(fake_protocol, movie)

        self.assertIsNotNone(fake_protocol.command_args)
        self.assertIn(
            '-Gain "' + os.path.abspath(relative_gain) + '"',
            fake_protocol.command_args,
        )
