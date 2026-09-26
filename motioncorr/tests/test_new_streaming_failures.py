from unittest import TestCase
from unittest.mock import patch

from motioncorr.protocols.protocol_motioncorr_ns import ProtMotionCorrNewStreaming


class _Value:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _NonEmptyOutput:
    def __len__(self):
        return 1


class _NewStreamingHarness(ProtMotionCorrNewStreaming):
    def __init__(self):
        self.failedMovies = []
        self.sRate = 1.0
        self.extraParams2 = _Value("")
        self.runCalls = []
        self.closed = False
        self.warnings = []
        self.publishedFailedMovies = None
        self.movies = _NonEmptyOutput()

    def _getMcArgs(self):
        return {}

    def _getResultMicFn(self, movieFName, suffix=""):
        return "/tmp/output.mrc"

    def _getExtraPath(self, *paths):
        return "/tmp" if not paths else "/tmp/" + "/".join(paths)

    def _getInputFormat(self, movieFName, absPath=False):
        return ""

    def runJob(self, *args, **kwargs):
        self.runCalls.append(args)
        if len(self.runCalls) == 1:
            raise RuntimeError("motioncor failed")

    def _saveAlignmentPlots(self, *args, **kwargs):
        pass

    def _closeOutputSet(self):
        self.closed = True

    def warning(self, message):
        self.warnings.append(message)

    def _createFailedMoviesOutput(self):
        self.publishedFailedMovies = list(self.failedMovies)


class TestMotionCorrNewStreamingFailures(TestCase):
    def test_FailedMovieDoesNotAbortFollowingMovieAndFinishesWithWarning(self):
        protocol = _NewStreamingHarness()

        module = "motioncorr.protocols.protocol_motioncorr_ns"
        with patch(module + ".Plugin.getProgram", return_value="MotionCor"), \
             patch(module + ".Plugin.getEnviron", return_value={}):
            protocol.processMovieStep("/data/movie_001.mrcs")
            protocol.processMovieStep("/data/movie_002.mrcs")

        self.assertEqual(
            2,
            len(protocol.runCalls),
            "A failed movie must not prevent the following movie from running.",
        )
        self.assertIn("/data/movie_001.mrcs", protocol.failedMovies)

        protocol.closeOutputSetStep(["movies"])

        self.assertTrue(
            protocol.closed,
            "Partial movie failures must not make the whole protocol fail.",
        )
        self.assertEqual(
            ["/data/movie_001.mrcs"],
            protocol.publishedFailedMovies,
            "Failed movies must be exposed explicitly as an output.",
        )
        self.assertTrue(
            any("failed" in message.lower() for message in protocol.warnings),
            "Partial failures must be reported with a warning.",
        )

    def test_AllMoviesFailedStillFailsProtocol(self):
        protocol = _NewStreamingHarness()
        protocol.failedMovies = [
            "/data/movie_001.mrcs",
            "/data/movie_002.mrcs",
        ]
        del protocol.movies

        with self.assertRaisesRegex(RuntimeError, "failed.*movie|movie.*failed"):
            protocol.closeOutputSetStep(["movies"])

        self.assertFalse(
            protocol.closed,
            "A run with no useful output must not close as a successful protocol.",
        )
