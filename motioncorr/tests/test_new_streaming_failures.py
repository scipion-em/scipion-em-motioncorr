import sqlite3
import threading
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


class _OutputItem:
    def __init__(self, obj_id):
        self._obj_id = obj_id

    def getObjId(self):
        return self._obj_id


class _OutputSet:
    def __init__(self, obj_ids):
        self._items = [_OutputItem(obj_id) for obj_id in obj_ids]

    def __iter__(self):
        return iter(self._items)

    def __len__(self):
        return len(self._items)


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

    def info(self, message):
        pass

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

    def test_ReadingOutputRequiresAllExpectedOutputsForResumeCheckpoint(self):
        protocol = _NewStreamingHarness()
        protocol.itemIdReadList = []
        protocol.splitEvenOdd = _Value(False)
        protocol.doApplyDoseFilter = _Value(False)

        class _InputSet:
            def getUniqueValues(self, field):
                if field != 'id':
                    raise AssertionError(f'Unexpected field: {field}')
                return [1, 2]

        protocol.getInputMovies = lambda *args, **kwargs: _InputSet()
        protocol.movies = _OutputSet([1, 2])
        protocol.micrographs = _OutputSet([1])
        protocol.info = lambda *args, **kwargs: None

        protocol.readingOutput()

        self.assertEqual(
            [1],
            protocol.itemIdReadList,
            "Resume must only consider a movie completed when every required "
            "persisted output contains that movie.",
        )

    def test_CreateOutputRetriesSqliteLockInsteadOfSwallowingIt(self):
        protocol = _NewStreamingHarness()
        protocol.doSaveMovie = _Value(False)
        protocol._lock = threading.RLock()

        attempts = []

        def locked_output_movies():
            attempts.append(1)
            raise sqlite3.OperationalError("database is locked")

        protocol._getOutputMovies = locked_output_movies

        with patch(
            "pyworkflow.utils.retry_streaming.time.sleep",
            return_value=None,
        ):
            with self.assertRaisesRegex(
                sqlite3.OperationalError,
                "database is locked",
            ):
                protocol.createOutputStep(
                    "/data/movie_001.mrcs",
                    None,
                )

        self.assertEqual(
            15,
            len(attempts),
            "createOutputStep must let SQLite lock errors reach "
            "retry_on_sqlite_lock so the decorator can retry them.",
        )

    def test_CreateOutputDoesNotSwallowProtocolStoreFailure(self):
        protocol = _NewStreamingHarness()
        protocol.doSaveMovie = _Value(False)
        protocol.doApplyDoseFilter = _Value(False)
        protocol.splitEvenOdd = _Value(False)
        protocol._lock = threading.RLock()

        class _InputMovie:
            def getObjId(self):
                return 7

        class _FakeMovie:
            def copyInfo(self, _):
                pass

            def setFileName(self, _):
                pass

            def setMicName(self, _):
                pass

            def getNumberOfFrames(self):
                return 1

            def setAlignment(self, _):
                pass

        class _OutputMovies:
            def __init__(self):
                self.writes = 0

            def __len__(self):
                return 0

            def __contains__(self, obj_id):
                return False

            def append(self, _):
                pass

            def update(self, _):
                pass

            def write(self):
                self.writes += 1

        outputMovies = _OutputMovies()
        protocol._getOutputMovies = lambda: outputMovies
        protocol.getMovieAlignment = lambda *args, **kwargs: object()

        def fail_protocol_store(*args, **kwargs):
            raise RuntimeError("protocol store failed")

        protocol._store = fail_protocol_store

        module = "motioncorr.protocols.protocol_motioncorr_ns"
        with patch(module + ".Movie", _FakeMovie):
            with self.assertRaisesRegex(
                RuntimeError,
                "protocol store failed",
            ):
                protocol.createOutputStep(
                    "/data/movie_001.mrcs",
                    _InputMovie(),
                )

        self.assertEqual(
            1,
            outputMovies.writes,
            "The test must reach the persisted-output checkpoint before "
            "simulating the protocol store failure.",
        )

    def test_RegisterMicsIsIdempotentForAlreadyPersistedMovie(self):
        protocol = _NewStreamingHarness()
        protocol.splitEvenOdd = _Value(False)
        protocol.doApplyDoseFilter = _Value(False)

        class _InputMovie:
            def getObjId(self):
                return 7

        class _FakeMicrograph:
            def __init__(self):
                self._obj_id = None

            def copyInfo(self, movie):
                self._obj_id = movie.getObjId()

            def getObjId(self):
                return self._obj_id

            def setFileName(self, _):
                pass

            def setSamplingRate(self, _):
                pass

        class _OutputMics:
            def __init__(self):
                self.appended_ids = []

            def __len__(self):
                return len(self.appended_ids)

            def __contains__(self, obj_id):
                return obj_id in self.appended_ids

            def append(self, item):
                self.appended_ids.append(item.getObjId())

            def update(self, _):
                pass

            def write(self):
                pass

        outputMics = _OutputMics()
        protocol._getOutputMics = lambda *args, **kwargs: outputMics
        protocol._store = lambda *args, **kwargs: None
        protocol.setMicPlotInfo = lambda *args, **kwargs: None
        protocol.setMicsEvenOdd = lambda *args, **kwargs: None

        module = "motioncorr.protocols.protocol_motioncorr_ns"
        with patch(module + ".Micrograph", _FakeMicrograph), \
             patch(module + ".setMRCSamplingRate", return_value=None):
            protocol._registerMics(
                "/data/movie_007.mrcs",
                _InputMovie(),
                "micrographs",
            )
            protocol._registerMics(
                "/data/movie_007.mrcs",
                _InputMovie(),
                "micrographs",
            )

        self.assertEqual(
            [7],
            outputMics.appended_ids,
            "Retrying a partially persisted movie must not append the same "
            "micrograph twice.",
        )

    def test_CreateOutputDoesNotDuplicateAlreadyPersistedMovie(self):
        protocol = _NewStreamingHarness()
        protocol.doSaveMovie = _Value(False)
        protocol.doApplyDoseFilter = _Value(False)
        protocol.splitEvenOdd = _Value(False)
        protocol._lock = threading.RLock()

        class _InputMovie:
            def getObjId(self):
                return 7

        class _FakeMovie:
            def __init__(self):
                self._obj_id = None

            def copyInfo(self, movie):
                self._obj_id = movie.getObjId()

            def getObjId(self):
                return self._obj_id

            def setFileName(self, _):
                pass

            def setMicName(self, _):
                pass

            def getNumberOfFrames(self):
                return 1

            def setAlignment(self, _):
                pass

        class _OutputMovies:
            def __init__(self):
                self.ids = [7]

            def __len__(self):
                return len(self.ids)

            def __contains__(self, obj_id):
                return obj_id in self.ids

            def append(self, item):
                self.ids.append(item.getObjId())

            def update(self, _):
                pass

            def write(self):
                pass

        outputMovies = _OutputMovies()
        protocol._getOutputMovies = lambda: outputMovies
        protocol.getMovieAlignment = lambda *args, **kwargs: object()
        protocol._store = lambda *args, **kwargs: None
        protocol.closeOutputsForStreaming = lambda: None

        registered = []
        protocol._registerMics = (
            lambda movieFName, inMovie, outputName, suffix="":
            registered.append(inMovie.getObjId())
        )

        module = "motioncorr.protocols.protocol_motioncorr_ns"
        with patch(module + ".Movie", _FakeMovie):
            protocol.createOutputStep(
                "/data/movie_007.mrcs",
                _InputMovie(),
            )

        self.assertEqual(
            [7],
            outputMovies.ids,
            "Retrying a partially persisted movie must not append it "
            "again to outputMovies.",
        )
        self.assertEqual(
            [7],
            registered,
            "An already persisted movie must still continue registering "
            "the missing downstream outputs.",
        )

    def test_ReadingOutputIgnoresCompletedIdsOutsideCurrentInput(self):
        protocol = _NewStreamingHarness()
        protocol.itemIdReadList = []
        protocol.splitEvenOdd = _Value(False)
        protocol.doApplyDoseFilter = _Value(False)

        class _Item:
            def __init__(self, obj_id):
                self._obj_id = obj_id

            def getObjId(self):
                return self._obj_id

        class _Set:
            def __init__(self, ids):
                self._items = [_Item(obj_id) for obj_id in ids]

            def __iter__(self):
                return iter(self._items)

            def __len__(self):
                return len(self._items)

            def getUniqueValues(self, field):
                if field != "id":
                    raise AssertionError(f"Unexpected field: {field}")
                return [item.getObjId() for item in self._items]

        input_movies = _Set([1])
        protocol.getInputMovies = lambda *args, **kwargs: input_movies

        setattr(
            protocol,
            protocol._possibleOutputs.movies.name,
            _Set([1, 99]),
        )
        setattr(
            protocol,
            protocol._possibleOutputs.micrographs.name,
            _Set([1, 99]),
        )

        protocol.readingOutput()

        self.assertEqual(
            [1],
            protocol.itemIdReadList,
            "Resume must ignore persisted output IDs that are not part "
            "of the current input set.",
        )

    def test_CreateOutputDoesNotSwallowMicrographStoreFailure(self):
        protocol = _NewStreamingHarness()
        protocol.doSaveMovie = _Value(False)
        protocol.doApplyDoseFilter = _Value(False)
        protocol.splitEvenOdd = _Value(False)
        protocol._lock = __import__("threading").RLock()

        class _InputMovie:
            def getObjId(self):
                return 7

        class _OutputMovies:
            def __len__(self):
                return 1

            def __contains__(self, obj_id):
                return obj_id == 7

        protocol._getOutputMovies = lambda: _OutputMovies()
        protocol.closeOutputsForStreaming = lambda: None

        def fail_micrograph_registration(*args, **kwargs):
            raise RuntimeError("micrograph store failed")

        protocol._registerMics = fail_micrograph_registration

        with self.assertRaisesRegex(RuntimeError, "micrograph store failed"):
            protocol.createOutputStep(
                "/data/movie_007.mrcs",
                _InputMovie(),
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
