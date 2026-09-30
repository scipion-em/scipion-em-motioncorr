
import threading
from unittest import TestCase
from unittest.mock import patch

import motioncorr.protocols.protocol_motioncorr_ns as motioncorrNs
from motioncorr.protocols.protocol_motioncorr_ns import (
    MotionCorrOutputs,
    ProtMotionCorrNewStreaming,
)


class _ValueStub:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _InputMovieStub:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId


class _OutputMovieStub:
    def __init__(self):
        self._objId = None
        self.fileName = None
        self.micName = None
        self.alignment = None

    def copyInfo(self, other):
        # Mirror Scipion copyInfo(): metadata is copied
        # without copying the object identity.
        return None

    def copyObjId(self, other):
        self._objId = other.getObjId()

    def getObjId(self):
        return self._objId

    def setFileName(self, fileName):
        self.fileName = fileName

    def setMicName(self, micName):
        self.micName = micName

    def getNumberOfFrames(self):
        return 10

    def setAlignment(self, alignment):
        self.alignment = alignment


class _OutputSetStub:
    def __init__(self):
        self.appended = None
        self.updated = None
        self.writeCalls = 0

    def __len__(self):
        return (
            0
            if self.appended is None
            else 1
        )

    def __contains__(self, objId):
        return (
            self.appended is not None
            and self.appended.getObjId() == objId
        )

    def append(self, item):
        self.appended = item

    def update(self, item):
        self.updated = item

    def write(self):
        self.writeCalls += 1


class _ProtocolStub:
    _possibleOutputs = MotionCorrOutputs

    def __init__(self):
        self.failedMovies = []
        self.doSaveMovie = _ValueStub(False)
        self.doApplyDoseFilter = _ValueStub(False)
        self.splitEvenOdd = _ValueStub(False)
        self._lock = threading.Lock()
        self.outputMoviesSet = _OutputSetStub()

    def _getOutputMovies(self):
        return self.outputMoviesSet

    def getMovieAlignment(self, movieFName, nFrames):
        return ("alignment", movieFName, nFrames)

    def _store(self, outputSet):
        return None

    def _registerMics(
            self,
            movieFName,
            inMovie,
            outputName,
            suffix="",
    ):
        return None

    def closeOutputsForStreaming(self):
        return None


class _OutputMicrographStub:
    def __init__(self):
        self._objId = None
        self.fileName = None
        self.samplingRate = None

    def copyInfo(self, other):
        # Mirror Scipion copyInfo(): metadata is copied
        # without copying the object identity.
        return None

    def copyObjId(self, other):
        self._objId = other.getObjId()

    def getObjId(self):
        return self._objId

    def setFileName(self, fileName):
        self.fileName = fileName

    def setSamplingRate(self, samplingRate):
        self.samplingRate = samplingRate


class _MicrographProtocolStub:
    def __init__(self):
        self.sRate = 1.5
        self.splitEvenOdd = _ValueStub(False)
        self.outputMicsSet = _OutputSetStub()

    def _getOutputMics(
            self,
            outputName,
            suffix="",
    ):
        return self.outputMicsSet

    def _getResultMicFn(
            self,
            movieFName,
            suffix="",
    ):
        return "/tmp/mic-37.mrc"

    def setMicPlotInfo(
            self,
            mic,
            movieFName,
    ):
        return None

    def setMicsEvenOdd(
            self,
            movieFName,
            mic,
    ):
        return None

    def _store(self, outputSet):
        return None

class _GeneratorMovieStub:
    def __init__(self, objId, fileName):
        self._objId = objId
        self._fileName = fileName

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return self._fileName

    def clone(self):
        return _GeneratorMovieStub(
            self._objId,
            self._fileName,
        )


class _GeneratorInputSetStub:
    def __init__(self, movie):
        self._movie = movie

    def getUniqueValues(self, attribute):
        if attribute != "id":
            raise AssertionError(
                "Unexpected attribute: %s" % attribute
            )
        return [self._movie.getObjId()]

    def isStreamOpen(self):
        return False

    def iterItems(self):
        return iter([self._movie])


class _GeneratorProtocolStub:
    def __init__(self):
        self._lock = threading.Lock()
        self.itemIdReadList = []
        self.insertedSteps = []
        self.inputSet = _GeneratorInputSetStub(
            _GeneratorMovieStub(
                37,
                "/tmp/movie-37.mrc",
            )
        )

    def _initialize(self):
        return None

    def getInputMovies(self):
        return self.inputSet

    def readingOutput(self):
        return None

    def _getOutputsToCheck(self):
        return [
            "outputMovies",
            "outputMicrographs",
        ]

    def _insertFunctionStep(
            self,
            function,
            *args,
            **kwargs
    ):
        stepId = len(self.insertedSteps) + 1

        self.insertedSteps.append(
            {
                "id": stepId,
                "function": function.__name__,
                "args": args,
                "prerequisites": kwargs.get(
                    "prerequisites"
                ),
            }
        )

        return stepId

    def _convertInputStep(self):
        return None

    def convertInputStep(self, movieFName):
        return None

    def processMovieStep(self, movieFName):
        return None

    def createOutputStep(
            self,
            movieFName,
            inMovie
    ):
        return None

    def closeOutputSetStep(self, outputs):
        return None

class _ExistingEmptyOutputSetStub:
    def __init__(self):
        self.appendEnabled = False

    def __len__(self):
        return 0

    def enableAppend(self):
        self.appendEnabled = True

class TestMotionCorrNewStreamingRuntime(TestCase):
    def test_OutputMoviesPreserveInputMovieIdentity(self):
        protocol = _ProtocolStub()
        inputMovie = _InputMovieStub(
            objId=37
        )

        with patch.object(
            motioncorrNs,
            "Movie",
            _OutputMovieStub,
        ):
            ProtMotionCorrNewStreaming.createOutputStep(
                protocol,
                "/tmp/movie-37.mrc",
                inputMovie,
            )

        outputMovie = (
            protocol
            .outputMoviesSet
            .appended
        )

        self.assertIsNotNone(
            outputMovie
        )

        self.assertEqual(
            outputMovie.getObjId(),
            inputMovie.getObjId(),
        )

    def test_OutputMicrographsPreserveInputMovieIdentity(self):
        protocol = _MicrographProtocolStub()
        inputMovie = _InputMovieStub(
            objId=37
        )

        with patch.object(
            motioncorrNs,
            "Micrograph",
            _OutputMicrographStub,
        ), patch.object(
            motioncorrNs,
            "setMRCSamplingRate",
            lambda *args, **kwargs: None,
        ):
            ProtMotionCorrNewStreaming._registerMics(
                protocol,
                "/tmp/movie-37.mrc",
                inputMovie,
                "micrographs",
                suffix="",
            )

        outputMic = (
            protocol
            .outputMicsSet
            .appended
        )

        self.assertIsNotNone(
            outputMic
        )

        self.assertEqual(
            outputMic.getObjId(),
            inputMovie.getObjId(),
        )

    def test_CreateOutputStepIsNotWrappedByStorageRetryDecorator(self):
        self.assertFalse(
            hasattr(
                ProtMotionCorrNewStreaming.createOutputStep,
                "__wrapped__",
            ),
            "Output persistence must not be wrapped by a "
            "storage-backend-specific retry decorator.",
        )

    def test_StreamGeneratorUsesOneSharedInputPreparationStep(self):
        protocol = _GeneratorProtocolStub()

        with patch.object(
            motioncorrNs.time,
            "sleep",
            lambda seconds: None,
        ):
            ProtMotionCorrNewStreaming.stepsGeneratorStep(
                protocol
            )

        stepNames = [
            step["function"]
            for step in protocol.insertedSteps
        ]

        self.assertEqual(
            stepNames,
            [
                "_convertInputStep",
                "processMovieStep",
                "createOutputStep",
                "closeOutputSetStep",
            ],
        )

        self.assertEqual(
            protocol.insertedSteps[1]["prerequisites"],
            protocol.insertedSteps[0]["id"],
        )

    def test_InputPreparationDoesNotCreateFilesystemCheckpoint(self):
        constants = (
            motioncorrNs
            .ProtMotionCorrBase
            ._convertInputStep
            .__code__
            .co_consts
        )

        self.assertNotIn(
            "DONE",
            constants,
            "Input preparation must not use a filesystem "
            "checkpoint to track streaming state.",
        )

    def test_OutputFactoriesReuseAlreadyDefinedEmptySets(self):
        existingMovies = _ExistingEmptyOutputSetStub()
        existingMics = _ExistingEmptyOutputSetStub()

        class ProtocolStub:
            _possibleOutputs = motioncorrNs.MotionCorrOutputs

            def __init__(self):
                setattr(
                    self,
                    self._possibleOutputs.movies.name,
                    existingMovies,
                )
                setattr(
                    self,
                    self._possibleOutputs.micrographs.name,
                    existingMics,
                )

            def getInputMovies(self, asPointer=False):
                return object()

            def _getPath(self):
                return "/tmp"

        protocol = ProtocolStub()

        with patch.object(
            motioncorrNs.SetOfMovies,
            "create",
            side_effect=AssertionError(
                "An already-defined empty movie output "
                "must be reused, not recreated."
            ),
        ):
            movies = (
                ProtMotionCorrNewStreaming
                ._getOutputMovies(
                    protocol
                )
            )

        with patch.object(
            motioncorrNs.SetOfMicrographs,
            "create",
            side_effect=AssertionError(
                "An already-defined empty micrograph output "
                "must be reused, not recreated."
            ),
        ):
            mics = (
                ProtMotionCorrNewStreaming
                ._getOutputMics(
                    protocol,
                    protocol._possibleOutputs.micrographs.name,
                )
            )

        self.assertIs(
            movies,
            existingMovies,
        )
        self.assertTrue(
            existingMovies.appendEnabled
        )

        self.assertIs(
            mics,
            existingMics,
        )
        self.assertTrue(
            existingMics.appendEnabled
        )
