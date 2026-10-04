
import threading
from unittest import TestCase
from unittest.mock import patch
import pyworkflow.object as pwobj

import motioncorr.protocols.protocol_motioncorr_ns as motioncorrNs
import motioncorr.protocols.protocol_motioncorr_tasks as motioncorrTasks
from motioncorr.protocols.protocol_motioncorr_ns import (
    MotionCorrOutputs,
    ProtMotionCorrNewStreaming,
)
import motioncorr.protocols.protocol_motioncorr as motioncorrLegacy


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

    def getSize(self):
        return 0


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
        return True

    def _hasValidDose(self):
        return True

    def closeOutputsForStreaming(self):
        return None

    def _updateOutputMoviesOptics(self, outputMovies):
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
        self.failedMovies = []

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

    def test_CreateOutputStepSkipsDoseWeightedOutputWhenDoseUnavailable(self):
        # Regression test: createOutputStep used to request the DW
        # (dose-weighted) output purely based on the doApplyDoseFilter
        # form flag, regardless of whether a usable dose was actually
        # available. _getMcArgs only asks MotionCor2 to dose-weight
        # when _hasValidDose() is True, so MotionCor2 never writes a
        # DW file when it is False - createOutputStep must request the
        # plain (non-DW) output in that case instead of looking for a
        # file that was never produced.
        protocol = _ProtocolStub()
        protocol.doApplyDoseFilter = _ValueStub(True)
        protocol._hasValidDose = lambda: False
        registerCalls = []

        def _recordRegisterMics(movieFName, inMovie, outputName, suffix=""):
            registerCalls.append((outputName, suffix))
            return True

        protocol._registerMics = _recordRegisterMics
        inputMovie = _InputMovieStub(objId=37)

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

        self.assertEqual(1, len(registerCalls))
        outputName, suffix = registerCalls[0]
        self.assertEqual("", suffix)
        self.assertEqual(
            MotionCorrOutputs.micrographs.name,
            outputName,
        )

    def test_CreateOutputStepStopsWhenRegisterMicsFails(self):
        # createOutputStep must not proceed to splitEvenOdd registration
        # or closeOutputsForStreaming for a movie whose primary
        # micrograph registration already failed and was routed to
        # failedMovies.
        protocol = _ProtocolStub()
        protocol.splitEvenOdd = _ValueStub(True)
        protocol._registerMics = lambda *args, **kwargs: False
        closeCalls = []
        protocol.closeOutputsForStreaming = lambda: closeCalls.append(1)
        inputMovie = _InputMovieStub(objId=37)

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

        self.assertEqual([], closeCalls)

    def test_RegisterMicsRoutesMissingOutputFileToFailedMoviesInsteadOfCrashing(self):
        # Regression test: _registerMics never checked that the
        # motion-correction output file actually existed before
        # reading its header (setMRCSamplingRate) - the external tool
        # can return success without producing every expected output
        # for a given movie (e.g. too few frames for dose weighting),
        # and createOutputStep's SET OF MICROGRAPHS block deliberately
        # re-raises on failure for output-Set-integrity reasons, so a
        # single missing micrograph file crashed the whole protocol
        # instead of being routed to failedMovies like every other
        # per-movie failure in this class.
        protocol = _MicrographProtocolStub()
        inputMovie = _InputMovieStub(objId=37)

        def _failSetMRCSamplingRate(*args, **kwargs):
            raise FileNotFoundError(
                "[Errno 2] No such file or directory: '/tmp/mic-37.mrc'"
            )

        with patch.object(
            motioncorrNs,
            "Micrograph",
            _OutputMicrographStub,
        ), patch.object(
            motioncorrNs,
            "setMRCSamplingRate",
            _failSetMRCSamplingRate,
        ):
            result = ProtMotionCorrNewStreaming._registerMics(
                protocol,
                "/tmp/movie-37.mrc",
                inputMovie,
                "micrographs",
                suffix="",
            )

        self.assertFalse(result)
        self.assertEqual(["/tmp/movie-37.mrc"], protocol.failedMovies)
        self.assertIsNone(protocol.outputMicsSet.appended)

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

    def test_OutputChecksUseUnweightedMicrographsWhenDoseIsUnavailable(self):
        protocol = _ProtocolStub()
        protocol.doApplyDoseFilter = _ValueStub(True)
        protocol._hasValidDose = lambda: False

        outputs = ProtMotionCorrNewStreaming._getOutputsToCheck(
            protocol
        )

        self.assertIn(
            MotionCorrOutputs.micrographs.name,
            outputs,
        )
        self.assertNotIn(
            MotionCorrOutputs.micrographsDW.name,
            outputs,
        )

    def test_ReadingOutputUsesUnweightedMicrographsWhenDoseIsUnavailable(self):
        class ItemStub:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

        class SetStub:
            def __init__(self, objIds):
                self._items = [
                    ItemStub(objId)
                    for objId in objIds
                ]

            def __iter__(self):
                return iter(self._items)

            def __len__(self):
                return len(self._items)

        class InputSetStub:
            def getUniqueValues(self, attribute):
                if attribute != "id":
                    raise AssertionError(
                        "Unexpected attribute: %s" % attribute
                    )
                return [1]

        class ProtocolStub:
            _possibleOutputs = MotionCorrOutputs

            def __init__(self):
                self.itemIdReadList = []
                self.splitEvenOdd = _ValueStub(False)
                self.doApplyDoseFilter = _ValueStub(True)
                self.inputSet = InputSetStub()
                self._hasValidDose = lambda: False

                setattr(
                    self,
                    self._possibleOutputs.movies.name,
                    SetStub([1]),
                )
                setattr(
                    self,
                    self._possibleOutputs.micrographs.name,
                    SetStub([1]),
                )

            def getInputMovies(self):
                return self.inputSet

            def info(self, message):
                return None

        protocol = ProtocolStub()

        ProtMotionCorrNewStreaming.readingOutput(
            protocol
        )

        self.assertEqual(
            [1],
            protocol.itemIdReadList,
            "Resume must require the same micrograph output that "
            "createOutputStep actually publishes.",
        )

    def test_InitializeRejectsMissingDoseOnlyAfterFirstMovieExists(self):
        class FirstMovie:
            def getFileName(self):
                return "/tmp/movie-1.mrc"

        class InputMovies:
            def getGain(self):
                return None

            def getDark(self):
                return None

            def getSamplingRate(self):
                return 1.0

            def getFirstItem(self):
                return FirstMovie()

            def isStreamOpen(self):
                return True

        class ProtocolStub:
            def __init__(self):
                self._lock = threading.RLock()
                self.binFactor = _ValueStub(1.0)
                self.doApplyDoseFilter = _ValueStub(True)
                self.inputSet = InputMovies()
                self._hasValidDose = lambda: False

            def getInputMovies(self):
                return self.inputSet

        protocol = ProtocolStub()

        with self.assertRaisesRegex(
            RuntimeError,
            "dose",
        ):
            ProtMotionCorrNewStreaming._initialize(
                protocol
            )

    def test_PreparedGainAndDarkSurviveInputPropertyRefresh(self):
        class MutableInputSet:
            def __init__(self):
                self.originalGain = "/data/gain.dm4"
                self.originalDark = "/data/dark.dm4"
                self.gain = self.originalGain
                self.dark = self.originalDark

            def getGain(self):
                return self.gain

            def setGain(self, value):
                self.gain = value

            def getDark(self):
                return self.dark

            def setDark(self, value):
                self.dark = value

            def getSamplingRate(self):
                return 1.0

            def loadAllProperties(self):
                # Model the streaming refresh: persisted input Set
                # properties are authoritative and restore the original
                # correction-image paths.
                self.gain = self.originalGain
                self.dark = self.originalDark

        class Harness:
            def __init__(self):
                self.inputSet = MutableInputSet()
                self.isEER = False
                self.cropDimX = _ValueStub(0)
                self.cropDimY = _ValueStub(0)
                self.patchX = _ValueStub(0)
                self.patchY = _ValueStub(0)
                self.cropOffsetX = 0
                self.cropOffsetY = 0
                self.binFactor = _ValueStub(1.0)
                self.tol = _ValueStub(0.5)
                self.doSaveMovie = False
                self.doApplyDoseFilter = False
                self.group = _ValueStub(1)
                self.groupLocal = _ValueStub(1)
                self.splitEvenOdd = False
                self.defectFile = _ValueStub("")
                self.defectMap = _ValueStub("")
                self.doMagCor = False
                self.gainRot = _ValueStub(0)
                self.gainFlip = _ValueStub(0)

                setattr(
                    self,
                    "_ProtMotionCorrBase__convertCorrectionImage",
                    lambda image: {
                        "/data/gain.dm4": "/extra/gain.mrc",
                        "/data/dark.dm4": "/extra/dark.mrc",
                    }.get(image, image),
                )

            def getInputMovies(self):
                return self.inputSet

            def _prepareEERFiles(self):
                return None

            def _getFramesRange(self):
                return 1, 10

            def _getNumberOfFrames(self):
                return 10

            def _getCachedAcquisitionValues(self):
                return 300.0, 0.0, 1.0

            def _getExtraPath(self, *parts):
                return "/extra/" + "/".join(parts)

            def getAttributeValue(self, name):
                return None

        protocol = Harness()

        motioncorrNs.ProtMotionCorrBase._convertInputStep(
            protocol
        )

        self.assertEqual(
            "/extra/gain.mrc",
            protocol.inputSet.getGain(),
        )
        self.assertEqual(
            "/extra/dark.mrc",
            protocol.inputSet.getDark(),
        )

        protocol.inputSet.loadAllProperties()

        self.assertEqual(
            "/data/gain.dm4",
            protocol.inputSet.getGain(),
            "The regression must model the input Set refresh restoring "
            "its persisted correction-image metadata.",
        )

        args = motioncorrNs.ProtMotionCorrBase._getMcArgs(
            protocol
        )

        self.assertEqual(
            '"/extra/gain.mrc"',
            args["-Gain"],
            "Movies scheduled after loadAllProperties() must keep using "
            "the gain prepared by the shared input-preparation step.",
        )
        self.assertEqual(
            "/extra/dark.mrc",
            args["-Dark"],
            "Movies scheduled after loadAllProperties() must keep using "
            "the dark prepared by the shared input-preparation step.",
        )

    def test_NewStreamingBuildsMovieOpticsOnlyAfterFirstOutputAppend(self):
        import inspect

        get_output_source = inspect.getsource(
            motioncorrNs.ProtMotionCorrNewStreaming._getOutputMovies
        )
        create_output_source = inspect.getsource(
            motioncorrNs.ProtMotionCorrNewStreaming.createOutputStep
        )

        self.assertNotIn(
            "OpticsGroups.fromImages(outputMovies)",
            get_output_source,
            "NewStreaming must not build Relion optics from a newly created "
            "empty outputMovies Set.",
        )

        optics_call = "self._updateOutputMoviesOptics(outputMovies)"
        append_call = "outputMovies.append(outMovie)"

        self.assertIn(optics_call, create_output_source)
        self.assertIn(append_call, create_output_source)
        self.assertLess(
            create_output_source.index(append_call),
            create_output_source.index(optics_call),
            "Relion optics must be built only after the first output movie "
            "has populated the Set dimensions.",
        )

class TestMotionCorrTasksStreamingRuntime(TestCase):
    def test_ProcessAllMoviesUsesLogicalInputWithoutBackingFile(self):
        class LogicalMovie:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

        class LogicalInputSet:
            def __init__(self):
                self.movies = [
                    LogicalMovie(1),
                    LogicalMovie(2),
                ]
                self.iterCalls = 0

            def getFileName(self):
                raise AssertionError(
                    "ProtMotionCorrTasks must not require a physical "
                    "SQLite backing file for streaming input."
                )

            def iterItems(self):
                self.iterCalls += 1
                return iter(self.movies)

            def getUniqueValues(self, attribute):
                if attribute != "id":
                    raise AssertionError(
                        "Unexpected attribute: %s" % attribute
                    )
                return [
                    movie.getObjId()
                    for movie in self.movies
                ]

            def isStreamOpen(self):
                return False

        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class BatchManagerStub:
            def __init__(
                    self,
                    batchSize,
                    moviesIter,
                    tmpPath,
            ):
                self.items = list(moviesIter)

            def generate(self):
                return iter(())

        class PipelineNode:
            outputQueue = None

        class PipelineStub:
            def addGenerator(self, generator):
                # Force the logical iterator/batch source to be created.
                list(generator())
                return PipelineNode()

            def addProcessor(
                    self,
                    inputQueue,
                    processor,
                    outputQueue=None,
            ):
                return PipelineNode()

            def run(self):
                return None

        class Harness:
            def __init__(self):
                self.inputSet = LogicalInputSet()
                self.streamingSleepOnWait = ValueStub(0)
                self.streamingBatchSize = ValueStub(10)
                self.lock = None
                self._firstTimeOutput = None

            def getInputMovies(self):
                return self.inputSet

            def _createOutputMovies(self):
                return True

            def _createOutputMicrographs(self):
                return False

            def _createOutputWeightedMicrographs(self):
                return False

            def _doSplitEvenOdd(self):
                return False

            def error(self, message):
                return None

            def debug(self, message):
                return None

            def _getCmd(self):
                return "motioncor"

            def _getOutputItemIds(self, outputSet):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getOutputItemIds(outputSet)
                )

            def _getRequiredOutputNamesForResume(self):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getRequiredOutputNamesForResume(self)
                )

            def _getPersistedOutputMovieIds(self):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getPersistedOutputMovieIds(self)
                )

            def _iterLogicalInputMovies(self, waitSecs):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._iterLogicalInputMovies(self, waitSecs)
                )

            def getGpuList(self):
                return []

            def _getTmpPath(self):
                return "/tmp"

            def _updateOutputSets(self, movies, streamState):
                return None

            def _moveBatchOutput(self, batch):
                return batch
            _getOutputProcessor = motioncorrTasks.ProtMotionCorrTasks._getOutputProcessor

        protocol = Harness()

        with patch.object(
            motioncorrTasks,
            "BatchManager",
            BatchManagerStub,
        ), patch.object(
            motioncorrTasks,
            "Pipeline",
            PipelineStub,
        ), patch.object(
            motioncorrTasks.Plugin,
            "getProgram",
            return_value="MotionCor2",
        ):
            motioncorrTasks.ProtMotionCorrTasks._processAllMoviesStep(
                protocol,
                "now",
            )

        self.assertGreater(
            protocol.inputSet.iterCalls,
            0,
            "Tasks streaming must discover input Movies through the "
            "logical Set API.",
        )

    def test_LogicalResumeRequiresAllConfiguredOutputsToBePersisted(self):
        class MovieStub:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

            def clone(self):
                return MovieStub(self._objId)

        class InputSetStub:
            def __init__(self):
                self.movies = [
                    MovieStub(1),
                    MovieStub(2),
                    MovieStub(3),
                ]

            def iterItems(self):
                return iter(self.movies)

            def isStreamOpen(self):
                return False

        class OutputSetStub:
            def __init__(self, objIds):
                self.objIds = list(objIds)

            def getUniqueValues(self, attribute):
                if attribute != "id":
                    raise AssertionError(
                        "Unexpected attribute: %s" % attribute
                    )
                return list(self.objIds)

        class Harness:
            def __init__(self):
                self.inputSet = InputSetStub()
                self.outputMovies = OutputSetStub([1, 2])
                self.outputMicrographs = OutputSetStub([1])

            def getInputMovies(self):
                return self.inputSet

            def _getOutputItemIds(self, outputSet):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getOutputItemIds(outputSet)
                )

            def _getRequiredOutputNamesForResume(self):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getRequiredOutputNamesForResume(self)
                )

            def _getPersistedOutputMovieIds(self):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getPersistedOutputMovieIds(self)
                )

            def _iterLogicalInputMovies(self, waitSecs):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._iterLogicalInputMovies(self, waitSecs)
                )

            def _createOutputMovies(self):
                return True

            def _createOutputMicrographs(self):
                return True

            def _createOutputWeightedMicrographs(self):
                return False

            def _doSplitEvenOdd(self):
                return False

        protocol = Harness()

        pendingIds = [
            movie.getObjId()
            for movie in protocol._iterLogicalInputMovies(0)
        ]

        self.assertEqual(
            [2, 3],
            pendingIds,
            "A Movie is complete only when every output required by the "
            "current Tasks configuration has been persisted. outputMovies "
            "alone must not advance the resume checkpoint.",
        )

    def test_OutputPublicationDoesNotUseSqliteOutputLoader(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

            def __bool__(self):
                return bool(self.value)

        class MovieStub:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

        class Harness:
            def __init__(self):
                self.doSaveMovie = ValueStub(False)
                self._firstTimeOutput = False

            def _createOutputMovies(self):
                return True

            def _createOutputMicrographs(self):
                return False

            def _createOutputWeightedMicrographs(self):
                return False

            def _doSplitEvenOdd(self):
                return False

            def getAttributeValue(self, name, default=None):
                value = getattr(self, name, default)
                getter = getattr(value, "get", None)
                return getter() if callable(getter) else value

            def _updateLogicalOutputMovieSet(
                    self,
                    newDone,
                    streamMode,
            ):
                self.logicalMovieUpdates = (
                    getattr(
                        self,
                        "logicalMovieUpdates",
                        0,
                    )
                    + 1
                )

            def _loadOutputSet(self, *args, **kwargs):
                raise AssertionError(
                    "ProtMotionCorrTasks output publication must not "
                    "load or create SQLite-backed output Sets."
                )

        protocol = Harness()

        motioncorrTasks.ProtMotionCorrTasks._updateOutputSets(
            protocol,
            [MovieStub(1)],
            pwobj.Set.STREAM_OPEN,
        )

        self.assertEqual(
            1,
            protocol.logicalMovieUpdates,
            "Tasks must route output publication through the logical "
            "output path.",
        )

    def test_DoseWeightedOutputRequiresUsableDose(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

            def __bool__(self):
                return bool(self.value)

        class Harness:
            def __init__(self):
                self.doApplyDoseFilter = ValueStub(True)
                self.doSaveUnweightedMic = ValueStub(False)

            def _hasValidDose(self):
                return False

            def _createOutputWeightedMicrographs(self):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._createOutputWeightedMicrographs(self)
                )

            def _createOutputMicrographs(self):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._createOutputMicrographs(self)
                )

        protocol = Harness()

        self.assertFalse(
            protocol._createOutputWeightedMicrographs(),
            "Tasks must not require a dose-weighted output when dose "
            "metadata is unavailable, even if the user enabled dose "
            "filtering.",
        )
        self.assertTrue(
            protocol._createOutputMicrographs(),
            "When dose weighting cannot actually be applied, Tasks must "
            "fall back to the regular micrograph output.",
        )

    def test_ExistingLogicalOutputIsReusedWithoutDuplicateAppend(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class MovieStub:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

            def getFileName(self):
                return "/tmp/movie-%d.mrc" % self._objId

        class OutputMoviesStub:
            def __init__(self):
                self.ids = {1}
                self.enableAppendCalls = 0
                self.appendCalls = 0
                self.updateCalls = 0
                self.writeCalls = 0
                self.streamStates = []

            def enableAppend(self):
                self.enableAppendCalls += 1

            def __contains__(self, objId):
                return objId in self.ids

            def append(self, movie):
                self.appendCalls += 1
                self.ids.add(movie.getObjId())

            def update(self, movie):
                self.updateCalls += 1

            def setStreamState(self, state):
                self.streamStates.append(state)

            def write(self):
                self.writeCalls += 1

        class Harness:
            def __init__(self):
                self.outputMovies = OutputMoviesStub()
                self.doSaveMovie = ValueStub(False)
                self.storeCalls = []

            def getAttributeValue(self, name, default=None):
                value = getattr(self, name, default)
                getter = getattr(value, "get", None)
                return getter() if callable(getter) else value

            def _getLogicalOutputSet(
                    self,
                    outputName,
                    setClass,
                    suffix="",
                    fixSampling=True,
            ):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._getLogicalOutputSet(
                        self,
                        outputName,
                        setClass,
                        suffix=suffix,
                        fixSampling=fixSampling,
                    )
                )

            def _persistLogicalOutputSet(self, outputSet, streamMode):
                return (
                    motioncorrTasks.ProtMotionCorrTasks
                    ._persistLogicalOutputSet(
                        self,
                        outputSet,
                        streamMode,
                    )
                )

            def _createOutputMovie(self, movie):
                raise AssertionError(
                    "An already persisted Movie must not be rebuilt."
                )

            def _defineOutputs(self, **kwargs):
                raise AssertionError(
                    "An existing logical output Set must be reused, "
                    "not redefined."
                )

            def _store(self, outputSet):
                self.storeCalls.append(outputSet)
            @staticmethod
            def _getOutputItemIds(outputSet):
                # These legacy unit-test stubs model persisted
                # membership through __contains__ instead of the
                # logical Set iteration/getUniqueValues API. Adapt
                # that old fixture contract without changing
                # production behavior.
                class PersistedIds:
                    def __init__(self, persistedOutputSet):
                        self.outputSet = persistedOutputSet
                        self.added = set()

                    def __contains__(self, objId):
                        return (objId in self.added or
                                objId in self.outputSet)

                    def add(self, objId):
                        self.added.add(objId)

                return PersistedIds(outputSet)

        protocol = Harness()
        originalSet = protocol.outputMovies

        motioncorrTasks.ProtMotionCorrTasks._updateLogicalOutputMovieSet(
            protocol,
            [MovieStub(1)],
            pwobj.Set.STREAM_OPEN,
        )

        self.assertIs(
            originalSet,
            protocol.outputMovies,
            "Tasks must reuse the existing logical output Set.",
        )
        self.assertEqual(
            1,
            originalSet.enableAppendCalls,
            "An existing output Set must be reopened for append.",
        )
        self.assertEqual(
            0,
            originalSet.appendCalls,
            "Resume/replay must not append an already persisted objId.",
        )
        self.assertEqual(
            0,
            originalSet.updateCalls,
            "Resume/replay must not update an already persisted objId.",
        )
        self.assertEqual(
            [pwobj.Set.STREAM_OPEN],
            originalSet.streamStates,
        )
        self.assertEqual(
            1,
            originalSet.writeCalls,
            "The reused Set must still persist its current stream state.",
        )
        self.assertEqual(
            [originalSet],
            protocol.storeCalls,
        )

    def test_NewLogicalOutputUsesCanonicalSetPublishedByDefineOutputs(self):
        class InputSetStub:
            def getSamplingRate(self):
                return 1.5

        class OutputSetStub:
            def __init__(self, name):
                self.name = name
                self.copyInfoCalls = 0
                self.samplingRates = []
                self.streamStates = []
                self.writeCalls = 0

            def copyInfo(self, inputSet):
                self.copyInfoCalls += 1

            def setSamplingRate(self, value):
                self.samplingRates.append(value)

            def setStreamState(self, state):
                self.streamStates.append(state)

            def write(self):
                self.writeCalls += 1

        provisional = OutputSetStub("provisional")
        canonical = OutputSetStub("canonical")

        class SetClassStub:
            @classmethod
            def create(
                    cls,
                    outputPath,
                    template=None,
                    suffix=None,
            ):
                return provisional

        class Harness:
            def __init__(self):
                self.inputSet = InputSetStub()
                self.inputMovies = object()
                self.definedOutputs = []
                self.relations = []

            def getInputMovies(self):
                return self.inputSet

            def _getPath(self):
                return "/tmp/protocol"

            def _getBinFactor(self):
                return 1.0

            def _defineOutputs(self, **kwargs):
                self.definedOutputs.append(kwargs)
                # Model the runtime adapter replacing the provisional Set
                # with its canonical persisted representation.
                self.outputMovies = canonical

            def _defineTransformRelation(self, source, target):
                self.relations.append((source, target))

        protocol = Harness()

        outputSet, created = (
            motioncorrTasks.ProtMotionCorrTasks
            ._getLogicalOutputSet(
                protocol,
                "outputMovies",
                SetClassStub,
                fixSampling=True,
            )
        )

        self.assertTrue(created)
        self.assertIs(
            canonical,
            outputSet,
            "After _defineOutputs(), Tasks must continue with the "
            "canonical output Set published by the runtime rather than "
            "the provisional Set used to create it.",
        )
        self.assertIs(
            canonical,
            protocol.outputMovies,
        )

    def test_MovieAppendPersistenceFailureIsNotSwallowed(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class AlignmentStub:
            def getShifts(self):
                return ([0.0, 1.0], [0.0, 1.0])

        class MovieStub:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

            def getFileName(self):
                return "/tmp/movie-%d.mrc" % self._objId

            def copyObjId(self, other):
                self._objId = other.getObjId()

            def getAlignment(self):
                return AlignmentStub()

        class OutputMoviesStub:
            def enableAppend(self):
                pass

            def __contains__(self, objId):
                return False

            def append(self, movie):
                raise RuntimeError(
                    "simulated outputMovies append persistence failure"
                )

            def update(self, movie):
                raise AssertionError(
                    "update() must not be reached after append() fails"
                )

        class Harness:
            def __init__(self):
                self.outputMovies = OutputMoviesStub()
                self.doSaveMovie = ValueStub(False)
                self.loggedErrors = []

            def getAttributeValue(self, name, default=None):
                value = getattr(self, name, default)
                getter = getattr(value, "get", None)
                return getter() if callable(getter) else value

            def _getLogicalOutputSet(
                    self,
                    outputName,
                    setClass,
                    suffix="",
                    fixSampling=True,
            ):
                return self.outputMovies, False

            def _createOutputMovie(self, movie):
                return MovieStub(movie.getObjId())

            def _persistLogicalOutputSet(self, outputSet, streamMode):
                raise AssertionError(
                    "Set persistence must not continue after append() fails."
                )

            def error(self, message):
                self.loggedErrors.append(message)
            @staticmethod
            def _getOutputItemIds(outputSet):
                # These legacy unit-test stubs model persisted
                # membership through __contains__ instead of the
                # logical Set iteration/getUniqueValues API. Adapt
                # that old fixture contract without changing
                # production behavior.
                class PersistedIds:
                    def __init__(self, persistedOutputSet):
                        self.outputSet = persistedOutputSet
                        self.added = set()

                    def __contains__(self, objId):
                        return (objId in self.added or
                                objId in self.outputSet)

                    def add(self, objId):
                        self.added.add(objId)

                return PersistedIds(outputSet)

        protocol = Harness()

        with self.assertRaisesRegex(
            RuntimeError,
            "simulated outputMovies append persistence failure",
        ):
            motioncorrTasks.ProtMotionCorrTasks._updateLogicalOutputMovieSet(
                protocol,
                [MovieStub(7)],
                pwobj.Set.STREAM_OPEN,
            )

    def test_GetCmdWaitsForFirstMovieWhenLogicalStreamStartsEmpty(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class MovieStub:
            def getFileName(self):
                return "/tmp/movie-1.mrc"

        class InputSetStub:
            def __init__(self):
                self.refreshed = False
                self.refreshCalls = 0

            def getFirstItem(self):
                return MovieStub() if self.refreshed else None

            def isStreamOpen(self):
                return True

            def loadAllProperties(self):
                self.refreshCalls += 1
                self.refreshed = True

        class Harness:
            def __init__(self):
                self.inputSet = InputSetStub()
                self.extraParams2 = ValueStub("")
                self.streamingSleepOnWait = ValueStub(0)

            def getInputMovies(self):
                return self.inputSet

            def _getMcArgs(self):
                return {}

        protocol = Harness()

        with patch.object(
            motioncorrTasks.time,
            "sleep",
            return_value=None,
        ):
            command = motioncorrTasks.ProtMotionCorrTasks._getCmd(
                protocol,
            )

        self.assertIn("-InMrc ./", command)
        self.assertEqual(
            1,
            protocol.inputSet.refreshCalls,
            "Tasks must refresh an open logical input instead of "
            "dereferencing None when it starts empty.",
        )

    def test_BatchFailureAfterLastRetryPropagatesWithoutFinalSleep(self):
        class LockStub:
            def __enter__(self):
                return self

            def __exit__(self, excType, excValue, traceback):
                return False

        class Harness:
            def __init__(self):
                self.command = "motioncor -Gpu #"
                self.program = "MotionCor3"
                self.lock = LockStub()
                self.processed = 0
                self.runJobCalls = 0
                self._batchFailures = []
                self._failedBatchIds = set()

            def runJob(self, program, command, cwd=None):
                self.runJobCalls += 1
                raise RuntimeError("simulated MotionCor batch failure")

            def debug(self, message):
                return None

            def error(self, message):
                return None

            def batch_str(self, batch):
                return "batch"

        protocol = Harness()
        batch = {
            "id": "batch-1",
            "index": 0,
            "items": [object()],
            "path": "/tmp/motioncorr-tasks-retry-test",
        }

        processor = (
            motioncorrTasks.ProtMotionCorrTasks
            ._getMcProcessor(protocol, 0)
        )

        with patch.object(
            motioncorrTasks.os.path,
            "exists",
            return_value=True,
        ), patch.object(
            motioncorrTasks.time,
            "sleep",
            return_value=None,
        ) as sleepMock:
            returnedBatch = processor(batch)

        self.assertIs(
            batch,
            returnedBatch,
            "A worker failure must stay inside the emtools Pipeline thread "
            "so its output queue can be closed normally.",
        )
        self.assertEqual(
            2,
            protocol.runJobCalls,
            "Tasks must try a failed MotionCor batch exactly twice.",
        )
        self.assertEqual(
            1,
            sleepMock.call_count,
            "Tasks must sleep only between retries, never after the "
            "final failed attempt.",
        )
        self.assertEqual(1, len(protocol._batchFailures))
        self.assertRegex(
            str(protocol._batchFailures[0]),
            "simulated MotionCor batch failure",
        )

    def test_ProcessAllMoviesRaisesWorkerFailureAfterPipelineDrains(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class BatchManagerStub:
            def __init__(self, batchSize, moviesIter, tmpPath):
                pass

            def generate(self):
                return iter(())

        class PipelineNode:
            outputQueue = None

        protocol = None

        class PipelineStub:
            def addGenerator(self, generator):
                return PipelineNode()

            def addProcessor(
                    self,
                    inputQueue,
                    processor,
                    outputQueue=None,
            ):
                return PipelineNode()

            def run(self):
                protocol._batchFailures = [
                    RuntimeError(
                        "simulated worker-thread MotionCor failure"
                    )
                ]

        class Harness:
            def __init__(self):
                self.streamingSleepOnWait = ValueStub(0)
                self.streamingBatchSize = ValueStub(10)
                self.closedOutputs = 0

            def error(self, message):
                return None

            def debug(self, message):
                return None

            def _getCmd(self):
                return "motioncor"

            def _iterLogicalInputMovies(self, waitSecs):
                return iter(())

            def _getTmpPath(self):
                return "/tmp"

            def getGpuList(self):
                return [0]

            def _getMcProcessor(self, gpu):
                return lambda batch: batch

            def _moveBatchOutput(self, batch):
                return batch

            def _updateOutputSets(self, movies, streamState):
                self.closedOutputs += 1
            _getOutputProcessor = motioncorrTasks.ProtMotionCorrTasks._getOutputProcessor

        protocol = Harness()

        with patch.object(
            motioncorrTasks,
            "BatchManager",
            BatchManagerStub,
        ), patch.object(
            motioncorrTasks,
            "Pipeline",
            PipelineStub,
        ), patch.object(
            motioncorrTasks.Plugin,
            "getProgram",
            return_value="MotionCor3",
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "simulated worker-thread MotionCor failure",
            ):
                motioncorrTasks.ProtMotionCorrTasks._processAllMoviesStep(
                    protocol,
                    "now",
                )

        self.assertEqual(
            0,
            protocol.closedOutputs,
            "Tasks must not close outputs as successful after a worker "
            "failure; the main Scipion step must fail instead.",
        )

    def test_LegacyMotionCorrKeepsDoneCompatibilityWithoutTasksSidecar(self):
        class Harness:
            def _getExtraPath(self, *parts):
                return "/tmp/" + "/".join(parts)

        legacy = Harness()
        tasks = Harness()

        with patch.object(
            motioncorrLegacy.ProtMotionCorrBase,
            "_convertInputStep",
            autospec=True,
            return_value=None,
        ) as baseConvert, patch.object(
            motioncorrLegacy.pwutils,
            "makePath",
        ) as makePath:
            motioncorrLegacy.ProtMotionCorr._convertInputStep(
                legacy,
            )

            baseConvert.assert_called_once_with(legacy)
            makePath.assert_called_once_with("/tmp/DONE")

        with patch.object(
            motioncorrLegacy.ProtMotionCorrBase,
            "_convertInputStep",
            autospec=True,
            return_value=None,
        ) as baseConvert, patch.object(
            motioncorrLegacy.pwutils,
            "makePath",
        ) as makePath:
            motioncorrTasks.ProtMotionCorrTasks._convertInputStep(
                tasks,
            )

            baseConvert.assert_called_once_with(tasks)
            makePath.assert_not_called()

    def test_BatchCommandFiltersSerialInputsByMovieExtension(self):
        class ValueStub:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class MovieStub:
            def __init__(self, fileName):
                self.fileName = fileName

            def getFileName(self):
                return self.fileName

        class InputSetStub:
            def __init__(self, fileName):
                self.movie = MovieStub(fileName)

            def getFirstItem(self):
                return self.movie

            def isStreamOpen(self):
                return False

        class Harness:
            def __init__(self, fileName):
                self.inputSet = InputSetStub(fileName)
                self.extraParams2 = ValueStub("")

            def getInputMovies(self):
                return self.inputSet

            def _getMcArgs(self):
                # Deliberately empty: the suffix must be added by
                # ProtMotionCorrTasks._getCmd itself.
                return {}

        protocol = Harness("/tmp/movie.tif")

        cmd = motioncorrTasks.ProtMotionCorrTasks._getCmd(
            protocol,
        )

        self.assertIn("-InTiff ./", cmd)
        self.assertIn(
            "-InSuffix .tif",
            cmd,
            "ProtMotionCorrTasks._getCmd must add an input suffix "
            "filter for MotionCor serial directory processing.",
        )

    def test_OutputMovieUpdateDoesNotUseSetContainsMapperLookup(self):
        import inspect

        source = inspect.getsource(
            motioncorrTasks.ProtMotionCorrTasks._updateLogicalOutputMovieSet
        )

        self.assertNotIn(
            "objId in outputMovies",
            source,
            "Logical output idempotency must not call Set.__contains__, "
            "because that forces the physical mapper backend.",
        )
        self.assertIn(
            "_getOutputItemIds(outputMovies)",
            source,
            "Logical output idempotency must derive persisted ids through "
            "the mapper-agnostic helper.",
        )

    def test_OutputMicrographUpdateDoesNotUseSetContainsMapperLookup(self):
        import inspect

        source = inspect.getsource(
            motioncorrTasks.ProtMotionCorrTasks._updateLogicalOutputMicSet
        )

        self.assertNotIn(
            "objId in outputMics",
            source,
            "Logical micrograph idempotency must not call Set.__contains__, "
            "because that forces the physical mapper backend.",
        )
        self.assertIn(
            "_getOutputItemIds(outputMics)",
            source,
            "Logical micrograph idempotency must derive persisted ids through "
            "the mapper-agnostic helper.",
        )

    def test_OutputProcessorRecordsPublicationFailureForMainStep(self):
        import threading

        error = RuntimeError("output publication failed")
        batch = {"id": 7, "items": []}

        class Harness:
            def __init__(self):
                self.lock = threading.Lock()
                self._batchFailures = []

            def _moveBatchOutput(self, currentBatch):
                raise error

        protocol = Harness()

        processor = motioncorrTasks.ProtMotionCorrTasks._getOutputProcessor(
            protocol
        )
        returnedBatch = processor(batch)

        self.assertIs(returnedBatch, batch)
        self.assertEqual(1, len(protocol._batchFailures))
        self.assertIs(error, protocol._batchFailures[0])


class _RefreshRequiredTasksOutput:
    def __init__(self, ids=None):
        self.ids = set(ids or [])
        self.loaded = False
        self.appendEnabled = False

    def loadAllProperties(self):
        self.loaded = True

    def enableAppend(self):
        if not self.loaded:
            raise AssertionError(
                "Existing logical output must be refreshed before enableAppend()."
            )
        self.appendEnabled = True

    def getUniqueValues(self, attribute):
        if not self.loaded:
            raise AssertionError(
                "Persisted logical output must be refreshed before reading ids."
            )
        if attribute != "id":
            raise AssertionError("Unexpected attribute: %s" % attribute)
        return sorted(self.ids)


class TestMotionCorrTasksLogicalOutputRestore(TestCase):
    def testResumeRefreshesPersistedOutputsBeforeReadingCompletedIds(self):
        class Harness:
            def __init__(self):
                self.outputMovies = _RefreshRequiredTasksOutput({1, 2})

            def _getRequiredOutputNamesForResume(self):
                return ["outputMovies"]

            _getOutputItemIds = staticmethod(
                motioncorrTasks.ProtMotionCorrTasks._getOutputItemIds
            )

        protocol = Harness()

        completed = (
            motioncorrTasks.ProtMotionCorrTasks
            ._getPersistedOutputMovieIds(protocol)
        )

        self.assertTrue(protocol.outputMovies.loaded)
        self.assertEqual({1, 2}, completed)

    def testExistingLogicalOutputIsRefreshedBeforeReuse(self):
        class Harness:
            def __init__(self):
                self.outputMovies = _RefreshRequiredTasksOutput({1})

        protocol = Harness()

        outputSet, created = (
            motioncorrTasks.ProtMotionCorrTasks
            ._getLogicalOutputSet(
                protocol,
                "outputMovies",
                object,
            )
        )

        self.assertFalse(created)
        self.assertIs(protocol.outputMovies, outputSet)
        self.assertTrue(outputSet.loaded)
        self.assertTrue(outputSet.appendEnabled)
