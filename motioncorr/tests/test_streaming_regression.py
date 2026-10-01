# **************************************************************************
# *
# * Regression tests for MotionCorr streaming/resume behaviour.
# *
# **************************************************************************

import threading
import unittest
from unittest.mock import patch

from motioncorr.protocols.protocol_motioncorr import ProtMotionCorr
from motioncorr.protocols.protocol_base import ProtMotionCorrBase


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value

    def __bool__(self):
        # Real pyworkflow Scalar-backed params (e.g. Boolean) are
        # truthy/falsy based on their held value - production code
        # relies on that (e.g. "if self.doApplyDoseFilter:" without
        # .get()), so this fake must match it.
        return bool(self._value)


class _InputMovies:
    def getSamplingRate(self):
        return 1.0


class _Movie:
    def getFileName(self):
        return "movie_001.mrc"


class _MotionCorrHarness:
    def __init__(self):
        self.errors = []
        self.isEER = False
        self.defectFile = _Value(None)
        self.defectMap = _Value(None)
        self.extraParams2 = _Value("")

    def getInputMovies(self):
        return _InputMovies()

    def _getOutputMovieFolder(self, movie):
        return "/tmp"

    def _getOutputMicName(self, movie):
        return "out.mrc"

    def _getMcArgs(self):
        raise ValueError("simulated bad acquisition metadata")

    def _getExtraPath(self, *parts):
        return "/tmp/extra"

    def _getInputFormat(self, inputFn, absPath=False):
        return ProtMotionCorr._getInputFormat(inputFn, absPath=absPath)

    def runJob(self, *args, **kwargs):
        raise AssertionError(
            "runJob must not be called when argument building already "
            "failed."
        )

    def _useWorkerThread(self):
        return False

    def error(self, message):
        self.errors.append(message)


class TestMotionCorrStreamingRegression(unittest.TestCase):

    def testArgBuildingFailureIsToleratedInsteadOfCrashingProtocol(self):
        # Regression test: _getMcArgs()/_getInputFormat() and the rest
        # of the argument-building code used to run BEFORE the try/
        # except. processMovieStep (the pwem base class step calling
        # _processMovie) has no exception boundary of its own, so a
        # single movie failing here (e.g. bad acquisition metadata, or
        # an unsupported file extension raising ValueError from
        # _getInputFormat) would crash the whole protocol instead of
        # being logged and skipped like the runJob failure already was.
        protocol = _MotionCorrHarness()

        ProtMotionCorr._processMovie(protocol, _Movie())

        self.assertEqual(1, len(protocol.errors))
        self.assertIn(
            "simulated bad acquisition metadata",
            protocol.errors[0],
        )


class _Acquisition:
    def __init__(self, doseInitial=None, dosePerFrame=None, voltage=300.0):
        self._doseInitial = doseInitial
        self._dosePerFrame = dosePerFrame
        self._voltage = voltage

    def getDoseInitial(self):
        return self._doseInitial

    def getDosePerFrame(self):
        return self._dosePerFrame

    def getVoltage(self):
        return self._voltage


class _DoseInputMovies:
    def __init__(self, acquisition):
        self._acquisition = acquisition

    def getAcquisition(self):
        return self._acquisition

    def getFramesRange(self):
        return [1, 10, 0]


class _DoseHarness(ProtMotionCorrBase):
    def __init__(self, acquisition):
        self._inputMovies = _DoseInputMovies(acquisition)
        # _getCachedAcquisitionValues() serializes its first population with
        # self._lock (see its docstring) - this harness bypasses the
        # real Protocol.__init__ (which sets up a threading.RLock), so
        # it must provide one itself.
        self._lock = threading.RLock()

    def getInputMovies(self):
        return self._inputMovies

    def _getNumberOfFrames(self):
        return 10


class TestMotionCorrCorrectedDoseRegression(unittest.TestCase):
    # Regression tests: acquisition metadata may not always carry dose
    # values (e.g. an import that didn't set them, or a movie whose
    # acquisition dose is not visible yet at processing time). The
    # tomo branch of _getCorrectedDose already guarded against a
    # missing dose ("dose / ... if dose else 0.0"), but the non-tomo
    # (else) branch did an unguarded "dose * (firstFrame - 1)",
    # crashing with TypeError - both in processMovieStep (via
    # _getMcArgs) and in createOutputStep (via calcFrameMotion) for
    # every single movie once the acquisition's dose is unset.

    def testGetCorrectedDoseDefaultsToZeroWhenDoseValuesMissing(self):
        harness = _DoseHarness(
            _Acquisition(doseInitial=None, dosePerFrame=None))

        preExp, dose = ProtMotionCorrBase._getCorrectedDose(harness)

        self.assertEqual(0.0, preExp)
        self.assertEqual(0.0, dose)

    def testGetCorrectedDoseKeepsRealDoseValues(self):
        harness = _DoseHarness(
            _Acquisition(doseInitial=1.0, dosePerFrame=2.0))

        preExp, dose = ProtMotionCorrBase._getCorrectedDose(harness)

        # firstFrame == 1 -> preExp += dose * (1 - 1) == 0
        self.assertEqual(1.0, preExp)
        self.assertEqual(2.0, dose)


class TestMotionCorrHasValidDoseRegression(unittest.TestCase):
    # Regression tests: _validate() blocks a direct launch with
    # doApplyDoseFilter on and no dose, but that validation is not
    # necessarily re-enforced on every launch path (e.g. a protocol
    # started as part of a resumed/chained workflow). Asking MotionCor2
    # to dose-weight with an unusable dose (-FmDose 0) makes it
    # silently skip writing the dose-weighted output altogether, which
    # then crashed createOutputStep looking for a file that was never
    # produced. _hasValidDose() must be checked at the point the args
    # are built and the point the output is registered, not only at
    # _validate() time.

    def testHasValidDoseIsFalseWhenDoseMissing(self):
        harness = _DoseHarness(
            _Acquisition(doseInitial=None, dosePerFrame=None))

        self.assertFalse(ProtMotionCorrBase._hasValidDose(harness))

    def testHasValidDoseIsFalseWhenDoseIsZero(self):
        harness = _DoseHarness(
            _Acquisition(doseInitial=0.0, dosePerFrame=0.0))

        self.assertFalse(ProtMotionCorrBase._hasValidDose(harness))

    def testHasValidDoseIsTrueWhenDoseIsPresent(self):
        harness = _DoseHarness(
            _Acquisition(doseInitial=1.0, dosePerFrame=2.0))

        self.assertTrue(ProtMotionCorrBase._hasValidDose(harness))

    def testHasValidDoseResolvesInputMoviesOnlyOnceAcrossManyCalls(self):
        # Regression test: _hasValidDose() is called multiple times per
        # movie from three different places (_getMcArgs,
        # createOutputStep, setMicPlotInfo). Re-resolving
        # getInputMovies().getAcquisition().getDosePerFrame() fresh on
        # every single call multiplies the exposure to any
        # inconsistency in how the input Set gets reconstructed across
        # calls - which is exactly what let the same protocol run
        # inconsistently split its movies between a dose-weighted and
        # a non-dose-weighted output. Caching after the first call
        # eliminates that exposure entirely for the rest of the run.
        harness = _DoseHarness(
            _Acquisition(doseInitial=1.0, dosePerFrame=2.0))
        harness._inputMovies.getAcquisitionCalls = 0
        realGetAcquisition = harness._inputMovies.getAcquisition

        def _countingGetAcquisition():
            harness._inputMovies.getAcquisitionCalls += 1
            return realGetAcquisition()

        harness._inputMovies.getAcquisition = _countingGetAcquisition

        for _ in range(5):
            ProtMotionCorrBase._hasValidDose(harness)

        self.assertEqual(
            1,
            harness._inputMovies.getAcquisitionCalls,
            "_hasValidDose() must resolve the input Set's Acquisition "
            "only once per protocol instance, not once per call.",
        )

    def testHasValidDoseAndGetCorrectedDoseShareTheSameCachedAcquisition(self):
        # Regression test for a real production failure: within a
        # single run, _hasValidDose() cached True (so dose-weighting
        # keeps being requested for every movie), but _getCorrectedDose
        # (called from _getMcArgs/calcFrameMotion, once per movie) kept
        # resolving getInputMovies().getAcquisition() fresh on its own
        # and occasionally got back dose=None/voltage=None for a movie
        # in the middle of the run - MotionCor2 then silently skipped
        # writing the dose-weighted output for that movie (and every
        # one after it) while still being asked for one, making
        # createOutputStep crash looking for a _DW.mrc that was never
        # produced. Both call sites must resolve the Acquisition
        # through the same cache, not independently.
        harness = _DoseHarness(
            _Acquisition(doseInitial=1.0, dosePerFrame=2.0))
        harness._inputMovies.getAcquisitionCalls = 0
        realGetAcquisition = harness._inputMovies.getAcquisition

        def _countingGetAcquisition():
            harness._inputMovies.getAcquisitionCalls += 1
            return realGetAcquisition()

        harness._inputMovies.getAcquisition = _countingGetAcquisition

        ProtMotionCorrBase._hasValidDose(harness)
        for _ in range(3):
            ProtMotionCorrBase._getCorrectedDose(harness)

        self.assertEqual(
            1,
            harness._inputMovies.getAcquisitionCalls,
            "_hasValidDose() and _getCorrectedDose() must resolve the "
            "input Set's Acquisition through the same shared cache.",
        )


class _ValidateInputMovies:
    def __init__(self, acquisition, gain=None):
        from pwem.objects import Movie

        self._acquisition = acquisition
        self._gain = gain
        self._movie = Movie()
        self._movie.setFileName("movie_001.mrc")

    def getFirstItem(self):
        return self._movie

    def getAcquisition(self):
        return self._acquisition

    def getGain(self):
        return self._gain


class _DoseValidateHarness(ProtMotionCorrBase):
    def __init__(self, acquisition, doApplyDoseFilter):
        self._inputMovies = _ValidateInputMovies(acquisition)
        self.doApplyDoseFilter = _Value(doApplyDoseFilter)
        self.alignFrame0 = _Value(1)
        self.alignFrameN = _Value(0)

    def getInputMovies(self):
        return self._inputMovies

    def _getNumberOfFrames(self):
        return 10


class TestMotionCorrMissingDoseValidationRegression(unittest.TestCase):
    # Regression test: a missing dose was briefly turned into a
    # non-blocking _warnings() notice (letting the protocol launch
    # with dose treated as 0), but that surfaced worse downstream
    # failures - MotionCor2 itself silently skips writing a
    # dose-weighted output when given -FmDose 0, so the protocol later
    # crashed in createOutputStep looking for a _DW.mrc file that was
    # never produced. Missing dose with doApplyDoseFilter on must go
    # back to blocking the launch with a clear _validate() error.

    def testMissingDoseBlocksLaunchWhenDoseFilterIsOn(self):
        harness = _DoseValidateHarness(
            _Acquisition(doseInitial=None, dosePerFrame=None),
            doApplyDoseFilter=True,
        )

        module = "motioncorr.protocols.protocol_base"
        with patch(module + ".exists", return_value=True):
            errors = ProtMotionCorrBase._validate(harness)

        self.assertEqual(1, len(errors))
        self.assertIn("dose", errors[0].lower())

    def testMissingDoseDoesNotBlockLaunchWhenDoseFilterIsOff(self):
        harness = _DoseValidateHarness(
            _Acquisition(doseInitial=None, dosePerFrame=None),
            doApplyDoseFilter=False,
        )

        module = "motioncorr.protocols.protocol_base"
        with patch(module + ".exists", return_value=True):
            errors = ProtMotionCorrBase._validate(harness)

        self.assertEqual([], errors)

    def testPresentDoseDoesNotBlockLaunch(self):
        harness = _DoseValidateHarness(
            _Acquisition(doseInitial=1.0, dosePerFrame=2.0),
            doApplyDoseFilter=True,
        )

        module = "motioncorr.protocols.protocol_base"
        with patch(module + ".exists", return_value=True):
            errors = ProtMotionCorrBase._validate(harness)

        self.assertEqual([], errors)

    def test_OpenEmptyStreamingInputDefersDoseValidationUntilMetadataArrives(self):
        class OpenEmptyInputMovies:
            def getFirstItem(self):
                return None

            def isStreamOpen(self):
                return True

            def getAcquisition(self):
                return _Acquisition(
                    doseInitial=None,
                    dosePerFrame=None,
                )

            def getGain(self):
                return None

        class Harness(ProtMotionCorrBase):
            def __init__(self):
                self._inputMovies = OpenEmptyInputMovies()
                self.doApplyDoseFilter = _Value(True)
                self.alignFrame0 = _Value(1)
                self.alignFrameN = _Value(0)

            def getInputMovies(self):
                return self._inputMovies

            def _getNumberOfFrames(self):
                raise AssertionError(
                    "Frame metadata must not be required before the "
                    "first streaming movie exists."
                )

        harness = Harness()

        module = "motioncorr.protocols.protocol_base"
        with patch(module + ".exists", return_value=True):
            errors = ProtMotionCorrBase._validate(harness)

        self.assertEqual(
            [],
            errors,
            "An open, still-empty streaming input must not be rejected "
            "as if its final acquisition metadata were already known.",
        )

    def test_HasValidDoseFallsBackToFirstMovieAcquisitionWhenSetMetadataIsIncomplete(self):
        setAcquisition = _Acquisition(
            doseInitial=None,
            dosePerFrame=None,
            voltage=None,
        )
        movieAcquisition = _Acquisition(
            doseInitial=0.5,
            dosePerFrame=1.25,
            voltage=300.0,
        )

        class FirstMovie:
            def getAcquisition(self):
                return movieAcquisition

        class InputMovies:
            def getAcquisition(self):
                return setAcquisition

            def getFirstItem(self):
                return FirstMovie()

            def getFramesRange(self):
                return [1, 10, 0]

        class Harness(ProtMotionCorrBase):
            def __init__(self):
                self._inputMovies = InputMovies()
                self._lock = threading.RLock()

            def getInputMovies(self):
                return self._inputMovies

            def _getNumberOfFrames(self):
                return 10

        harness = Harness()

        self.assertTrue(
            ProtMotionCorrBase._hasValidDose(harness),
            "A temporarily incomplete Set-level Acquisition must not hide "
            "valid acquisition metadata already present on the first Movie.",
        )
        self.assertEqual(
            (300.0, 0.5, 1.25),
            ProtMotionCorrBase._getCachedAcquisitionValues(harness),
        )


class _FrameMotionHarness:
    # Minimal harness shared by both calcFrameMotion regression tests
    # below (classic ProtMotionCorr and ProtMotionCorrNewStreaming
    # share the same calcFrameMotion shape and bug).
    def __init__(self, dose):
        self.isEER = False
        self.sRate = 1.0
        self._dose = dose

    def _getMovieShifts(self, movie):
        return [0.0, 1.0, 2.0], [0.0, 1.0, 2.0]

    def _getFramesRange(self):
        return 1, 3

    def _getCorrectedDose(self):
        return 0.0, self._dose

    def getSamplingRate(self):
        return self.sRate


class TestMotionCorrCalcFrameMotionZeroDoseRegression(unittest.TestCase):
    # Regression test: once a missing dose became a non-blocking
    # _warnings() notice instead of a hard _validate() error (see
    # TestMotionCorrMissingDoseWarningRegression above),
    # _getCorrectedDose legitimately returns dose == 0.0 for a movie
    # whose acquisition never had a dose per frame. calcFrameMotion's
    # "cutoff = (4 - preExp) // dose" had no guard for that, crashing
    # with ZeroDivisionError - in protocol_motioncorr_ns.py this
    # happens inside createOutputStep's SET OF MICROGRAPHS block,
    # which deliberately re-raises (fail-loud for output Set
    # integrity), so it crashed the whole protocol rather than just
    # that one micrograph.

    def testNewStreamingCalcFrameMotionDoesNotCrashOnZeroDose(self):
        from motioncorr.protocols.protocol_motioncorr_ns import (
            ProtMotionCorrNewStreaming,
        )

        harness = _FrameMotionHarness(dose=0.0)

        total, early, late = ProtMotionCorrNewStreaming.calcFrameMotion(
            harness, "movie_001.mrc")

        self.assertGreater(total, 0.0)
        # With no known dose, every frame is treated as "early".
        self.assertEqual(total, early)
        self.assertEqual(0.0, late)

    def testClassicCalcFrameMotionDoesNotCrashOnZeroDose(self):
        harness = _FrameMotionHarness(dose=0.0)

        total, early, late = ProtMotionCorr.calcFrameMotion(
            harness, "movie_001.mrc")

        self.assertGreater(total, 0.0)
        self.assertEqual(total, early)
        self.assertEqual(0.0, late)


if __name__ == "__main__":
    unittest.main()
