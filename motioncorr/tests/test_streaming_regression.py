# **************************************************************************
# *
# * Regression tests for MotionCorr streaming/resume behaviour.
# *
# **************************************************************************

import unittest

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
    def __init__(self, doseInitial=None, dosePerFrame=None):
        self._doseInitial = doseInitial
        self._dosePerFrame = dosePerFrame

    def getDoseInitial(self):
        return self._doseInitial

    def getDosePerFrame(self):
        return self._dosePerFrame


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


class _DoseValidateHarness(ProtMotionCorrBase):
    def __init__(self, acquisition, doApplyDoseFilter):
        self._inputMovies = _DoseInputMovies(acquisition)
        self.doApplyDoseFilter = _Value(doApplyDoseFilter)

    def getInputMovies(self):
        return self._inputMovies


class TestMotionCorrMissingDoseWarningRegression(unittest.TestCase):
    # A missing dose no longer blocks the protocol from being
    # launched at all (_getCorrectedDose now degrades gracefully to
    # 0.0) - it must instead surface as a non-blocking _warnings()
    # message the user can approve past, not a hard _validate() error.

    def testMissingDoseIsAWarningNotAValidationError(self):
        harness = _DoseValidateHarness(
            _Acquisition(doseInitial=None, dosePerFrame=None),
            doApplyDoseFilter=True,
        )

        warnings = ProtMotionCorrBase._warnings(harness)

        self.assertEqual(1, len(warnings))
        self.assertIn("dose", warnings[0].lower())

    def testNoWarningWhenDoseFilterIsOff(self):
        harness = _DoseValidateHarness(
            _Acquisition(doseInitial=None, dosePerFrame=None),
            doApplyDoseFilter=False,
        )

        warnings = ProtMotionCorrBase._warnings(harness)

        self.assertEqual([], warnings)

    def testNoWarningWhenDoseIsPresent(self):
        harness = _DoseValidateHarness(
            _Acquisition(doseInitial=1.0, dosePerFrame=2.0),
            doApplyDoseFilter=True,
        )

        warnings = ProtMotionCorrBase._warnings(harness)

        self.assertEqual([], warnings)


if __name__ == "__main__":
    unittest.main()
