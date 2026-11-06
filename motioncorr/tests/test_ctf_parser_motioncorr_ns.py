import tempfile
import unittest
from pathlib import Path

import numpy as np

from pwem.objects import CTFModel

from ..protocols.protocol_motioncorr_ns import ProtMotionCorrNewStreaming


class TestMotionCorrNsCtfParser(unittest.TestCase):
    def _new_protocol(self):
        # Avoid full protocol initialization; _readCtfModel does not depend on it.
        return ProtMotionCorrNewStreaming.__new__(ProtMotionCorrNewStreaming)

    def _write_file(self, path: Path, content: str):
        path.write_text(content, encoding='utf-8')

    def _assert_standardized_angle(self, ctf_angle: float, raw_angle: float, places: int = 2):
        """setStandardDefocus normalizes astigmatism angle to [0, 180)."""
        expected = raw_angle % 180.0
        self.assertAlmostEqual(ctf_angle, expected, places=places)

    def test_read_ctf_model_six_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ctf_txt = tmp / 'movie_aligned_mic_Ctf.txt'
            psd_mrc = tmp / 'movie_aligned_mic_Ctf.mrc'
            psd_mrc.touch()

            self._write_file(
                ctf_txt,
                "# header\n"
                "11744.11 11479.96 -50.57 0.10 0.06715 17.2973\n",
            )

            prot = self._new_protocol()
            ctf = prot._readCtfModel(CTFModel(), str(ctf_txt), str(psd_mrc))

            self.assertAlmostEqual(ctf.getDefocusU(), 11744.11, places=2)
            self.assertAlmostEqual(ctf.getDefocusV(), 11479.96, places=2)
            self._assert_standardized_angle(ctf.getDefocusAngle(), -50.57)
            self.assertAlmostEqual(ctf.getFitQuality(), 0.06715, places=5)
            self.assertAlmostEqual(ctf.getResolution(), 17.2973, places=4)
            self.assertAlmostEqual(ctf.getPhaseShift(), np.rad2deg(0.10), places=5)
            self.assertEqual(ctf.getPsdFile(), str(psd_mrc))

    def test_read_ctf_model_with_index_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ctf_txt = tmp / 'movie_aligned_mic_Ctf.txt'

            self._write_file(
                ctf_txt,
                "# header\n"
                "1 12153.57 11175.66 -52.43 0.00 0.09779 13.0612\n",
            )

            prot = self._new_protocol()
            ctf = prot._readCtfModel(CTFModel(), str(ctf_txt), None)

            self.assertAlmostEqual(ctf.getDefocusU(), 12153.57, places=2)
            self.assertAlmostEqual(ctf.getDefocusV(), 11175.66, places=2)
            self._assert_standardized_angle(ctf.getDefocusAngle(), -52.43)
            self.assertAlmostEqual(ctf.getFitQuality(), 0.09779, places=5)
            self.assertAlmostEqual(ctf.getResolution(), 13.0612, places=4)
            self.assertIsNone(ctf.getPhaseShift())

    def test_read_ctf_model_invalid_values(self):
        # A negative defocus, too few values or an empty file are no CTF: the
        # movie is left out instead of being published with made-up values
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            prot = self._new_protocol()
            for content in ("# header\n-1 11175.66 -52.43 0.00 0.09779 13.0612\n",
                            "1 12153.57 11175.66\n",
                            "# only a header\n"):
                ctf_txt = tmp / 'movie_aligned_mic_Ctf.txt'
                self._write_file(ctf_txt, content)
                with self.assertRaises(Exception):
                    prot._readCtfModel(CTFModel(), str(ctf_txt), None)

            with self.assertRaises(Exception):
                prot._readCtfModel(CTFModel(), str(tmp / 'missing_Ctf.txt'), None)

    def test_ctf_file_names_follow_the_micrograph(self):
        prot = self._new_protocol()
        prot._getExtraPath = lambda name: '/run/extra/' + name
        self.assertEqual(prot._getResultCtfFn('/data/movie_01.mrc.eer'),
                         '/run/extra/movie_01_aligned_mic_Ctf.txt')
        self.assertEqual(prot._getResultPsdFn('/data/movie_01.tif'),
                         '/run/extra/movie_01_aligned_mic_Ctf.mrc')


class _Flag:
    """A boolean form parameter."""

    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class TestMotionCorrNsCtfPublication(unittest.TestCase):
    """A movie whose CTF cannot be read reaches no output."""

    def test_unreadable_ctf_leaves_the_movie_out_of_every_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ctf_txt = Path(tmpdir) / 'movie_aligned_mic_Ctf.txt'
            ctf_txt.write_text("# MotionCor3 gave up\n", encoding='utf-8')

            prot = ProtMotionCorrNewStreaming.__new__(ProtMotionCorrNewStreaming)
            prot.failedMovies = []
            prot._missingResultFiles = lambda movieFName: []
            prot.estimateCtf = _Flag(True)
            prot._getResultCtfFn = lambda movieFName: str(ctf_txt)
            prot._getResultPsdFn = lambda movieFName: None

            def published(*args, **kwargs):
                raise AssertionError("nothing may be published")

            prot._getOutputMovies = published
            prot._registerMics = published
            prot._registerCtf = published

            prot.createOutputStep('/data/movie.tif', object())

            self.assertEqual(prot.failedMovies, ['/data/movie.tif'])


class TestMotionCorrCtfArguments(unittest.TestCase):
    """The CTF is estimated only when the protocol offers it and it is chosen."""

    def test_estimate_ctf_switch(self):
        from ..protocols.protocol_base import ProtMotionCorrBase

        prot = ProtMotionCorrNewStreaming.__new__(ProtMotionCorrNewStreaming)
        prot._cachedCtfAcquisitionValues = (2.7, 0.1)

        prot._estimatesCtf = lambda: False
        self.assertFalse(ProtMotionCorrBase._estimatesCtf(prot))
        self.assertEqual(prot._getCtfAcquisitionValues(), (2.7, 0.1))

        prot.estimateCtf = _Flag(True)
        self.assertTrue(ProtMotionCorrNewStreaming._estimatesCtf(prot))
        prot.estimateCtf = _Flag(False)
        self.assertFalse(ProtMotionCorrNewStreaming._estimatesCtf(prot))


if __name__ == '__main__':
    unittest.main()
