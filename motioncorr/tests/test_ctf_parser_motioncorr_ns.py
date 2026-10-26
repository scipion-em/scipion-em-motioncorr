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
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ctf_txt = tmp / 'movie_aligned_mic_Ctf.txt'

            self._write_file(
                ctf_txt,
                "# header\n"
                "-1 11175.66 -52.43 0.00 0.09779 13.0612\n",
            )

            prot = self._new_protocol()
            ctf = prot._readCtfModel(CTFModel(), str(ctf_txt), None)

            self.assertEqual(ctf.getDefocusU(), -999)
            self.assertEqual(ctf.getDefocusV(), -1)
            self.assertEqual(ctf.getDefocusAngle(), -999)
            self.assertEqual(ctf.getFitQuality(), -999)
            self.assertEqual(ctf.getResolution(), -999)


if __name__ == '__main__':
    unittest.main()
