
from .test_protocols_motioncor2 import TestMotioncor2AlignMovies

try:
    from .test_protocols_tomo import TestMotioncor2TiltSeriesAlignMovies
except Exception as e:
    if not "'tomo'" in str(e):
        raise e
