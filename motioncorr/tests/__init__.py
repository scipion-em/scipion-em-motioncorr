from pyworkflow.utils import weakImport

from .test_protocols_motioncor_ns import TestMotioncorNSAlignMovies

with weakImport('tomo'):
    from .test_protocols_tomo import TestMotioncorTiltSeriesM
