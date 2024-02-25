from pyworkflow.utils import weakImport

from .test_protocols_motioncor import TestMotioncorAlignMovies

with weakImport('tomo'):
    from .test_protocols_tomo import TestMotioncorTiltSeriesAlignMovies
