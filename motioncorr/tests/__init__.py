from pyworkflow.utils import weakImport

from .test_protocols_motioncor2 import TestMotioncor2AlignMovies

with weakImport('tomo'):
    from .test_protocols_tomo import TestMotioncor2TiltSeriesAlignMovies
