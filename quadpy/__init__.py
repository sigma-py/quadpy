# -*- coding: utf-8 -*-
#
from __future__ import print_function

from .__about__ import (
    __author__,
    __email__,
    __copyright__,
    __credits__,
    __license__,
    __version__,
    __maintainer__,
    __status__,
    )

from . import helpers

from . import ball
from . import circle
from . import disk
from . import e1r
from . import e2r
from . import e3r
from . import enr
from . import e1r2
from . import e2r2
from . import e3r2
from . import enr2
from . import hexahedron
from . import line_segment
from . import nball
from . import ncube
from . import nsimplex
from . import pyramid
from . import quadrilateral
from . import sphere
from . import triangle
from . import tetrahedron
from . import wedge

from . import tools

try:
    import pipdate
except ImportError:
    pass
else:
    if pipdate.needs_checking(__name__):
        print(pipdate.check(__name__, __version__), end='')
