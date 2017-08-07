# -*- coding: utf-8 -*-
#
from __future__ import print_function

import pipdated

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

from . import circle
from . import disk
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

if pipdated.needs_checking(__name__):
    print(pipdated.check(__name__, __version__))
