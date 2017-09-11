# -*- coding: utf-8 -*-
#
from .chebyshev_gauss import ChebyshevGauss1, ChebyshevGauss2
from .clenshaw_curtis import ClenshawCurtis
from .fejer import Fejer1, Fejer2
from .gauss_hermite import GaussHermite
from .gauss_kronrod import GaussKronrod
from .gauss_laguerre import GaussLaguerre
from .gauss_legendre import GaussLegendre
from .gauss_lobatto import GaussLobatto
from .gauss_patterson import GaussPatterson
from .gauss_radau import GaussRadau
from .midpoint import Midpoint
from .newton_cotes import NewtonCotesOpen, NewtonCotesClosed
from .trapezoidal import Trapezoidal

from .tanh_sinh import tanh_sinh_quadrature

# import pylint: disable=wildcard-import
from .tools import *
