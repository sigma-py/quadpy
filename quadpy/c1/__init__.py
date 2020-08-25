from ._adaptive import integrate_adaptive
from ._chebyshev_gauss import chebyshev_gauss_1, chebyshev_gauss_2
from ._clenshaw_curtis import clenshaw_curtis
from ._fejer import fejer_1, fejer_2
from ._gauss_jacobi import gauss_jacobi
from ._gauss_kronrod import gauss_kronrod
from ._gauss_legendre import gauss_legendre
from ._gauss_lobatto import gauss_lobatto
from ._gauss_patterson import gauss_patterson
from ._gauss_radau import gauss_radau
from ._midpoint import midpoint
from ._newton_cotes import newton_cotes_closed, newton_cotes_open
from ._trapezoidal import trapezoidal

__all__ = [
    "chebyshev_gauss_1",
    "chebyshev_gauss_2",
    "clenshaw_curtis",
    "fejer_1",
    "fejer_2",
    "gauss_jacobi",
    "gauss_kronrod",
    "gauss_legendre",
    "gauss_lobatto",
    "gauss_patterson",
    "gauss_radau",
    "midpoint",
    "newton_cotes_open",
    "newton_cotes_closed",
    "trapezoidal",
    "integrate_adaptive",
]
