from ..helpers import article
from ._hammer_wymore import hammer_wymore
from ._helpers import C3Scheme, register

source = article(
    authors=["V.L.N. Sarma", "A. H. Stroud"],
    title="Eberlein Measure and Mechanical Quadrature Formulae. II. Numerical Results",
    journal="Mathematics of Computation",
    volume="23",
    number="108",
    month="oct",
    year="1969",
    pages="781-784",
    url="https://doi.org/10.2307/2004963",
)


def sarma_stroud():
    # Hammer-Wymore is a one-parameter family of schemes, and the parameter lambda is
    # chosen to minimize the standard deviation of Sarma's error functional. The
    # particular value of lambda is not explicitly given in the article, but computed
    # from the specified values. Note that it is only given in single precision.
    hw = hammer_wymore(lmbda=1.0329785305)
    return C3Scheme("Sarma-Stroud", hw.symmetry_data, hw.degree, source)


register([sarma_stroud])
