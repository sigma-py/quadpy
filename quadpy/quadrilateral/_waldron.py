import numpy

from ..helpers import techreport
from ._helpers import QuadrilateralScheme

citation = techreport(
    authors=["Shayne Waldron"],
    title="Symmetries of linear functionals",
    institution="University of Wisconsin-Madison, Center for Mathematical Sciences",
    note="CMS Technical Summary Report",
    month="oct",
    year="1994",
    url="http://ftp.cs.wisc.edu/Approx/symmetries.pdf",
)


def waldron(r, alpha):
    assert r ** 2 >= 1 / 3

    R = r / numpy.sqrt(3 * r ** 2 - 1)

    beta = alpha + numpy.pi / 2

    points = numpy.array(
        [
            [+r * numpy.cos(alpha), +r * numpy.sin(alpha)],
            [-r * numpy.cos(alpha), -r * numpy.sin(alpha)],
            [+R * numpy.cos(beta), +R * numpy.sin(beta)],
            [-R * numpy.cos(beta), -R * numpy.sin(beta)],
        ]
    )

    weights = numpy.array(
        [2 / 3 / r ** 2, 2 / 3 / r ** 2, 2 / 3 / R ** 2, 2 / 3 / R ** 2]
    )
    return QuadrilateralScheme("Waldron", weights, points, 3, citation)
