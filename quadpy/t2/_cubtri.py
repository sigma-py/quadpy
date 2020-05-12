from ..helpers import article
from ._helpers import T2Scheme, concat, s1, s2, s3

citation = article(
    authors=["D.P. Laurie"],
    title="Algorithm 584: CUBTRI: Automatic Cubature over a Triangle",
    journal="ACM Trans. Math. Softw.",
    month="jun",
    year="1982",
    url="https://doi.org/10.1145/355993.356001",
)


def cubtri():
    """
    Se also
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    """
    weights, points = concat(
        s3(0.0378610912003147),
        s2(
            [0.0376204254131829, 0.1012865073234563],
            [0.0783573522441174, 0.4701420641051151],
            [0.1162714796569659, 0.2321023267750504],
            [0.0134442673751655, 0.0294808608844396],
        ),
        s1([0.0375097224552317, 0.7384168123405100, 0.2321023267750504]),
    )
    return T2Scheme("CUBTRI", weights, points, 8, citation)
