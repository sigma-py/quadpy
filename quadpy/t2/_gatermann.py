from ..helpers import article
from ._helpers import T2Scheme, rot_ab

source = article(
    authors=["Karin Gatermann"],
    title="The construction of symmetric cubature formulas for the square and the triangle",
    journal="Computing",
    month="sep",
    year="1988",
    volume="40",
    number="3",
    pages="229–240",
    url="https://doi.org/10.1007/BF02251251",
)


def gatermann():
    weights, points = rot_ab(
        [0.2651702815743450e-01, 0.6238226509439084e-01, 0.6751786707392436e-01],
        [0.4388140871444811e-01, 0.5522545665692000e-01, 0.3215024938520156],
        [0.2877504278497528e-01, 0.3432430294509488e-01, 0.6609491961867980],
        [0.6749318700980879e-01, 0.5158423343536001, 0.2777161669764050],
    )
    weights *= 2
    return T2Scheme("Gatermann", weights, points, 7, source, 2.004e-13)
