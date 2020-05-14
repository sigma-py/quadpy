from ..helpers import article
from ._helpers import C2Scheme, concat, pm, pm2

source = article(
    authors=["Robert Piessens", "Ann Haegemans"],
    title="Cubature Formulas of Degree Nine for Symmetric Planar Regions",
    journal="Mathematics of Computation",
    volume="29",
    number="11",
    month="jul",
    year="1975",
    pages="810-815",
    url="https://doi.org/10.2307/2005291",
)


def piessens_haegemans_1():
    weights, points = concat(
        pm2(
            [0.68416522462309305679e-1, 0.87980721399752853896, 0.92797961509268528861],
            [0.27903384209687301395, 0.50445910315479838456, 0.75347199103161505380],
            [0.16806533822999587126, 0.91531235408227324183, 0.42299357094876513066],
        ),
        pm(
            [0.40927359555433144329, 0.57882826011929170546, 0],
            [0.10648011781560231854, 0, 0.97700090158004246059],
            [0.45321488105170985638, 0, 0.39364057271848893512],
        ),
    )
    weights /= 4
    return C2Scheme("Piessens-Haegemans 1", weights, points, 9, source)


def piessens_haegemans_2():
    weights, points = concat(
        pm2(
            [0.42853317248897088536e-1, 0.93742666622066710914, 0.94145119299928430974],
            [0.25788406360659644304, 0.57077001686857404415, 0.79214654516847247531],
            [0.19397744037003970872, 0.89774224179848572970, 0.40001733897633692860],
        ),
        pm(
            [0.45212398131214854997, 0.49471787965159623409, 0],
            [0.10243215270991495821, 0, 0.98085697194664054422],
            [0.45601422352687001122, 0, 0.48311469619727965642],
        ),
    )
    weights /= 4
    return C2Scheme("Piessens-Haegemans 2", weights, points, 9, source)
