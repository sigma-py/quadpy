import pathlib

from ...helpers import article
from .._helpers import C2Scheme, _read, concat, s4, zero

source = article(
    authors=["R. Cools", "Ann Haegemans"],
    title="Another Step Forward in Searching for Cubature Formulae with a Minimal Number of Knots for the Square",
    journal="Computing",
    volume="40",
    pages="139-146",
    year="1988",
    url="https://doi.org/10.1007/BF02247942",
)

this_dir = pathlib.Path(__file__).resolve().parent


def cools_haegemans_1988_1():
    return _read(this_dir / "cools_haegemans_1988_1.json", source)


def cools_haegemans_1988_2():
    return _read(this_dir / "cools_haegemans_1988_2.json", source)
    weights, points = concat(
        s4(
            [
                0.29991838864499131666e-01,
                0.77880971155441942252e00,
                0.98348668243987226379e00,
            ],
            [
                0.38174421317083669640e-01,
                0.95729769978630736566e00,
                0.85955600564163892859e00,
            ],
            [
                0.60424923817749980681e-01,
                0.13818345986246535375e00,
                0.95892517028753485754e00,
            ],
            [
                0.77492738533105339358e-01,
                0.94132722587292523695e00,
                0.39073621612946100068e00,
            ],
            [
                0.11884466730059560108e00,
                0.47580862521827590507e00,
                0.85007667369974857597e00,
            ],
            [
                0.12976355037000271129e00,
                0.75580535657208143627e00,
                0.64782163718701073204e00,
            ],
            [
                0.21334158145718938943e00,
                0.69625007849174941396e00,
                0.70741508996444936217e-01,
            ],
            [
                0.25687074948196783651e00,
                0.34271655604040678941e00,
                0.40930456169403884330e00,
            ],
        ),
        zero(0.30038211543122536139e00),
    )
    weights /= 4
    return C2Scheme("Cools-Haegemans 1988-2", weights, points, 13, source)
