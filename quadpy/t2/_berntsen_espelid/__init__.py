import pathlib
from mpmath import mp

from ...helpers import article, techreport
from .._helpers import T2Scheme, concat, s1, s2, s3, _read

source = techreport(
    authors=["J. Berntsen", "T.O. Espelid"],
    title="Degree 13 symmetric quadrature rules for the triangle",
    # Reports in Informatics,
    institution="Dept. of Informatics, University of Bergen",
    year="1990",
)

# This first scheme was published separately as
c2 = article(
    authors=["J. Berntsen", "T.O. Espelid"],
    title="Algorithm 706: DCUTRI: An Algorithm for Adaptive Cubature over a Collection of Triangles",
    journal="ACM Trans. Math. Softw.",
    month="sep",
    year="1992",
    url="https://dl.acm.org/source.cfm?id=131772",
)

this_dir = pathlib.Path(__file__).resolve().parent


def dcutri(mpmath=False):
    out = berntsen_espelid_1(mpmath)
    out.source = c2
    return out


def berntsen_espelid_1(mpmath=False):
    return _read(this_dir / "berntsen_espelid_1.json", source)
    flt = mp.mpf if mpmath else float

    mp.dps = 30
    weights, points = concat(
        s3(flt("0.051739766065744133555179145422")),
        s2(
            [
                flt("0.008007799555564801597804123460"),
                flt("0.024862168537947217274823955239"),
            ],
            [
                flt("0.046868898981821644823226732071"),
                flt("0.414192542538082326221847602214"),
            ],
            [
                flt("0.046590940183976487960361770070"),
                flt("0.230293878161404779868453507244"),
            ],
            [
                flt("0.031016943313796381407646220131"),
                flt("0.113919981661733719124857214943"),
            ],
            [
                flt("0.010791612736631273623178240136"),
                flt("0.495457300025082323058213517632"),
            ],
            [
                flt("0.032195534242431618819414482205"),
                flt("0.468861354847056503251458179727"),
            ],
        ),
        s1(
            [
                flt("0.015445834210701583817692900053"),
                flt("0.022076289653624405142446876931"),
                flt("0.851306504174348550389457672223"),
            ],
            [
                flt("0.017822989923178661888748319485"),
                flt("0.018620522802520968955913511549"),
                flt("0.689441970728591295496647976487"),
            ],
            [
                flt("0.037038683681384627918546472190"),
                flt("0.096506481292159228736516560903"),
                flt("0.635867859433372768286976979827"),
            ],
        ),
    )
    return T2Scheme(
        "Berntsen-Espelid 1 (DCUTRI)", weights, points, 13, source, 1.516e-12
    )


def berntsen_espelid_2(mpmath=False):
    return _read(this_dir / "berntsen_espelid_2.json", source)


def berntsen_espelid_3(mpmath=False):
    return _read(this_dir / "berntsen_espelid_3.json", source)


def berntsen_espelid_4(mpmath=False):
    return _read(this_dir / "berntsen_espelid_4.json", source)
