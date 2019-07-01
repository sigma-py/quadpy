# -*- coding: utf-8 -*-
#
from mpmath import mp

from ..helpers import article, techreport
from ._helpers import TriangleScheme, concat, s1, s2, s3

citation = techreport(
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
    url="https://dl.acm.org/citation.cfm?id=131772",
)


def dcutri(mpmath=False):
    out = berntsen_espelid_1(mpmath)
    out.citation = c2
    return out


def berntsen_espelid_1(mpmath=False):
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
    return TriangleScheme("Berntsen-Espelid 1 (DCUTRI)", weights, points, 13, citation)


def berntsen_espelid_2(mpmath=False):
    flt = mp.mpf if mpmath else float

    mp.dps = 30
    weights, points = concat(
        s3(flt("0.058696079612719031799193912788")),
        s2(
            [
                flt("0.007850768296100080327451819370"),
                flt("0.024607188643230218187849951620"),
            ],
            [
                flt("0.050668953175886963421095258589"),
                flt("0.420308753101194683716920517937"),
            ],
            [
                flt("0.050080326090509066160747699962"),
                flt("0.227900255506160619646298949779"),
            ],
            [
                flt("0.031647114592298319035326893473"),
                flt("0.116213058883517905247155308064"),
            ],
            [flt("0.005356903791090860889118181848"), flt("0.5")],
            [
                flt("0.031492563075968795690055730726"),
                flt("0.476602980049079152951254192421"),
            ],
        ),
        s1(
            [
                flt("0.015802532215260751359123743555"),
                flt("0.022797894538248612547720754462"),
                flt("0.851775587145410469734660000132"),
            ],
            [
                flt("0.015981637780928405322919308674"),
                flt("0.016275770991088540943703616092"),
                flt("0.692797317566660854594116271938"),
            ],
            [
                flt("0.036551502224097295256193503655"),
                flt("0.089733060451605359079629076100"),
                flt("0.637955883864209538412552781228"),
            ],
        ),
    )
    return TriangleScheme("Berntsen-Espelid 2", weights, points, 13, citation)


def berntsen_espelid_3(mpmath=False):
    flt = mp.mpf if mpmath else float

    mp.dps = 30
    weights, points = concat(
        s2(
            [flt("-4.438917939249711e-15"), flt("-1.097321247106281159287766916114")],
            [
                flt("0.023875084055169335843543623613"),
                flt("0.488287850733405315708960134736"),
            ],
            [
                flt("0.063189783598782833129430995388"),
                flt("0.271000295524474716503595027679"),
            ],
            [
                flt("0.008045069816524589599830031859"),
                flt("0.024788431033661361058352074973"),
            ],
            [
                flt("0.027856097113552204952023523591"),
                flt("0.107120353118147709346761786284"),
            ],
            [
                flt("0.050685061067025767745642589150"),
                flt("0.440323874478061332339068546065"),
            ],
        ),
        s1(
            [
                flt("0.014867088321983380610493967543"),
                flt("0.020821520846631616958730687380"),
                flt("0.850459062644356742678494398953"),
            ],
            [
                flt("0.021575699816275772518477875728"),
                flt("0.022919482804812809947480096117"),
                flt("0.683758575887968213394629723103"),
            ],
            [
                flt("0.043398330702882367361429063273"),
                flt("0.115458022821994138042223116054"),
                flt("0.631364930935447484201224031403"),
            ],
        ),
    )
    return TriangleScheme("Berntsen-Espelid 3", weights, points, 13, citation)


def berntsen_espelid_4(mpmath=False):
    flt = mp.mpf if mpmath else float

    mp.dps = 30
    weights, points = concat(
        s3(flt("0.055141401445961668095892272765")),
        s2(
            [flt("0.000011142520455322162070507537"), flt("0")],
            [
                flt("0.008019330681470505488363363198"),
                flt("0.024978640633391274114293084881"),
            ],
            [
                flt("0.033429216779221783453803543232"),
                flt("0.474489920436516855163277733910"),
            ],
            [
                flt("0.046966588930899169431852266167"),
                flt("0.230836272600280459320993940175"),
            ],
            [
                flt("0.031079169485602998741093672276"),
                flt("0.114080598593243463483923394518"),
            ],
            [
                flt("0.048947942555161210000640851464"),
                flt("0.417965185286509715766771174230"),
            ],
            [flt("0.005884459601338707440236321752"), flt("0.5")],
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
                flt("0.635867859433872768286976979827"),
            ],
        ),
    )
    return TriangleScheme("Berntsen-Espelid 4", weights, points, 13, citation)
