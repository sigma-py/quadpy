from ..helpers import article
from ._helpers import C2Scheme, concat, symm_r0, symm_s, symm_s_t, zero

source = article(
    authors=["Philip Rabinowitz", "Nira Richter"],
    title="Perfectly Symmetric Two-Dimensional Integration Formulas with Minimal Numbers of Points",
    journal="Mathematics of Computation",
    volume="23",
    number="108",
    month="oct",
    year="1969",
    pages="765-779",
    url="https://doi.org/10.2307/2004962",
)


def rabinowitz_richter_1():
    name = "Rabinowitz-Richter 1"
    degree = 9
    weights, points = concat(
        symm_r0(
            [0.0716134247098111, 0.9845398119422523],
            [0.4540903525515453, 0.4888863428423724],
        ),
        symm_s([0.0427846154667780, 0.9395672874215217]),
        symm_s_t([0.2157558036359328, 0.8367103250239890, 0.5073767736746132]),
    )
    weights /= 4
    return C2Scheme(name, weights, points, degree, source)


def rabinowitz_richter_2():
    name = "Rabinowitz-Richter 2"
    degree = 11
    weights, points = concat(
        zero(0.3653795255859022),
        symm_r0(
            [0.2442720577517539, 0.7697990683966493],
            [0.0277561655642043, 1.044402915409813],
        ),
        symm_s(
            [0.3089930361337136, 0.4134919534491139],
            [0.0342651038512293, 0.9357870124405403],
        ),
        symm_s_t([0.1466843776513117, 0.5756535958404649, 0.8830255085256902]),
    )
    weights /= 4
    return C2Scheme(name, weights, points, degree, source)


def rabinowitz_richter_3():
    name = "Rabinowitz-Richter 3"
    degree = 11
    weights, points = concat(
        symm_r0(
            [0.0176679598882646, 0.8989737240828844],
            [0.2322248008989674, 0.7632367891419969],
        ),
        symm_s(
            [0.0715516745178401, 0.8949648832822285],
            [0.2192868905662522, 0.6322452037101431],
            [0.2965842326220580, 0.2797353125538562],
        ),
        symm_s_t([0.0813422207533089, 0.9602661668053869, 0.4347413023856830]),
    )
    weights /= 4
    return C2Scheme(name, weights, points, degree, source)


def rabinowitz_richter_4():
    name = "Rabinowitz-Richter 4"
    degree = 13
    weights, points = concat(
        zero(0.2995235559387052),
        symm_r0(
            [0.0331100668669073, 0.9909890363004588],
            [0.1802214941550577, 0.6283940712305196],
        ),
        symm_s(
            [0.0391672789603492, 0.9194861553393097],
            [0.1387748348777338, 0.6973201917871096],
            [0.2268881207335663, 0.3805687186904865],
        ),
        symm_s_t(
            [0.0365739576550950, 0.9708504361720127, 0.6390348393207334],
            [0.1169047000557597, 0.8623637916722844, 0.3162277660168378],
        ),
    )
    weights /= 4
    return C2Scheme(name, weights, points, degree, source)


def rabinowitz_richter_5():
    name = "Rabinowitz-Richter 5"
    degree = 15
    weights, points = concat(
        symm_r0(
            [-0.40980941939297e-5, 1.315797935069747],
            [0.0414134647558384, 0.9796158388578564],
            [0.1837583771750436, 0.6375456844500517],
        ),
        symm_s(
            [0.0280217865486269, 0.9346799288936658],
            [0.0948146979601645, 0.7662665721615083],
            [0.1688054053337613, 0.5138362475917853],
            [0.1898474000367674, 0.2211895845055072],
        ),
        symm_s_t(
            [0.0331477474104121, 0.9769495664551867, 0.6375975639376926],
            [0.1135237357315838, 0.8607803779721935, 0.3368688874716777],
        ),
    )
    weights /= 4
    return C2Scheme(name, weights, points, degree, source)


def rabinowitz_richter_6():
    name = "Rabinowitz-Richter 6"
    degree = 15
    weights, points = concat(
        symm_r0(
            [0.0301245207981210, 0.9915377816777667],
            [0.0871146840209092, 0.8020163879230440],
            [0.1250080294351494, 0.5648674875232742],
        ),
        symm_s(
            [0.0267651407861666, 0.9354392392539896],
            [0.0959651863624437, 0.7624563338825799],
            [0.1750832998343375, 0.2156164241427213],
        ),
        symm_s_t(
            [0.0283136372033274, 0.9769662659711761, 0.6684480048977932],
            [0.0866414716025093, 0.8937128379503403, 0.3735205277617582],
            [0.1150144605755996, 0.6122485619312083, 0.4078983303613935],
        ),
    )
    weights /= 4
    return C2Scheme(name, weights, points, degree, source)
