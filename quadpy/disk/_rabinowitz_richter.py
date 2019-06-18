# -*- coding: utf-8 -*-
#
from .helpers import _s40, _s8, _s4, _z, DiskScheme
from ..helpers import untangle, article


_citation = article(
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
    data = [
        (0.1851958765246450, _s40(0.8377170225998396)),
        (0.2930225148631698, _s40(0.3924393142315810)),
        (0.2296152967863584, _s4(0.5505609906724360)),
        (0.0387822376116376, _s8(0.4249164962326038, 0.9112013890413142)),
    ]
    points, weights = untangle(data)
    return DiskScheme("Rabinowitz-Richter 1", weights, points, 9, _citation)


def rabinowitz_richter_2():
    data = [
        (0.0043173954188430, _z()),
        (0.0734867016303473, _s40(0.9499490053854548)),
        (0.3295210136662689, _s40(0.4184300297249359)),
        (0.0046091399966757, _s4(0.8485281374238570)),
        (0.1883509796247228, _s8(0.3830079234911947, 0.7409163950514299)),
    ]
    points, weights = untangle(data)
    return DiskScheme("Rabinowitz-Richter 2", weights, points, 9, _citation)


def rabinowitz_richter_3():
    data = [
        (0.0513100527123566, _s40(0.9619017737816972)),
        (0.1136282065100473, _s40(0.7745966692414835)),
        (0.2083682752319387, _s40(0.3287526591967855)),
        (0.0790679779683282, _s8(0.4683708939890903, 0.8112421851755608)),
        (0.1269778365032246, _s8(0.3375826402485671, 0.5847102846637651)),
    ]
    points, weights = untangle(data)
    return DiskScheme("Rabinowitz-Richter 3", weights, points, 11, _citation)


def rabinowitz_richter_4():
    data = [
        (0.0478396326404247, _s40(0.9669004345445009)),
        (0.1597003917456590, _s40(0.7226054070052285)),
        (0.2016322022034297, _s40(0.3233163607428629)),
        (0.0165089733783664, _s4(0.7036534680827588)),
        (0.1801837855454157, _s4(0.4638891735186042)),
        (0.0897665889420765, _s8(0.4135214625627066, 0.8138386408455507)),
    ]
    points, weights = untangle(data)
    return DiskScheme("Rabinowitz-Richter 4", weights, points, 11, _citation)


def rabinowitz_richter_5():
    data = [
        (0.1604310638138027, _z()),
        (0.1424323658922069, _s40(0.3879803784555729)),
        (0.0556845391070962, _s40(0.9358527527678654)),
        (0.1114444717392537, _s40(0.7134059509780893)),
        (0.0449789946826613, _s4(0.6759153919798939)),
        (0.1347199228191621, _s4(0.3835039628013994)),
        (0.0316618826416774, _s8(0.3464101615137754, 0.9066277008560241)),
        (0.0963531689601313, _s8(0.7106593341863341, 0.3816598192059473)),
    ]
    points, weights = untangle(data)
    return DiskScheme("Rabinowitz-Richter 5", weights, points, 13, _citation)


def rabinowitz_richter_6():
    # Myskovskih has the same, but with analytical points and weights.
    data = [
        (0.1252902085642858, _s40(0.2528637970912295)),
        (0.1095003911263660, _s40(0.5777289284448234)),
        (0.0167126254970435, _s40(0.9897468025114907)),
        (0.0662374557963763, _s40(0.8738369566448817)),
        (0.1274283726817204, _s4(0.3754168246261542)),
        (0.0261028601843605, _s4(0.6892993807911362)),
        (0.0660009346611046, _s4(0.5976143046671968)),
        (0.0425230658266824, _s8(0.3657908004006625, 0.8830971113185893)),
        (0.0815395916164132, _s8(0.2930307227106603, 0.7074387449600663)),
    ]
    points, weights = untangle(data)
    return DiskScheme("Rabinowitz-Richter 6", weights, points, 15, _citation)
