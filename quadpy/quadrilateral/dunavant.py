# -*- coding: utf-8 -*-
#
import numpy


class Dunavant(object):
    '''
    D.A. Dunavant,
    Economical symmetrical quadrature rules for complete polynomials over a
    square domain,
    Numerical Methods in Engineering, Volume 21, Issue 10, October 1985,
    Pages 1777–1784.

    Abstract:
    It is of interest in numerical analysis to develop symmetrical quadrature
    rules for integration of complete polynomial functions over a square domain
    with minimum computational effort. Gaussian product quadrature rules
    integrate such functions with maximum effort. Symmetrical quadrature rules
    are developed and presented for integration of complete polynomial
    functions up to 21st order with minimum computational effort.
    '''
    def __init__(self, index):
        if index == 0:
            self.degree = 1
            self.weights = numpy.array([4.0])
            self.points = numpy.array([[0.0, 0.0]])
        elif index == 1:
            self.degree = 3
            self.weights = numpy.full(4, 1.0)
            self.points = _symm_s(numpy.sqrt(1.0/3.0))
        elif index == 2:
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(4, 40.0/49.0),
                numpy.full(4, 9.0/49.0),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(numpy.sqrt(7.0/15.0)),
                _symm_s(numpy.sqrt(7.0/9.0)),
                ])
        elif index == 3:
            self.degree = 7
            self.weights = numpy.concatenate([
                numpy.full(4, 98.0 / 405.0),
                numpy.full(4, 0.237431774690630),
                numpy.full(4, 0.520592916667394),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(numpy.sqrt(6.0/7.0)),
                _symm_s(0.805979782918599),
                _symm_s(0.380554433208316),
                ])
        elif index == 4:
            self.degree = 9
            self.weights = numpy.concatenate([
                numpy.full(4, 0.018475842507491),
                numpy.full(4, 0.390052939160735),
                numpy.full(4, 0.083095178026482),
                numpy.full(8, 0.254188020152646),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(1.121225763866564),
                _symm_r_0(0.451773049920657),
                _symm_s(0.891849420851512),
                _symm_s_t(0.824396370749276, 0.411623426336542),
                ])
        elif index == 5:
            self.degree = 11
            self.weights = numpy.concatenate([
                numpy.full(1, 0.365379525585903),
                numpy.full(4, 0.027756165564204),
                numpy.full(4, 0.244272057751754),
                numpy.full(4, 0.034265103851229),
                numpy.full(4, 0.308993036133713),
                numpy.full(8, 0.146684377651312),
                ])
            self.points = numpy.concatenate([
                _c(),
                _symm_r_0(1.044402915409813),
                _symm_r_0(0.769799068396649),
                _symm_s(0.935787012440540),
                _symm_s(0.413491953449114),
                _symm_s_t(0.883025508525690, 0.575653595840465),
                ])
        elif index == 6:
            self.degree = 13
            self.weights = numpy.concatenate([
                numpy.full(4, 0.005656169693764),
                numpy.full(4, 0.192443867470396),
                numpy.full(4, 0.005166832979773),
                numpy.full(4, 0.200302559622138),
                numpy.full(4, 0.228125175912536),
                numpy.full(8, 0.117496926974491),
                numpy.full(8, 0.066655770186205),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(1.086056158573971),
                _symm_r_0(0.658208197042585),
                _symm_s(1.001300602991729),
                _symm_s(0.584636168775946),
                _symm_s(0.246795612720261),
                _symm_s_t(0.900258815287201, 0.304720678579870),
                _symm_s_t(0.929866705560780, 0.745052720131169),
                ])
        elif index == 7:
            self.degree = 15
            self.weights = numpy.concatenate([
                numpy.full(1, -0.001768979827207),
                numpy.full(4, 0.012816726617512),
                numpy.full(4, 0.119897873101347),
                numpy.full(4, 0.210885452208801),
                numpy.full(4, 0.006392720128215),
                numpy.full(4, 0.104415680788580),
                numpy.full(8, 0.168053047203816),
                numpy.full(8, 0.076169694452294),
                numpy.full(8, 0.028794154400064),
                ])
            self.points = numpy.concatenate([
                _c(),
                _symm_r_0(1.027314357719367),
                _symm_r_0(0.856766776147643),
                _symm_r_0(0.327332998189723),
                _symm_s(0.967223740028505),
                _symm_s(0.732168901749711),
                _symm_s_t(0.621974427996805, 0.321696694921009),
                _symm_s_t(0.928618480068352, 0.455124178121179),
                _symm_s_t(0.960457474887516, 0.809863684081217),
                ])
        elif index == 8:
            self.degree = 17
            self.weights = numpy.concatenate([
                numpy.full(4, 0.020614915919991),
                numpy.full(4, 0.128025716179910),
                numpy.full(4, 0.005511739534032),
                numpy.full(4, 0.039207712457142),
                numpy.full(4, 0.076396945079863),
                numpy.full(8, 0.141513729949972),
                numpy.full(8, 0.083903279363798),
                numpy.full(8, 0.060394163649685),
                numpy.full(8, 0.057387752969213),
                numpy.full(8, 0.021922559481864),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(0.989353074512600),
                _symm_r_0(0.376285207157973),
                _symm_s(0.978848279262233),
                _symm_s(0.885794729164116),
                _symm_s(0.171756123838348),
                _symm_s_t(0.590499273806002, 0.319505036634574),
                _symm_s_t(0.799079131916863, 0.597972451929457),
                _symm_s_t(0.803743962958745, 0.058344481776551),
                _symm_s_t(0.936506276127495, 0.347386316166203),
                _symm_s_t(0.981321179805452, 0.706000287798646),
                ])
        elif index == 9:
            # TODO the article claims 19
            self.degree = 16
            self.weights = numpy.concatenate([
                numpy.full(4, 0.038205406871462),
                numpy.full(4, 0.135368502976521),
                numpy.full(4, 0.005773708558664),
                numpy.full(4, 0.067460759759473),
                numpy.full(4, 0.140899115227892),
                numpy.full(4, 0.047466627685662),
                numpy.full(8, 0.078619467342982),
                numpy.full(8, 0.094979169511394),
                numpy.full(8, 0.022331162356015),
                numpy.full(8, 0.055594877793785),
                numpy.full(8, 0.006049054506376),
                numpy.full(8, 0.024839207949609),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(0.943962831808239),
                _symm_r_0(0.536918434376013),
                _symm_s(0.973981076394170),
                _symm_s(0.742995535327609),
                _symm_s(0.285010052188916),
                _symm_s(0.068354569272491),
                _symm_s_t(0.802952004398543, 0.203345534163332),
                _symm_s_t(0.634244672807882, 0.426572172992877),
                _symm_s_t(0.978350706908227, 0.295830776620995),
                _symm_s_t(0.901672714410389, 0.541983037327871),
                _symm_s_t(1.007018449383116, 0.669414798783936),
                _symm_s_t(0.945161453573471, 0.829501421477824),
                ])
        else:
            assert index == 10
            # TODO the article claims 21
            self.degree = 20
            self.weights = numpy.concatenate([
                numpy.full(4, 0.019503841092684),
                numpy.full(4, 0.089012127744268),
                numpy.full(4, 0.114568584702749),
                numpy.full(4, 0.007463627359106),
                numpy.full(4, 0.050585943594705),
                numpy.full(4, 0.074613865184212),
                numpy.full(8, 0.023501091310143),
                numpy.full(8, 0.011588562644144),
                numpy.full(8, 0.023073245798171),
                numpy.full(8, 0.001570221774472),
                numpy.full(8, 0.049102258016277),
                numpy.full(8, 0.042512352239126),
                numpy.full(8, 0.067270936863160),
                numpy.full(8, 0.103507336515645),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(0.980883148832881),
                _symm_r_0(0.678152700336576),
                _symm_r_0(0.240599282275864),
                _symm_s(0.965176994929162),
                _symm_s(0.749698539312765),
                _symm_s(0.568983925500818),
                _symm_s_t(0.971086142843168, 0.355832132274584),
                _symm_s_t(0.983453947854968, 0.645588139196562),
                _symm_s_t(0.933927707027213, 0.821920249234369),
                _symm_s_t(1.014086498915039, 0.862185099566557),
                _symm_s_t(0.877914842155496, 0.168914072450263),
                _symm_s_t(0.882246882640128, 0.568113580166780),
                _symm_s_t(0.741324453314596, 0.371360260002223),
                _symm_s_t(0.469570217710647, 0.237333359193547),
                ])

        return


def _c():
    return numpy.array([[0.0, 0.0]])


def _symm_r_0(r):
    return numpy.array([
        [+r, 0.0],
        [-r, 0.0],
        [0.0, +r],
        [0.0, -r],
        ])


def _symm_s(s):
    return numpy.array([
        [+s, +s],
        [-s, +s],
        [+s, -s],
        [-s, -s],
        ])


def _symm_s_t(s, t):
    return numpy.array([
        [+s, +t],
        [-s, +t],
        [+s, -t],
        [-s, -t],
        [+t, +s],
        [-t, +s],
        [+t, -s],
        [-t, -s],
        ])
