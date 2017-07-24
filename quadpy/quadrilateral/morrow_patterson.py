# -*- coding: utf-8 -*-
#
import warnings

import numpy


class MorrowPatterson(object):
    '''
    C.R. Morrow and T.N.L. Patterson,
    The Construction of Algebraic Cubature Formulae by the Distribution of
    Nodes Along Selected Lines,
    SIAM J. Numer. Anal., 22(6), 1178–1190, 13 pages,
    <https://doi.org/10.1137/0722071>.

    Abstract:
    A new technique for developing general cubature rules is discussed. By
    distributing the nodes strategically along selected lines the solution is
    reduced to dealing with a series of one-dimensional moment problems.
    Several analytical and numerical illustrative examples are provided. The
    extension to higher dimensions is also considered.
    '''
    def __init__(self, index):
        warnings.warn(
            'The Morrow-Patterson schemes are only single-precision.'
            )
        self.name = 'MP({})'.format(index)
        if index == 1:
            self.degree = 11
            self.weights = numpy.concatenate([
                numpy.full(4, 0.1627661292),
                numpy.full(4, 0.3051478053),
                numpy.full(4, 0.0140154409),
                numpy.full(4, 0.2089530426),
                numpy.full(2, 0.2755861791),
                numpy.full(4, -0.1118414644e-03),
                numpy.full(4, 0.0541846008),
                numpy.full(4, 0.1172517331),
                ])
            self.points = numpy.concatenate([
                _s4(0.2386191861, 0.8611363116),
                _s4(0.2386191861, 0.3399810436),
                _s4(0.6612093865, 1.112553370),
                _s4(0.6612093865, 0.7017715829),
                _s20(0.6612093865),
                _s4(0.9324695142, 1.532186887),
                _s4(0.9324695142, 0.8829273433),
                _s4(0.9324695142, 0.3592258245),
                ])
        else:
            assert index == 2
            # The article claims degree 31.
            # TODO check for errors
            self.degree = 20
            self.weights = numpy.concatenate([
                numpy.full(4, 0.153974824894e-1),
                numpy.full(4, 0.342239043212e-1),
                numpy.full(4, 0.493728555246e-1),
                numpy.full(4, 0.591743444190e-1),
                numpy.full(2, 0.625640474012e-1),
                numpy.full(4, 0.219687383042e-1),
                numpy.full(4, 0.422241228318e-1),
                numpy.full(4, 0.557059403767e-1),
                numpy.full(4, 0.626761415564e-1),
                numpy.full(4, 0.284719756769e-4),
                numpy.full(4, 0.107864593517e-1),
                numpy.full(4, 0.288090801678e-1),
                numpy.full(4, 0.449057826672e-1),
                numpy.full(4, 0.552567682967e-1),
                numpy.full(2, 0.587968625942e-1),
                numpy.full(4, -0.238570094309e-8),
                numpy.full(4, 0.152203943218e-2),
                numpy.full(4, 0.186940270679e-1),
                numpy.full(4, 0.334515581891e-1),
                numpy.full(4, 0.448929333192e-1),
                numpy.full(4, 0.510354307986e-1),
                numpy.full(4, 0.941251999907e-11),
                numpy.full(4, 0.839123482817e-2),
                numpy.full(4, 0.215480412155e-1),
                numpy.full(4, 0.328880560591e-1),
                numpy.full(4, 0.403387733018e-1),
                numpy.full(2, 0.429295886871e-1),
                numpy.full(4, -0.192849256780e-5),
                numpy.full(4, -0.211309155229e-12),
                numpy.full(4, 0.373264085082e-2),
                numpy.full(4, 0.108083624686e-1),
                numpy.full(4, 0.206089282780e-1),
                numpy.full(4, 0.280439642185e-1),
                numpy.full(4, 0.319642570619e-1),
                numpy.full(4, 0.132307491853e-13),
                numpy.full(4, 0.358804444921e-6),
                numpy.full(4, 0.377203097373e-2),
                numpy.full(4, 0.694302275949e-2),
                numpy.full(4, 0.651880370677e-2),
                numpy.full(4, 0.151222484324e-1),
                numpy.full(4, 0.194753613784e-1),
                numpy.full(2, 0.208441813388e-1),
                numpy.full(4, -0.155019627860e-14),
                numpy.full(4, -0.339816478433e-7),
                numpy.full(4, 0.102363309226e-2),
                numpy.full(4, 0.296184898206e-2),
                numpy.full(4, 0.463894901111e-2),
                numpy.full(4, 0.350689331377e-2),
                numpy.full(4, 0.655751993234e-2),
                numpy.full(4, 0.846360999901e-2),
                numpy.full(4, 0.203393765743e-15),
                numpy.full(4, 0.508117261709e-8),
                ])
            self.points = numpy.concatenate([
                _s4(0.950125098376e-1, 0.968160239507),
                _s4(0.950125098376e-1, 0.836031107326),
                _s4(0.950125098376e-1, 0.613371432700),
                _s4(0.950125098376e-1, 0.324253423403),
                _s20(0.950125098376e-1),
                _s4(0.281603550779, 0.950894801068),
                _s4(0.281603550779, 0.771492557752),
                _s4(0.281603550779, 0.500482824502),
                _s4(0.281603550779, 0.173233147858),
                _s4(0.281603550779, 0.114815304093e+1),
                _s4(0.458016777657, 0.976231599013),
                _s4(0.458016777657, 0.860064479607),
                _s4(0.458016777657, 0.639857912633),
                _s4(0.458016777657, 0.340595575643),
                _s20(0.458016777657),
                _s4(0.458016777657, 0.166045668443e+1),
                _s4(0.617876244403, 0.101356615914e+1),
                _s4(0.617876244403, 0.936005992420),
                _s4(0.617876244403, 0.760621522239),
                _s4(0.617876244403, 0.496236604823),
                _s4(0.617876244403, 0.172327489281),
                _s4(0.617876244403, 0.217081790217e+1),
                _s4(0.755404408355, 0.974824439241),
                _s4(0.755404408355, 0.854811906217),
                _s4(0.755404408355, 0.634285445805),
                _s4(0.755404408355, 0.337501230340),
                _s20(0.755404408355),
                _s4(0.755404408355, 0.116365304999e+1),
                _s4(0.755404408355, 0.262105563727e+1),
                _s4(0.865631202388, 0.984091441674),
                _s4(0.865631202388, 0.912693883002),
                _s4(0.865631202388, 0.746762989660),
                _s4(0.865631202388, 0.488326550840),
                _s4(0.865631202388, 0.169703337599),
                _s4(0.865631202388, 0.298557375978e+1),
                _s4(0.865631202388, 0.123267231547e+1),
                _s4(0.944575023073, 0.976647356726),
                _s4(0.944575023073, 0.882226287758),
                _s4(0.944575023073, 0.786493198166),
                _s4(0.944575023073, 0.610059075287),
                _s4(0.944575023073, 0.327509702777),
                _s20(0.944575023073),
                _s4(0.944575023073, 0.324790845348e+1),
                _s4(0.944575023073, 0.131269184417e+1),
                _s4(0.989400934992, 0.985756803095),
                _s4(0.989400934992, 0.914705094008),
                _s4(0.989400934992, 0.769393763416),
                _s4(0.989400934992, 0.616728021409),
                _s4(0.989400934992, 0.442872401546),
                _s4(0.989400934992, 0.158216672941),
                _s4(0.989400934992, 0.339722181604e+1),
                _s4(0.989400934992, 0.135173206143e+1),
                ])
        return


def _s20(a):
    return numpy.array([
        [+a, 0.0],
        [-a, 0.0],
        ])


def _s4(a, b):
    return numpy.array([
        [+a, +b],
        [+a, -b],
        [-a, +b],
        [-a, -b],
        ])
