# -*- coding: utf-8 -*-
#
import numpy


class CoolsHaegemans(object):
    '''
    R. Cools, A. Haegemans,
    Construction of minimal cubature formulae for the square and the triangle
    using invariant theory,
    Department of Computer Science, K.U.Leuven,
    TW Reports vol:TW96, Sept. 1987,
    <https://lirias.kuleuven.be/handle/123456789/131869>.
    '''
    def __init__(self, index):
        self.name = 'CH(%d)' % index
        assert index == 1
        self.weights = 2.0 * numpy.concatenate([
            numpy.full(3, 0.16058343856681218798E-09),
            numpy.full(3, 0.26530624434780379347E-01),
            numpy.full(3, 0.29285717640155892159E-01),
            numpy.full(3, 0.43909556791220782402E-01),
            numpy.full(3, 0.66940767639916174192E-01),
            ])
        self.bary = numpy.concatenate([
            _r3(
                0.34579201116826902882E+00,
                0.36231682215692616667E+01
                ),
            _r3(
                0.65101993458939166328E-01,
                0.87016510156356306078E+00
                ),
            _r3(
                0.65177530364879570754E+00,
                0.31347788752373300717E+00
                ),
            _r3(
                0.31325121067172530696E+00,
                0.63062143431895614010E+00
                ),
            _r3(
                0.51334692063945414949E+00,
                0.28104124731511039057E+00
                ),
            ])
        self.degree = 8
        # elif index == 2:
        #     self.weights = 2.0 * numpy.concatenate([
        #         numpy.full(3, 0.15319130036758557631E-06),
        #         numpy.full(3, 0.13260526227928785221E-01),
        #         numpy.full(3, 0.15646439344539042136E-01),
        #         numpy.full(3, 0.21704258224807323311E-01),
        #         numpy.full(3, 0.21797613600129922367E-01),
        #         numpy.full(3, 0.38587913508193459468E-01),
        #         numpy.full(3, 0.39699584282594413022E-01),
        #         numpy.full(1, 0.47910534861520060665E-01),
        #         ])
        #     self.bary = numpy.concatenate([
        #         _r3(
        #             +0.58469201683584513031E-01,
        #             -0.54887778772527519316E+00
        #             ),
        #         _r3(
        #             0.50849285064031410705E-01,
        #             0.90799059794957813439E+00
        #             ),
        #         _r3(
        #             0.51586732419949574487E+00,
        #             0.46312452842927062902E+00
        #             ),
        #         _r3(
        #             0.24311033191739048230E+00,
        #             0.72180595182371959467E-00
        #             ),
        #         _r3(
        #             0.75397765920922660134E-00,
        #             0.20647569839132397633E+00
        #             ),
        #         _r3(
        #             0.42209207910846960294E-00,
        #             0.12689533413411127327E+00
        #             ),
        #         _r3(
        #             0.19823878346663354068E+00,
        #             0.62124412566393319745E+00
        #             ),
        #         numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]])
        #         ])
        #     self.degree = 10

        self.points = self.bary[:, 1:]
        return


def _r3(a, b):
    # rotation group R_3
    c = 1.0 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])
