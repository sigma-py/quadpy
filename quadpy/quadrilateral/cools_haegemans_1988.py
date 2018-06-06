# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle


class CoolsHaegemans1988(object):
    """
    R. Cools and Ann Haegemans,
    Another Step Forward in Searching for Cubature Formulae with a Minimal
    Number of Knots for the Square,
    Computing 40, 139- 146 (1988),
    <https://doi.org/10.1007/BF02247942>.

    Abstract:
    The knots and weights of a cubature formula are determined by a system of
    nonlinear equations. The number of equations and unknowns can be reduced by
    imposing some structure on the formula.  We are concerned with the
    construction of cubature formulae which are invariant under rotations.
    Using invariant theory, we obtain a smaller system of algebraically
    independent equations.  This is used to construct cubature formulae for the
    square. One of the results is a 24-point formula of degree 11.
    """

    def __init__(self, index):
        self.name = "CH88(%d)" % index
        if index == 1:
            self.degree = 11
            data = [
                (
                    0.0480207633507238145627631759775806,
                    _s4(
                        0.982639223540855472952491497004009,
                        0.698076104549567564776469806174958,
                    ),
                ),
                (
                    0.0660713291645505956736350808495464,
                    _s4(
                        0.825775835902963937302274585289940,
                        0.939486382816736907206432362169896,
                    ),
                ),
                (
                    0.0973867773586681641961204397995472,
                    _s4(
                        0.188586138718641954600324568182745,
                        0.953539528201532015845004266823976,
                    ),
                ),
                (
                    0.211736349998948600503931661356261,
                    _s4(
                        0.812520548304813100489382581912299,
                        0.315623432915254195985609716402104,
                    ),
                ),
                (
                    0.225626061728863387403158016208490,
                    _s4(
                        0.525320250364547762341631887140024,
                        0.712001913075336306549065895123759,
                    ),
                ),
                (
                    0.351158718398245437660391625808574,
                    _s4(
                        0.0416580719120223682735468045377018,
                        0.424847248848669250615430111511957,
                    ),
                ),
            ]
        else:
            assert index == 2
            self.degree = 13
            data = [
                (
                    0.29991838864499131666e-01,
                    _s4(0.77880971155441942252e+00, 0.98348668243987226379e+00),
                ),
                (
                    0.38174421317083669640e-01,
                    _s4(0.95729769978630736566e+00, 0.85955600564163892859e+00),
                ),
                (
                    0.60424923817749980681e-01,
                    _s4(0.13818345986246535375e+00, 0.95892517028753485754e+00),
                ),
                (
                    0.77492738533105339358e-01,
                    _s4(0.94132722587292523695e+00, 0.39073621612946100068e+00),
                ),
                (
                    0.11884466730059560108e+00,
                    _s4(0.47580862521827590507e+00, 0.85007667369974857597e+00),
                ),
                (
                    0.12976355037000271129e+00,
                    _s4(0.75580535657208143627e+00, 0.64782163718701073204e+00),
                ),
                (
                    0.21334158145718938943e+00,
                    _s4(0.69625007849174941396e+00, 0.70741508996444936217e-01),
                ),
                (
                    0.25687074948196783651e+00,
                    _s4(0.34271655604040678941e+00, 0.40930456169403884330e+00),
                ),
                (0.30038211543122536139e+00, _s1()),
            ]

        self.points, self.weights = untangle(data)
        return


def _s1():
    return numpy.array([[0, 0]])


def _s4(a, b):
    return numpy.array([[+a, +b], [-a, -b], [-b, +a], [+b, -a]])
