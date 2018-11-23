# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from ..line_segment import GaussLegendre


class RathodNagarajaVenkatesudu(object):
    """
    H.T. Rathod, K.V. Nagaraja, B. Venkatesudu,
    Symmetric Gauss Legendre quadrature formulas for composite numerical integration
    over a triangular surface,
    Applied Mathematics and Computation 188 (2007) 865–876
    <https://doi.org/10.1016/j.amc.2006.10.041>.

    Abstract:
    This paper first presents a Gauss Legendre quadrature method for numerical
    integration of , where f(x, y) is an analytic function in x, y and T is the standard
    triangular surface: {(x, y)∣0 ⩽ x, y ⩽ 1, x + y ⩽ 1} in the Cartesian two
    dimensional (x, y) space. We then use a transformation x = x(ξ, η), y = y(ξ, η) to
    change the integral I to an equivalent integral , where S is now the 2-square in (ξ,
    η) space: {(ξ, η)∣ − 1 ⩽ ξ, η ⩽ 1}. We then apply the one dimensional Gauss Legendre
    quadrature rules in ξ and η variables to arrive at an efficient quadrature rule with
    new weight coefficients and new sampling points. We then propose the discretisation
    of the standard triangular surface T into n2 right isosceles triangular surfaces Ti
    (i = 1(1)n2) each of which has an area equal to 1/(2n2) units. We have again shown
    that the use of affine transformation over each Ti and the use of linearity property
    of integrals lead to the result: where  and x = xi(X, Y) and y = yi(X, Y) refer to
    affine transformations which map each Ti in (x, y) space into a standard triangular
    surface T in (X, Y) space. We can now apply Gauss Legendre quadrature formulas which
    are derived earlier for I to evaluate the integral . We observe that the above
    procedure which clearly amounts to Composite Numerical Integration over T and it
    converges to the exact value of the integral , for sufficiently large value of n,
    even for the lower order Gauss Legendre quadrature rules. We have demonstrated this
    aspect by applying the above explained Composite Numerical Integration method to
    some typical integrals.
    """

    def __init__(self, index, symbolic=False):
        self.name = "RathodNagarajaVenkatesudu({})".format(index)

        gl = GaussLegendre(index)

        p1 = numpy.multiply.outer(1 + gl.points, 1 + gl.points) / 4
        p2 = numpy.multiply.outer(1 + gl.points, 1 - gl.points) / 4
        w = numpy.multiply.outer((1 + gl.points) * gl.weights, gl.weights) / 8

        self.points = numpy.array([p1.flatten(), p2.flatten()]).T
        self.weights = 2 * w.flatten()

        self.degree = 2 * (index - 1)
        return
