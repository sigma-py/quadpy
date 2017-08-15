# -*- coding: utf-8 -*-
#
from math import pi

from ..nball import integrate
from .. import helpers


def show(
        scheme,
        backend='mpl'
        ):
    '''Displays scheme for 3D ball quadrature.
    '''
    helpers.backend_to_function[backend](
            scheme.points,
            scheme.weights,
            volume=4.0/3.0*pi,
            edges=[]
            )
    return
