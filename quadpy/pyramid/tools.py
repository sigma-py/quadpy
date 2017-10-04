# -*- coding: utf-8 -*-
#
import numpy

from .felippa import Felippa
from .. import helpers


def integrate(f, pyra, scheme, sumfun=helpers.kahan_sum):
    flt = numpy.vectorize(float)

    xi = flt(scheme.points).T
    x = _transform(xi, pyra)
    det = _get_det_J(pyra, xi)

    return sumfun(flt(scheme.weights) * f(x) * abs(det.T), axis=-1)


def _transform(xi, pyra):
    mo = numpy.multiply.outer
    return (
        + mo(0.125*(1.0-xi[0])*(1.0-xi[1])*(1-xi[2]), pyra[0])
        + mo(0.125*(1.0+xi[0])*(1.0-xi[1])*(1-xi[2]), pyra[1])
        + mo(0.125*(1.0+xi[0])*(1.0+xi[1])*(1-xi[2]), pyra[2])
        + mo(0.125*(1.0-xi[0])*(1.0+xi[1])*(1-xi[2]), pyra[3])
        + mo(0.500*(1.0+xi[2]), pyra[4])
        ).T


def _get_det_J(pyra, xi):
    J0 = (
        - numpy.multiply.outer(0.125*(1.0-xi[1])*(1-xi[2]), pyra[0])
        + numpy.multiply.outer(0.125*(1.0-xi[1])*(1-xi[2]), pyra[1])
        + numpy.multiply.outer(0.125*(1.0+xi[1])*(1-xi[2]), pyra[2])
        - numpy.multiply.outer(0.125*(1.0+xi[1])*(1-xi[2]), pyra[3])
        ).T
    J1 = (
        - numpy.multiply.outer(0.125*(1.0-xi[0])*(1-xi[2]), pyra[0])
        - numpy.multiply.outer(0.125*(1.0+xi[0])*(1-xi[2]), pyra[1])
        + numpy.multiply.outer(0.125*(1.0+xi[0])*(1-xi[2]), pyra[2])
        + numpy.multiply.outer(0.125*(1.0-xi[0])*(1-xi[2]), pyra[3])
        ).T
    J2 = (
        - numpy.multiply.outer(0.125*(1.0-xi[0])*(1.0-xi[1]), pyra[0])
        - numpy.multiply.outer(0.125*(1.0+xi[0])*(1.0-xi[1]), pyra[1])
        - numpy.multiply.outer(0.125*(1.0+xi[0])*(1.0+xi[1]), pyra[2])
        - numpy.multiply.outer(0.125*(1.0-xi[0])*(1.0+xi[1]), pyra[3])
        + numpy.multiply.outer(0.500*numpy.ones(1), pyra[4])
        ).T
    det = (
        + J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2]
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]
        )
    return det.T


def show(
        scheme,
        pyra=numpy.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0.5, 0.5, 1.0],
            ]),
        backend='mpl'
        ):
    edges = numpy.array([
        [pyra[0], pyra[1]],
        [pyra[1], pyra[2]],
        [pyra[2], pyra[3]],
        [pyra[3], pyra[0]],
        #
        [pyra[0], pyra[4]],
        [pyra[1], pyra[4]],
        [pyra[2], pyra[4]],
        [pyra[3], pyra[4]],
        ])
    edges = numpy.moveaxis(edges, 1, 2)

    helpers.backend_to_function[backend](
            _transform(scheme.points.T, pyra).T,
            scheme.weights,
            integrate(lambda x: 1.0, pyra, Felippa(1)),
            edges
            )
    return
