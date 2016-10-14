# -*- coding: utf-8 -*-
#
import numpy


def create_monomial_exponents3(degree):
    '''Returns a list of all monomial exponents of degree :degree:.
    '''
    return [
        numpy.array([degree-i-j, i, j])
        for i in range(degree+1)
        for j in range(degree-i+1)
        ]
