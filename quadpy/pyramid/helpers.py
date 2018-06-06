# -*- coding: utf-8 -*-
#


def _s4(a, z):
    return [[+a, +a, z], [-a, +a, z], [+a, -a, z], [-a, -a, z]]


def _s4_0(a, z):
    return [[+a, 0.0, z], [-a, 0.0, z], [0.0, +a, z], [0.0, -a, z]]
