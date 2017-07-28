# -*- coding: utf-8 -*-
#
import itertools
import numpy


def _z(n):
    return numpy.zeros((1, n))


def _fs11(n, r, s):
    '''Get all permutations of [+-r, +-s, ..., +-s] of length n.
    len(out) == n * 2**n.
    '''
    # <https://stackoverflow.com/a/45321972/353337>
    k1 = 1
    k2 = n-1
    idx = itertools.combinations(range(k1 + k2), k1)
    vs = ((s if j not in i else r for j in range(k1 + k2)) for i in idx)
    return numpy.array(list(itertools.chain.from_iterable(
        itertools.product(*((+vij, -vij) for vij in vi)) for vi in vs
        )))


def _fsd(n, r, d):
    '''Get all permutations of [+-r, +-r, 0, ..., 0] of length n, where +-r
    occurs d times.
    len(out) == 2**d * (n over d).
    n==1:  2*n
    n==2:  2*n*(n-1)
    n==3:  4*n*(n-1)*(n-2) / 3
    '''
    assert 0 <= d <= n
    return _combine([[+r, -r]] * d + [[0.0]] * (n-d))


def _combine(pools):
    '''Given an input array with lists of options, e.g.,

    [[a, b], [c], [d]],

    this methods returns all combinations with one element from each
    subset, e.g.,

    [a, c, d], [a, d, c], [c, d, a], ...
    [b, c, d], [b, d, c], [c, d, b], ...
    '''
    # https://stackoverflow.com/a/45322199/353337
    return numpy.array(list(set(itertools.chain.from_iterable([
        itertools.permutations(x) for x in itertools.product(*pools)
        ]))))


def _s(n, a, b):
    '''Get all permutations of [a, b, ..., b] of length n.
    len(out) == n.
    '''
    out = numpy.full((n, n), b)
    numpy.fill_diagonal(out, a)
    return out


def _s2(n, a):
    '''Get all permutations of [a, a, 0, 0, ..., 0] of length n.
    len(out) == (n-1)*n / 2.
    '''
    return _s11(n, a, a)


def _s11(n, a, b):
    '''Get all permutations of [a, b, 0, 0, ..., 0] of length n.
    len(out) == (n-1)*n.
    '''
    s = [a, b] + (n-2) * [0]
    # First permutations, then set can be really inefficient if items are
    # repeated. Check out <https://stackoverflow.com/q/6284396/353337> for
    # improvements.
    return numpy.array(list(set(itertools.permutations(s, n))))


def _pm(n, a):
    '''Return all combinations of [+a, -a] with length n (with repetition).
    len(out) == 2**n.
    '''
    return numpy.array(list(itertools.product([+a, -a], repeat=n)))
