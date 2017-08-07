# -*- coding: utf-8 -*-
#
import itertools
import numpy


def z(n):
    return numpy.zeros((1, n))


def fsd(n, r, d):
    '''Get all permutations of [+-r, +-r, 0, ..., 0] of length n, where +-r
    occurs d times.
    len(out) == 2**d * (n over d).
    n==1:  2*n
    n==2:  2*n*(n-1)
    n==3:  4*n*(n-1)*(n-2) / 3
    '''
    assert 0 <= d <= n
    return combine([[+r, -r]] * d + [[0.0]] * (n-d))


def fsd2(n, r, s, i, j):
    '''Get all permutations of [+-r, +-r, +-s, +-s, 0, ..., 0] of length n,
    with i times the number r and and j times the number s.
    '''
    assert i+j <= n
    return combine([[+r, -r]] * i + [[+s, -s]] * j + [[0.0]] * (n-i-j))


def combine(pools):
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


def pm(n, a):
    '''Return all combinations of [+a, -a] with length n (with repetition).
    len(out) == 2**n.
    '''
    return numpy.array(list(itertools.product([+a, -a], repeat=n)))


def pm_array(v):
    '''Given an array `v = [v0, v1, ..., vn]`, this methods returns all
    combinations of [+-v0, +-v1, ..., +-vn].
    '''
    n = len(v)
    pm_one = numpy.array(list(itertools.product(*(n*[[+1, -1]]))))
    return pm_one * v


def partition(balls, boxes):
    '''Create all nonnegative tuples of length d which sum up to n.
    '''
    # <https://stackoverflow.com/a/36748940/353337>
    # See <https://stackoverflow.com/a/45348441/353337> for an alterantive
    # solution.
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                for x in rec(boxes - 1, i, parent + (balls - i,)):
                    yield x
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))
