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
    return combine(((+r, -r), d), ((0.0,), n-d))


def fsd2(n, r, s, i, j):
    '''Get all permutations of [+-r, +-r, +-s, +-s, 0, ..., 0] of length n,
    with i times the number r and and j times the number s.
    '''
    assert i+j <= n
    return combine(((+r, -r), i), ((+s, -s), j), ((0.0,), n-i-j))


def combine(*elems):
    '''Given an input array with lists of options, e.g.,

    ((+r, -r), 2), ((+s, -s), 2), ((0.0,), 1)

    this methods returns all combinations, i.e.,

    (+r, 0.0, -s, +r, +s),
    (0.0, -s, -r, +r, +s),
    ...
    '''
    # https://stackoverflow.com/a/45321972/353337

    # Could be replaced by Knuth's "Algorithm L"; see, e.g.,
    # <https://stackoverflow.com/a/4250183/353337>.
    def partitions(*sizes):
        if not sizes or all(s <= 0 for s in sizes):
            yield ()
        for i_size, size in enumerate(sizes):
            if size <= 0:
                continue
            next_sizes = \
                sizes[:i_size] + (sizes[i_size] - 1,) + sizes[i_size + 1:]
            for p in partitions(*next_sizes):
                yield (i_size,) + p

    values, sizes = zip(*elems)
    templates = partitions(*sizes)
    prod = [
        itertools.product(*(values[ti] for ti in t)) for t in templates
        ]
    out = numpy.array(list(itertools.chain.from_iterable(prod)))
    return out


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
