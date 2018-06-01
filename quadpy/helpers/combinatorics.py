# -*- coding: utf-8 -*-
#
import itertools
import numpy


def z(n):
    return numpy.zeros((1, n), dtype=int)


def rd(n, items):
    '''Items is an array of 2-tuples of type (value, count). This method
    returns all all permutations of

    [value1, value1, value2, 0, ..., 0]

    of length n, where value_i occurs count_i times.
    '''
    sum_counts = 0
    for item in items:
        _, count = item
        assert count > 0
        sum_counts += count
    assert 0 <= sum_counts <= n

    elems = [
        ((item[0],), item[1])
        for item in items
        ]
    elems += [
        ((0,), n-sum_counts)
        ]

    return combine(*elems)


def fsd(n, *tuples):
    '''tuples is a list of tuples (value, count). This method returns all
    permutations of [+-r, +-r, +-s, +-s, 0, ..., 0] of length n, with `count`
    times the number `value`.
    '''
    pm_tuples = [((+val, -val), count) for val, count in tuples]
    total_count = sum([item[1] for item in tuples])
    assert total_count <= n
    pm_tuples += [
        ((0,), n - total_count)
        ]
    return combine(*pm_tuples)


def fs_array(v):
    '''Given an array v = [v0, v1, ..., vn], this returns all permutation of v
    with plus-minus, i.e., [+v0, +v1, ..., +vn], [-v1, +v0, ..., +vn], etc.
    '''
    elems = [((+vi, -vi), 1) for vi in v]
    return combine(*elems)


def combine(*elems):
    '''Given an input array with lists of options, e.g.,

    ((+r, -r), 2), ((+s, -s), 2), ((0,), 1)

    this methods returns all combinations, i.e.,

    (+r, 0, -s, +r, +s),
    (0, -s, -r, +r, +s),
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


def pm_array0(n, v, idx):
    '''Like pm_array, but put the plus-minused values in a larger array of
    length n at indices idx with the rest filled up with zeros.
    '''
    pm_v = pm_array(v)
    out = numpy.zeros((n, len(pm_v)), dtype=pm_v.dtype)
    out[idx] = pm_v.T
    return out.T


def pm_roll(n, v):
    '''Returns `2**k * n` number of points of dimension `n` such that

    p[0] = [+-v[0], ..., +-v[k], 0, ..., 0]
    p[1] = [0, +-v[0], ..., +-v[k], 0, ..., 0]
    ...
    p[n-1] = [+-v[1], ..., +-v[k], 0, ..., 0, +-v[0]]

    with all +- configurations.
    '''
    k = len(v)
    assert k <= n

    pm_v = pm_array(v)

    r0 = numpy.zeros((len(pm_v), n), dtype=pm_v.dtype)
    r0[:, :k] = pm_v

    return numpy.concatenate([
        numpy.roll(r0, i, axis=1)
        for i in range(n)
        ])


# TODO remove
def partition(boxes, balls):
    '''Create all nonnegative tuples of length d which sum up to n.
    '''
    # <https://stackoverflow.com/a/36748940/353337>
    # See <https://stackoverflow.com/a/45348441/353337> for an alterative
    # solution.
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                for x in rec(boxes - 1, i, parent + (balls - i,)):
                    yield x
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))


def get_all_exponents(dim, max_degree):
    '''Get all exponent combinations of dimension `dim` and maximum degree
    `max_degree`. This method is actually meant for evaluating all polynomials
    with these exponents.

    This problem is similar to the weak_compositions, e.g.,
    <https://stackoverflow.com/a/36748940/353337>. The solution here, however,
    only ever adds 1, making it better suited for a possible extension with
    actual polynomial evaluation.
    '''
    def augment(exponents):
        '''This function takes the values and exponents of a given monomial
        level, e.g., [(1,0,0), (0,1,0), (0,0,1)], and augments them by one
        level, i.e., [(2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)].
        The methods works for all dimension and is based on the observation
        that the augmentation can happen by

          1. adding 1 to all exponents from the previous level (i.e.,
             multiplication with x[0]),
          2. adding 1 to the second exponent of all exponent tuples where
             ex[0]==0 (i.e., multiplication with x[1]),
          3. adding 1 to the third exponent of all exponent tuples where
             ex[0]==0, ex[1]=0 (i.e., multiplication with x[2]),

        etc. The function call is recursive.
        '''

        if len(exponents) == 0 or len(exponents[0]) == 0:
            return []

        idx_leading_zero = [
            k for k in range(len(exponents)) if exponents[k][0] == 0
            ]
        exponents_with_leading_zero = \
            [exponents[k][1:] for k in idx_leading_zero]
        # val1 = vals[idx_leading_zero]
        # x1 = x[1:]
        out1 = augment(exponents_with_leading_zero)

        # increment leading exponent by 1
        out = [[e[0]+1] + e[1:] for e in exponents]
        # vals0 = vals * x[0]

        out += [[0] + e for e in out1]
        # out_vals = numpy.concatenate([vals0, vals1])
        return out

    # dim = x.shape[0]

    # Initialization, level 0
    exponents = [dim * [0]]
    # vals = numpy.array(numpy.ones(x.shape[1:]))

    # all_vals = []
    all_exponents = []

    # all_vals.append(vals)
    all_exponents.append(exponents)
    for _ in range(max_degree):
        exponents = augment(exponents)
        # all_vals.append(vals)
        all_exponents.append(exponents)

    return all_exponents
