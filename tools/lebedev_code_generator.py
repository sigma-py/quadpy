# -*- coding: utf-8 -*-
#
'''
This little helper takes Lebedev point and weight data from [1] and produces
Python code compatible with this library.

[1] https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
'''
import numpy as np


def read(filename):
    data = np.loadtxt(filename)

    theta = data[:, 1] / 180.0 * np.pi
    phi = data[:, 0] / 180.0 * np.pi
    weights = data[:, 2]

    X = np.c_[
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
        ]
    return X, weights


def chunk_data(weights):
    # divide the weight vector into chunks of 6, 8, 12, 24, or 48
    chunks = []
    k = 0
    ref_weight = 0.0
    while k < len(weights):
        if len(chunks) > 0 and abs(weights[k] - ref_weight) < 1.0e-12:
            chunks[-1].append(k)
        else:
            chunks.append([k])
            ref_weight = weights[k]
        k += 1
    return chunks


def sort_into_symmetry_classes(weights, X):
    data = {
        'a1': [],
        'a2': [],
        'a3': [],
        'pq0': [],
        'llm': [],
        'rSW': [],
        }
    for c in chunks:
        if len(c) == 6:
            data['a1'].append(weights[c[0]])
        elif len(c) == 12:
            data['a2'].append(weights[c[0]])
        elif len(c) == 8:
            data['a3'].append(weights[c[0]])
        elif len(c) == 24 and min(X[c[0]]) < 1.0e-12:
            assert abs(np.linalg.norm(X[c[0]]) - 1.0) < 1.0e-12
            if abs(X[c[0], 0]) < 1.0e-12:
                p = X[c[0], 1]
                q = X[c[0], 2]
            elif abs(X[c[0], 1]) < 1.0e-12:
                p = X[c[0], 2]
                q = X[c[0], 0]
            elif abs(X[c[0], 2]) < 1.0e-12:
                p = X[c[0], 0]
                q = X[c[0], 1]
            else:
                raise ValueError('')
            data['pq0'].append((
                weights[c[0]],
                p, q
                ))
        elif len(c) == 24:
            if abs(X[c[0], 0] - X[c[0], 1]) < 1.0e-12:
                l = X[c[0], 1]
                m = X[c[0], 2]
            elif abs(X[c[0], 1] - X[c[0], 2]) < 1.0e-12:
                l = X[c[0], 2]
                m = X[c[0], 0]
            elif abs(X[c[0], 2] - X[c[0], 0]) < 1.0e-12:
                l = X[c[0], 0]
                m = X[c[0], 1]
            else:
                raise ValueError('')
            data['llm'].append((
                weights[c[0]],
                l, m
                ))
        elif len(c) == 48:
            data['rSW'].append((
                weights[c[0]],
                X[c[0], 0], X[c[0], 1], X[c[0], 2]
                ))
        else:
            raise RuntimeError('')

    return data


def generate_python_code(data):
    # generate the code à la
    # ```
    # self.weights = numpy.concatenate([
    #   0.00051306717973400001 * numpy.ones(6),
    #   64.0 / 2835.0 * numpy.ones(12),
    #   27.0 / 1280.0 * numpy.ones(8),
    #   14641.0 / 725760.0 * numpy.ones(24),
    #   ])
    # self.points = numpy.concatenate([
    #   self.a1(),
    #   self.a2(),
    #   self.a3(),
    #   self.llm(3.0151134457776357367e-01, 9.0453403373329088755e-01)
    #   ])
    # ```
    print('self.weights = numpy.concatenate([')
    for d in data['a1']:
        print('    %0.16e * numpy.ones(6),' % d)
    for d in data['a2']:
        print('    %0.16e * numpy.ones(12),' % d)
    for d in data['a3']:
        print('    %0.16e * numpy.ones(8),' % d)
    for d in data['pq0']:
        print('    %0.16e * numpy.ones(24),' % d[0])
    for d in data['llm']:
        print('    %0.16e * numpy.ones(24),' % d[0])
    for d in data['rSW']:
        print('    %0.16e * numpy.ones(48),' % d[0])
    print('    ])')
    # points
    print('self.points = numpy.concatenate([')
    for d in data['a1']:
        print('    self.a1(),')
    for d in data['a2']:
        print('    self.a2(),')
    for d in data['a3']:
        print('    self.a3(),')
    for d in data['pq0']:
        print('    self.pq0(%0.16e, %0.16e),' % (d[1], d[2]))
    for d in data['llm']:
        print('    self.llm(%0.16e, %0.16e),' % (d[1], d[2]))
    for d in data['rSW']:
        print('    self.rsw(%0.16e, %0.16e, %0.16e),' % (d[1], d[2], d[3]))
    print('    ])')
    return


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate code from Lebedev data.'
        )
    parser.add_argument(
            'filename',
            metavar='FILE',
            type=str,
            help='Lebedev data file'
            )
    args = parser.parse_args()

    X, weights = read(args.filename)
    chunks = chunk_data(weights)
    data = sort_into_symmetry_classes(weights, X)
    generate_python_code(data)
