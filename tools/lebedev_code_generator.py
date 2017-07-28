# -*- coding: utf-8 -*-
#
'''
This little helper takes Lebedev point and weight data from [1] and produces
Python code compatible with this library.

[1]
https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
'''
import numpy


def read(filename):
    data = numpy.loadtxt(filename)
    phi_theta = data[:, :2] / 180.0
    weights = data[:, 2]
    return phi_theta, weights


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


def sort_into_symmetry_classes(weights, phi_theta):
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
        elif len(c) == 24:
            if any(abs(phi_theta[c, 1] - 0.5) < 1.0e-12):
                # theta == pi/2   =>   X == [p, q, 0].
                # Find the smallest positive phi that's paired with `theta ==
                # pi/2`; the symmetry is fully characterized by that phi.
                k = numpy.where(abs(phi_theta[c, 1] - 0.5) < 1.0e-12)[0]
                assert len(k) == 8
                k2 = numpy.where(phi_theta[c, 0][k] > 0.0)[0]
                phi_min = numpy.min(phi_theta[c, 0][k][k2])
                data['pq0'].append((weights[c[0]], phi_min))
            else:
                # X = [l, l, m].
                # In this case, there must by exactly two phi with the value
                # pi/4. Take the value of the smaller corresponding `theta`;
                # all points are characterized by it.
                k = numpy.where(abs(phi_theta[c, 0] - 0.25) < 1.0e-12)[0]
                assert len(k) == 2
                k2 = numpy.where(phi_theta[c, 1][k] > 0.0)[0]
                theta_min = numpy.min(phi_theta[c, 1][k][k2])
                data['llm'].append((weights[c[0]], theta_min))
        else:
            assert len(c) == 48
            # This most general symmetry is characterized by two angles; one
            # could take any two here.
            # To make things easier later on, out of the 6 smallest theta, take
            # the one with the smallest positive phi.
            min_theta = numpy.min(phi_theta[c, 1])
            k = numpy.where(abs(phi_theta[c, 1] - min_theta) < 1.0e-12)[0]
            k2 = numpy.where(phi_theta[c, 0][k] > 0.0)[0]
            min_phi = numpy.min(phi_theta[c, 0][k][k2])
            data['rSW'].append((
                weights[c[0]], min_phi, min_theta
                ))

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
        print('    numpy.full(6, %0.16e),' % d)
    for d in data['a2']:
        print('    numpy.full(12, %0.16e),' % d)
    for d in data['a3']:
        print('    numpy.full(8, %0.16e),' % d)
    for d in data['pq0']:
        print('    numpy.full(24, %0.16e),' % d[0])
    for d in data['llm']:
        print('    numpy.full(24, %0.16e),' % d[0])
    for d in data['rSW']:
        print('    numpy.full(48, %0.16e),' % d[0])
    print('    ])')
    # points
    print('self.phi_theta = numpy.concatenate([')
    for d in data['a1']:
        print('    self.a1(),')
    for d in data['a2']:
        print('    self.a2(),')
    for d in data['a3']:
        print('    self.a3(),')
    for d in data['pq0']:
        print('    self.pq0(%0.16e),' % d[1])
    for d in data['llm']:
        print('    self.llm(%0.16e),' % d[1])
    for d in data['rSW']:
        print('    self.rsw(%0.16e, %0.16e),' % (d[1], d[2]))
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

    phi_theta, weights = read(args.filename)
    chunks = chunk_data(weights)
    data = sort_into_symmetry_classes(weights, phi_theta)
    generate_python_code(data)
