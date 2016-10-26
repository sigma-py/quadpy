'''
Generate code from cartesian coordinates.
'''
import mpmath
# set precision
mpmath.mp.dps = 32

# reference triangle:
t0 = [-1, -1/mpmath.sqrt(3)]
t1 = [+1, -1/mpmath.sqrt(3)]
t2 = [0, 2/mpmath.sqrt(3)]


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# read data from file
data = []
current_block = None
next_block = 'x'
with open('symq.txt') as f:
    for line in f:
        line = line \
            .replace(',', '') \
            .replace('&', '') \
            .replace('/', '') \
            .replace('D', 'e') \
            .strip()

        if is_float(line):
            if current_block is None:
                current_block = next_block
                if current_block == 'x':
                    data.append({})
            if current_block not in data[-1]:
                data[-1][current_block] = []
            data[-1][current_block].append(mpmath.mp.mpf(line))
        else:
            if current_block is not None:
                if current_block == 'x':
                    next_block = 'y'
                elif current_block == 'y':
                    next_block = 'w'
                elif current_block == 'w':
                    next_block = 'x'
                current_block = None


for i, scheme_data in enumerate(data):
    X = scheme_data['x']
    Y = scheme_data['y']
    W = scheme_data['w']

    print('elif index == %d:' % (i+1))
    print('    bary = numpy.concatenate([')

    # generate barycentric coordinate code
    XY = [[xx, yy] for xx, yy in zip(X, Y)]
    T = mpmath.matrix([
        [t1[0] - t0[0], t2[0] - t0[0]],
        [t1[1] - t0[1], t2[1] - t0[1]],
        ])
    tol = 1.0e-10
    multiplicities = []
    for k, xy in enumerate(XY):
        b = [xy[0] - t0[0], xy[1] - t0[1]]
        sol = mpmath.lu_solve(T, b)
        lmbda = [sol[0], sol[1], 1.0-sol[0]-sol[1]]
        assert abs(sum(lmbda) - 1.0) < tol
        # print('%.15e %.15e %.15e' % (lmbda[0], lmbda[1], lmbda[2]))
        diffs = [
            abs(lmbda[0] - lmbda[1]),
            abs(lmbda[1] - lmbda[2]),
            abs(lmbda[2] - lmbda[0]),
            ]
        if diffs[0] < tol and diffs[1] < tol and diffs[2] < tol:
            print('        _s3(),')
            multiplicities.append(1)
        elif diffs[0] < tol:
            print('        _s21(%s),' % lmbda[0])
            multiplicities.append(3)
        elif diffs[1] < tol:
            print('        _s21(%s),' % lmbda[1])
            multiplicities.append(3)
        elif diffs[2] < tol:
            print('        _s21(%s),' % lmbda[2])
            multiplicities.append(3)
        else:
            print('        _s111(\n            %s,\n            %s\n            ),' % (lmbda[0], lmbda[1]))
            multiplicities.append(6)
    print('        ])')
    print('    self.weights = numpy.concatenate([')
    # generate weight code
    print
    alpha = mpmath.mp.mpf('0.21934566882541541013653648363283')
    for weight, m in zip(W, multiplicities):
        # print('%.15e' % (w / 0.21934566882541541013653648363283 / m))
        print(
            '        %s * numpy.ones(%d),'
            % (weight / alpha / m, m)
            )
    print('        ])')
