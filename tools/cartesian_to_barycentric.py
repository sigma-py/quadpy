'''
Generate code from cartesian coordinates.
'''
import numpy

# reference triangle:
x0 = numpy.array([-1.0, -1.0/numpy.sqrt(3.0)])
x1 = numpy.array([+1.0, -1.0/numpy.sqrt(3.0)])
x2 = numpy.array([0.0, 2.0/numpy.sqrt(3.0)])

T = numpy.array([
    [x1[0] - x0[0], x2[0] - x0[0]],
    [x1[1] - x0[1], x2[1] - x0[1]],
    ])


x = numpy.array([
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    -0.69780686931593427582366555189730e+00,
    0.00000000000000000000000000000000e+00,
    -0.30903100009613455447142535490429e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00
    ])
y = numpy.array([
    -0.56063064349133316993278274030104e+00,
    0.10883996591237330573849553975247e+01,
    -0.51720719429149531106975531479846e+00,
    0.00000000000000000000000000000000e+00,
    -0.51225507594767383123122215407598e+00,
    0.51562570796758001502147445483555e+00,
    -0.32874839651011214404597359718045e+00
    ])
w = numpy.array([
    0.64438869372269120316756176676570e-02,
    0.42018026730627469477523487256964e-02,
    0.38116505989607355372802824839739e-01,
    0.18340560543312397707554978722912e-01,
    0.50983455788600478135692116331631e-01,
    0.51743930451848519084272818465418e-01,
    0.49515526441757000856785778879781e-01
    ])


# generate barycentric coordinate code
S = numpy.column_stack([x, y])
tol = 1.0e-10
multiplicities = numpy.empty(len(S), dtype=int)
for k, s in enumerate(S):
    sol = numpy.linalg.solve(T, numpy.array(s) - x0)
    lmbda = numpy.array([sol[0], sol[1], 1.0-sol[0]-sol[1]])
    assert abs(sum(lmbda) - 1.0) < tol
    # print('%.15e %.15e %.15e' % (lmbda[0], lmbda[1], lmbda[2]))
    diffs = numpy.array([
        abs(lmbda[0] - lmbda[1]),
        abs(lmbda[1] - lmbda[2]),
        abs(lmbda[2] - lmbda[0]),
        ])
    if all(diffs < tol):
        print('_s3(),')
        multiplicities[k] = 1
    elif diffs[0] < tol:
        print('_s21(%.15e),' % lmbda[0])
        multiplicities[k] = 3
    elif diffs[1] < tol:
        print('_s21(%.15e),' % lmbda[1])
        multiplicities[k] = 3
    elif diffs[2] < tol:
        print('_s21(%.15e),' % lmbda[2])
        multiplicities[k] = 3
    else:
        print('_s111(%.15e, %.15e)' % (lmbda[0], lmbda[1]))
        multiplicities[k] = 6

# generate weight code
print
for weight, m in zip(w, multiplicities):
    # print('%.15e' % (w / 0.21934566882541541013653648363283 / m))
    print(
        '%.15e * numpy.ones(%d),'
        % (weight / 0.21934566882541541013653648363283 / m, m)
        )
