# translate to barycentric coordinates
import numpy

# reference triangle:
x0 = numpy.array([-1.0, -1.0/numpy.sqrt(3.0)])
x1 = numpy.array([+1.0, -1.0/numpy.sqrt(3.0)])
x2 = numpy.array([0.0, 2.0/numpy.sqrt(3.0)])

T = numpy.array([
    [x1[0] - x0[0], x2[0] - x0[0]],
    [x1[1] - x0[1], x2[1] - x0[1]],
    ])

S = numpy.array([
    [0.0000000000000000E+00, .80383378476840441740523498549521],
    [0.0000000000000000E+00, -0.47391934523147540911429282520838],
    [0.0000000000000000E+00, -0.0],
    ])
# s = [-0.4104261923153453E+00. 0.0000000000000000E+00]
# s = [0.6961404780296310E+00. 0.0000000000000000E+00]

for s in S:
    sol = numpy.linalg.solve(T, numpy.array(s) - x0)
    lmbda = numpy.array([sol[0], sol[1], 1.0-sol[0]-sol[1]])
    print('%.15e %.15e %.15e' % (lmbda[0], lmbda[1], lmbda[2]))


print
W = numpy.array([
    0.82872641363789568276918148034109e-01,
    0.87120251975907374578897626781336e-01,
    0.49352775485718467280720708817387e-01,
    ])
M = [3, 3, 1]
for w, m in zip(W, M):
    print('%.15e' % (w / 0.21934566882541541013653648363283 / m))
