# translate to barycentric coordinates
import numpy

# triangle:
x0 = numpy.array([-0.5, -numpy.sqrt(3.0)/2.0])
x1 = numpy.array([-0.5, +numpy.sqrt(3.0)/2.0])
x2 = numpy.array([1.0, 0.0])

T = numpy.array([
    [x1[0] - x0[0], x2[0] - x0[0]],
    [x1[1] - x0[1], x2[1] - x0[1]],
    ])

# s = [0.0000000000000000E+00. 0.0000000000000000E+00]
# s = [-0.4104261923153453E+00. 0.0000000000000000E+00]
# s = [0.6961404780296310E+00. 0.0000000000000000E+00]

for s in S:
    sol = numpy.linalg.solve(T, numpy.array(s) - x0)
    lmbda = numpy.array([sol[0], sol[1], 1.0-sol[0]-sol[1]])
    print('%.15e %.15e %.15e' % (lmbda[0], lmbda[1], lmbda[2]))
