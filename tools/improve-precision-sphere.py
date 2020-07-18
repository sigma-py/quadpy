import json

import numpy
import orthopy
from scipy.optimize import minimize

from quadpy.u3._helpers import expand_symmetries_points_only


def improve_precision_sphere():
    filepath = "/home/nschloe/rcs/quadpy/quadpy/u3/_heo_xu/heo_xu_21a.json"
    with open(filepath, "r") as fh:
        content = json.load(fh)

    degree = content["degree"]
    # name = content["name"]
    # tol = content["test_tolerance"]
    keys = list(content["data"].keys())
    num_symm = [len(item[0]) for item in content["data"].values()]

    def f(x):
        # convert x to dictionary (without weights)
        x_split = numpy.split(x, splits)
        vals = [item.reshape(shape) for item, shape in zip(x_split, shapes)]
        d = dict(zip(keys, vals))
        points, len_symm = expand_symmetries_points_only(d)

        assert not numpy.any(numpy.isnan(points))

        # evaluate all orthogonal polynomials up to `degree` at all points
        evaluator = orthopy.u3.EvalCartesian(points, scaling="quantum mechanic")
        A2 = numpy.concatenate([next(evaluator) for _ in range(degree + 1)])

        assert sum(a * b for a, b in zip(len_symm, num_symm)) == A2.shape[1]

        # sum up all columns which belong to the same weight
        k = 0
        sums = []
        for lsym, nsym in zip(len_symm, num_symm):
            for i in range(nsym):
                # idx = numpy.arange(170)
                # ii  = idx[k + i : k + lsym * nsym : nsym]
                # print("xx", ii, len(ii))
                sums.append(
                    numpy.sum(A2[:, k + i : k + lsym * nsym : nsym], axis=1)
                )
            k += lsym * nsym
        A = numpy.column_stack(sums)

        # The exact values are 0 except for the first entry
        b = numpy.zeros(A.shape[0])
        b[0] = 1.0 / numpy.sqrt(4 * numpy.pi)

        w, res, rank, s = numpy.linalg.lstsq(A, b, rcond=None)
        assert numpy.all(numpy.abs(w.imag) < 1.0e-15)
        # assert res[0] < 1.0e-12
        w = w.real
        return numpy.sqrt(res[0])

    keys = list(content["data"].keys())
    values = list(content["data"].values())

    values_without_weights = [val[1:] for val in values]
    shapes = [numpy.array(val).shape for val in values_without_weights]
    sizes = [numpy.array(val).size for val in values_without_weights]
    splits = numpy.cumsum(sizes)[:-1]

    # collect and concatenate all point coords
    x0 = numpy.concatenate([numpy.array(val).flat for val in values_without_weights])

    out = minimize(f, x0, method="Nelder-Mead", tol=1.0e-17)
    print(out.status, out.nfev)
    print(out.message)
    assert out.success
    print()
    for x in out.x:
        print(f"{x:.15e}")


if __name__ == "__main__":
    improve_precision_sphere()
