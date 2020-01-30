import numpy

import modepy

from . import import_helpers


def _import_tri():
    vr = modepy.quadrature.vioreanu_rokhlin

    for k in range(20):
        q = vr.VioreanuRokhlinSimplexQuadrature(k, 2)
        # The reference triangle is simply
        #     [+1.0, -1.0],
        #     [-1.0, +1.0],
        #     [-1.0, -1.0],
        #
        # Sort the indices to make it easier to identify symmetry groups.
        idx = numpy.argsort(q.weights)
        sorted_weights = q.weights[idx]

        bary = (q.nodes + 1.0) / 2.0
        bary = numpy.vstack([bary, 1.0 - numpy.sum(bary, axis=0)]).T
        sorted_points = bary[idx]

        print(f"elif index == {k}:")
        print("    data = [")
        # identify groups of equal weights and put them out as numpy.full(x, y)
        tol = 1.0e-12
        count = 0

        kk = 0
        last_value = sorted_weights[0]
        for w in sorted_weights:
            if abs(last_value - w) < tol:
                count += 1
            else:
                pts = sorted_points[kk : kk + count]
                kk += count
                print(
                    8 * " "
                    + "(%.15e, %s),"
                    % (last_value, import_helpers.get_symmetry_code_tri(pts))
                )
                last_value = w
                count = 1

        pts = sorted_points[kk : kk + count]
        print(
            8 * " "
            + "({:.15e}, {}),".format(last_value, import_helpers.get_symmetry_code_tri(pts))
        )
        print(8 * " " + "]")

    return


def _import_tet():
    vr = modepy.quadrature.vioreanu_rokhlin

    for k in range(10):
        q = vr.VioreanuRokhlinSimplexQuadrature(k, 3)
        # Sort the indices to make it easier to identify symmetry groups.
        idx = numpy.argsort(q.weights)
        sorted_weights = q.weights[idx]

        bary = (q.nodes + 1.0) / 2.0
        bary = numpy.vstack([bary, 1.0 - numpy.sum(bary, axis=0)]).T
        sorted_points = bary[idx]

        print(f"elif index == {k}:")
        print("    data = [")
        # identify groups of equal weights
        tol = 1.0e-12
        count = 0
        last_value = sorted_weights[0]

        kk = 0
        for w in sorted_weights:
            if abs(last_value - w) < tol:
                count += 1
            else:
                pts = sorted_points[kk : kk + count]
                kk += count
                print(
                    8 * " "
                    + "(%.15e, %s),"
                    % (last_value, import_helpers.get_symmetry_code_tet(pts))
                )
                last_value = w
                count = 1

        pts = sorted_points[kk : kk + count]
        print(
            8 * " "
            + "({:.15e}, {}),".format(last_value, import_helpers.get_symmetry_code_tet(pts))
        )
        print(8 * " " + "]")

    return


if __name__ == "__main__":
    _import_tri()
    # _import_tet()
