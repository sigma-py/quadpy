import modepy
import numpy


def _get_symmetry_code_tri(pts):
    if len(pts) == 1:
        return '_s3()'
    elif len(pts) == 3:
        # Symmetry group [[a, a, b], [a, b, a], [b, a, a]].
        # Find the equal value `a`.
        tol = 1.0e-12
        beta = pts[0] - pts[0][0]
        ct = numpy.count_nonzero(abs(beta) < tol)
        assert ct in [1, 2], beta
        val = pts[0][0] if ct == 2 else pts[0][1]
        return '_s21({:.15e})'.format(val)

    # Symmetry group [[a, b, c], [c, a, b], ...].
    assert len(pts) == 6
    # Take the two largest value from a, b, c.
    pt0 = numpy.sort(pts[0])
    return '_s111({:.15e}, {:.15e})'.format(pt0[2], pt0[1])


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

        print('elif index == {}:'.format(k))
        print('    data = [')
        # identify groups of equal weights and put them out as numpy.full(x, y)
        tol = 1.0e-12
        count = 0

        kk = 0
        last_value = sorted_weights[0]
        for w in sorted_weights:
            if abs(last_value - w) < tol:
                count += 1
            else:
                pts = sorted_points[kk:kk+count]
                kk += count
                print(
                    8*' ' + '(%.15e, %s),' %
                    (last_value, _get_symmetry_code_tri(pts))
                    )
                last_value = w
                count = 1

        pts = sorted_points[kk:kk+count]
        print(
            8*' ' + '(%.15e, %s),' %
            (last_value, _get_symmetry_code_tri(pts))
            )
        print(8*' ' + ']')

    return


def _get_symmetry_code_tet(pts):
    if len(pts) == 1:
        return '_s4()'
    elif len(pts) == 4:
        # Symmetry group [[a, a, a, b], [a, a, b, a], ...].
        # Find the equal value `a`.
        tol = 1.0e-12
        beta = pts[0] - pts[0][0]
        ct = numpy.count_nonzero(abs(beta) < tol)
        if ct == 1:
            return '_s31({:.15e})'.format(pts[0][1])
        else:
            assert ct == 3
            return '_s31({:.15e})'.format(pts[0][0])
    elif len(pts) == 6:
        # Symmetry group [[a, a, b, b], [a, b, b, a], ...].
        return '_s22({:.15e})'.format(pts[0][0])
    elif len(pts) == 12:
        # Symmetry group [[a, a, b, c], [a, c, a, b], ...].
        # Find the double.
        pt0 = numpy.sort(pts[0])
        tol = 1.0e-12
        for k in range(4):
            beta = pt0 - pt0[k]
            ct = numpy.count_nonzero(abs(beta) < tol)
            if ct == 2:
                return '_s211({:.15e}, {:.15e})'.format(pt0[k], pt0[k-1])

    assert len(pts) == 24
    pt0 = numpy.sort(pts[0])
    return '_s1111({:.15e}, {:.15e}, {:.15e})'.format(
            pt0[-1], pt0[-2], pt0[-3]
            )


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

        print('elif index == {}:'.format(k))
        print('    data = [')
        # identify groups of equal weights
        tol = 1.0e-12
        count = 0
        last_value = sorted_weights[0]

        kk = 0
        for w in sorted_weights:
            if abs(last_value - w) < tol:
                count += 1
            else:
                pts = sorted_points[k:k+count]
                kk += count
                print(
                    8*' ' + '(%.15e, %s),'
                    % (last_value, _get_symmetry_code_tet(pts))
                    )
                last_value = w
                count = 1

        pts = sorted_points[kk:kk+count]
        print(
            8*' ' + '(%.15e, %s),'
            % (last_value, _get_symmetry_code_tet(pts))
            )
        print(8*' ' + ']')

    return


if __name__ == '__main__':
    _import_tri()
    # _import_tet()
