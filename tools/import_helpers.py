import numpy as np


def get_symmetry_code_tri(pts):
    if len(pts) == 1:
        return "_s3()"
    elif len(pts) == 3:
        # Symmetry group [[a, a, b], [a, b, a], [b, a, a]].
        # Find the equal value `a`.
        tol = 1.0e-12
        beta = pts[0] - pts[0][0]
        ct = np.count_nonzero(abs(beta) < tol)
        assert ct in [1, 2], beta
        val = pts[0][0] if ct == 2 else pts[0][1]
        return f"_s21({val:.15e})"

    # Symmetry group [[a, b, c], [c, a, b], ...].
    assert len(pts) == 6
    # Take the two largest value from a, b, c.
    pt0 = np.sort(pts[0])
    return f"_s111({pt0[2]:.15e}, {pt0[1]:.15e})"


def get_symmetry_code_tet(pts):
    if len(pts) == 1:
        return "_s4()"
    elif len(pts) == 4:
        # Symmetry group [[a, a, a, b], [a, a, b, a], ...].
        # Find the equal value `a`.
        tol = 1.0e-12
        beta = pts[0] - pts[0][0]
        ct = np.count_nonzero(abs(beta) < tol)
        if ct == 1:
            return f"_s31({pts[0][1]:.15e})"
        else:
            assert ct == 3
            return f"_s31({pts[0][0]:.15e})"
    elif len(pts) == 6:
        # Symmetry group [[a, a, b, b], [a, b, b, a], ...].
        return f"_s22({pts[0][0]:.15e})"
    elif len(pts) == 12:
        # Symmetry group [[a, a, b, c], [a, c, a, b], ...].
        # Find the double.
        pt0 = np.sort(pts[0])
        tol = 1.0e-12
        for k in range(4):
            beta = pt0 - pt0[k]
            ct = np.count_nonzero(abs(beta) < tol)
            if ct == 2:
                return f"_s211({pt0[k]:.15e}, {pt0[k - 1]:.15e})"

    print(len(pts))
    assert len(pts) == 24
    pt0 = np.sort(pts[0])
    return f"_s1111({pt0[-1]:.15e}, {pt0[-2]:.15e}, {pt0[-3]:.15e})"
