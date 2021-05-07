import sys


def optimize(content):
    if "domain" not in content:
        print('Missing key "domain".', file=sys.stderr)
        exit(1)

    domain = content["domain"]
    if domain.lower() == "t2":
        return _optimize_t2(content)
    elif domain.lower() == "s2":
        return _optimize_s2(content)
    elif domain.lower() == "c2":
        return _optimize_c2(content)
    elif domain.lower() == "u3":
        return _optimize_u3(content)

    print(f'Don\'t know how to optimize domain "{domain}".', file=sys.stderr)
    exit(1)


def _optimize_u3(content):
    import orthopy

    from .u3._helpers import (
        _scheme_from_dict,
        expand_symmetries,
        expand_symmetries_points_only,
    )

    return _optimize(
        content,
        expand_symmetries,
        expand_symmetries_points_only,
        _scheme_from_dict,
        get_evaluator=lambda points: orthopy.u3.EvalCartesian(
            points, scaling="quantum mechanic"
        ),
    )


def _optimize_c2(content):
    import orthopy

    from .c2._helpers import (
        _scheme_from_dict,
        expand_symmetries,
        expand_symmetries_points_only,
    )

    return _optimize(
        content,
        expand_symmetries,
        expand_symmetries_points_only,
        _scheme_from_dict,
        get_evaluator=lambda points: orthopy.cn.Eval(points),
    )


def _optimize_s2(content):
    import orthopy

    from .s2._helpers import (
        _scheme_from_dict,
        expand_symmetries,
        expand_symmetries_points_only,
    )

    return _optimize(
        content,
        expand_symmetries,
        expand_symmetries_points_only,
        _scheme_from_dict,
        get_evaluator=lambda points: orthopy.s2.zernike.Eval(points, scaling="normal"),
    )


def _optimize_t2(content):
    import orthopy

    from .helpers import expand_symmetries, expand_symmetries_points_only
    from .t2._helpers import _scheme_from_dict

    # return _optimize_weights_as_variables(
    return _optimize(
        content,
        expand_symmetries,
        expand_symmetries_points_only,
        _scheme_from_dict,
        get_evaluator=lambda points: orthopy.t2.Eval(points, scaling="normal"),
    )


def _optimize(
    content,
    expand_symmetries,
    expand_symmetries_points_only,
    scheme_from_dict,
    get_evaluator,
):
    """Compute the weights from the points via a least-squares problem. Only the point
    coordinates are variables.
    """
    import numpy as np
    from scipy.optimize import minimize

    degree = content["degree"]
    keys = list(content["data"].keys())
    num_symm = [len(item[0]) for item in content["data"].values()]

    def x_to_dict(x):
        # convert x to dictionary (without weights)
        x_split = np.split(x, splits)
        vals = [item.reshape(shape) for item, shape in zip(x_split, shapes)]
        d = dict(zip(keys, vals))
        return d

    def get_w_from_x(x):
        d = x_to_dict(x)
        points, len_symm = expand_symmetries_points_only(d, dim=2)

        if np.any(np.isnan(points)):
            # return some "large" residual value
            return None, None, None, 1.0

        # evaluate all orthogonal polynomials up to `degree` at all points
        evaluator = get_evaluator(points)
        A2 = np.concatenate([next(evaluator) for _ in range(degree + 1)])

        assert sum(a * b for a, b in zip(len_symm, num_symm)) == A2.shape[1]

        # sum up all columns which belong to the same weight
        k = 0
        sums = []
        for lsym, nsym in zip(len_symm, num_symm):
            for i in range(nsym):
                sums.append(np.sum(A2[:, k + i : k + lsym * nsym : nsym], axis=1))
            k += lsym * nsym

        A = np.column_stack(sums)

        # The exact values are 0 except for the first entry
        b = np.zeros(A.shape[0])
        b[0] = evaluator.int_p0
        # b[0] /= np.pi  # necessary for S2

        w, res, rank, s = np.linalg.lstsq(A, b, rcond=None)

        # spherical harmonics (for u3) are complex-valued
        assert np.all(np.abs(w.imag) < 1.0e-15)
        w = w.real

        if rank < min(A.shape):
            # return some "large" residual value
            return None, None, None, 1.0

        return A, b, w, np.sqrt(res[0])

    def f(x):
        _, _, _, res = get_w_from_x(x)
        # print(f"f(x) = {res:.6e}")
        return res

    keys = list(content["data"].keys())
    values = list(content["data"].values())

    values_without_weights = [val[1:] for val in values]
    shapes = [np.array(val).shape for val in values_without_weights]
    sizes = [np.array(val).size for val in values_without_weights]
    splits = np.cumsum(sizes)[:-1]

    # collect and concatenate all point coords
    x0 = np.concatenate([np.array(val).flat for val in values_without_weights])

    # TODO check initial residual with original weights
    r0 = f(x0)

    out = minimize(f, x0, method="Nelder-Mead", tol=1.0e-17)

    # Don't fail on `not out.success`. It could be because of
    # ```
    # Maximum number of function evaluations has been exceeded
    # ```
    # but the scheme could still have improved.

    # out.fun == f(out.x) exactly
    if r0 <= out.fun:
        raise RuntimeError("Optimization couldn't improve scheme.")

    # compute max(err)
    A, b, w, _ = get_w_from_x(out.x)
    # assert A is not None
    # max_res = np.max(np.abs(A @ w - b))

    # Compute max_res exactly like in the tests
    d = x_to_dict(out.x)
    # prepend weights
    k = 0
    for key, value in d.items():
        if len(value.shape) == 1:
            n = 1
            d[key] = [[w[k]]]
        else:
            n = value.shape[1]
            d[key] = np.column_stack([w[k : k + n], value.T]).T
        k += n
    # content["data"] = d
    # scheme = scheme_from_dict(content)
    # max_res = max(scheme.compute_residuals(degree))

    # print(max_res)

    return d, out.fun, np.linalg.cond(A)


def _optimize_weights_as_variables(
    content,
    expand_symmetries,
    expand_symmetries_points_only,
    scheme_from_dict,
    get_evaluator,
):
    """Treat the weights as variables.

    It turns out from numerical experiments that this is inferior to the above approach.
    TODO Information about the derivative + using Newton could improve it a lot.
    """
    import numpy as np
    from scipy.optimize import minimize

    degree = content["degree"]
    keys = list(content["data"].keys())
    A = None

    def x_to_dict(x):
        # convert x to dictionary (without weights)
        x_split = np.split(x, splits)
        vals = [item.reshape(shape) for item, shape in zip(x_split, shapes)]
        d = dict(zip(keys, vals))
        return d

    def get_system(x):
        d = x_to_dict(x)
        points, weights = expand_symmetries(d)

        if np.any(np.isnan(points)):
            # return some "large" residual value
            return 1.0

        # evaluate all orthogonal polynomials up to `degree` at all points
        evaluator = get_evaluator(points.T)
        A = np.concatenate([next(evaluator) for _ in range(degree + 1)])

        # The exact values are 0 except for the first entry
        b = np.zeros(A.shape[0])
        b[0] = evaluator.int_p0
        return A, weights, b

    def f(x):
        A, weights, b = get_system(x)
        res = A @ weights - b
        res_norm = np.sqrt(np.dot(res, res))
        return res_norm

    keys = list(content["data"].keys())
    values = list(content["data"].values())

    shapes = [np.array(val).shape for val in values]
    sizes = [np.array(val).size for val in values]
    splits = np.cumsum(sizes)[:-1]

    # collect and concatenate all point coords
    x0 = np.concatenate([np.array(val).flat for val in values])

    # TODO check initial residual with original weights
    r0 = f(x0)

    out = minimize(f, x0, method="Nelder-Mead", tol=1.0e-17)

    # Don't fail on `not out.success`. It could be because of
    # ```
    # Maximum number of function evaluations has been exceeded
    # ```
    # but the scheme could still have improved.

    # out.fun == f(out.x) exactly
    if r0 <= out.fun:
        raise RuntimeError("Optimization couldn't improve scheme.")

    d = x_to_dict(out.x)

    A, _, _ = get_system(out.x)

    return d, out.fun, np.linalg.cond(A)
