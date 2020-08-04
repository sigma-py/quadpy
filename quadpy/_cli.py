import sys

import numpy


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

    from .s2._helpers import expand_symmetries, expand_symmetries_points_only

    return _optimize(
        content,
        expand_symmetries,
        expand_symmetries_points_only,
        get_evaluator=lambda points: orthopy.s2.zernike.Eval(points, scaling="normal"),
    )


def _optimize_t2(content):
    import orthopy

    from .t2._helpers import (
        _scheme_from_dict,
        expand_symmetries,
        expand_symmetries_points_only,
    )

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
    import numpy
    from scipy.optimize import minimize

    degree = content["degree"]
    keys = list(content["data"].keys())
    num_symm = [len(item[0]) for item in content["data"].values()]

    def x_to_dict(x):
        # convert x to dictionary (without weights)
        x_split = numpy.split(x, splits)
        vals = [item.reshape(shape) for item, shape in zip(x_split, shapes)]
        d = dict(zip(keys, vals))
        return d

    def get_w_from_x(x):
        d = x_to_dict(x)
        points, len_symm = expand_symmetries_points_only(d)

        if numpy.any(numpy.isnan(points)):
            # return some "large" residual value
            return None, None, None, 1.0

        # evaluate all orthogonal polynomials up to `degree` at all points
        evaluator = get_evaluator(points)
        A2 = numpy.concatenate([next(evaluator) for _ in range(degree + 1)])

        assert sum(a * b for a, b in zip(len_symm, num_symm)) == A2.shape[1]

        # sum up all columns which belong to the same weight
        k = 0
        sums = []
        for lsym, nsym in zip(len_symm, num_symm):
            for i in range(nsym):
                sums.append(numpy.sum(A2[:, k + i : k + lsym * nsym : nsym], axis=1))
            k += lsym * nsym

        A = numpy.column_stack(sums)

        # The exact values are 0 except for the first entry
        b = numpy.zeros(A.shape[0])
        b[0] = evaluator.int_p0

        w, res, rank, s = numpy.linalg.lstsq(A, b, rcond=None)

        # spherical harmonics (for u3) are complex-valued
        assert numpy.all(numpy.abs(w.imag) < 1.0e-15)
        w = w.real

        if rank < min(A.shape):
            # return some "large" residual value
            return None, None, None, 1.0

        return A, b, w, numpy.sqrt(res[0])

    def f(x):
        _, _, _, res = get_w_from_x(x)
        # print(f"f(x) = {res:.6e}")
        return res

    keys = list(content["data"].keys())
    values = list(content["data"].values())

    values_without_weights = [val[1:] for val in values]
    shapes = [numpy.array(val).shape for val in values_without_weights]
    sizes = [numpy.array(val).size for val in values_without_weights]
    splits = numpy.cumsum(sizes)[:-1]

    # collect and concatenate all point coords
    x0 = numpy.concatenate([numpy.array(val).flat for val in values_without_weights])

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
    # max_res = numpy.max(numpy.abs(A @ w - b))

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
            d[key] = numpy.column_stack([w[k : k + n], value.T]).T
        k += n
    # content["data"] = d
    # scheme = scheme_from_dict(content)
    # max_res = max(scheme.compute_residuals(degree))

    # print(max_res)

    return d, out.fun, numpy.linalg.cond(A)


def main():
    import json

    import fjson

    args = _get_parser().parse_args()

    with open(args.infile, "r") as f:
        content = json.load(f)

    name = content["name"]
    try:
        new_data, max_err, cond_in_solution = optimize(content)
    except RuntimeError:
        print(f"{name}: Could not improve scheme any further.")
    else:
        if "weight factor" in content:
            w = content["weight factor"]
            for key, item in new_data.items():
                new_data[key][0] = (numpy.array(item)[0] / w).tolist()

        if "comments" in content:
            comments = content["comments"]
        else:
            comments = []

        text = "precision improved with quadpy-optimize"
        if text not in comments:
            content["comments"] = comments + [text]

        prev_tol = content["test_tolerance"]
        content["data"] = new_data
        content["test_tolerance"] = max_err

        # make sure that "data" is written last
        keys = list(content.keys())
        if keys[-1] != "data":
            keys.remove("data")
            keys.append("data")
            content = {key: content[key] for key in keys}

        with open(args.infile, "w") as f:
            # Write out 16 digits such that the numbers are preserved exactly when
            # reading the file back it. This makes sure the computed residuals don't
            # change, not even when they are in the range of machine precision.
            fjson.dump(content, f, indent=2, float_format=".16e")
            # for POSIX compliance:
            f.write("\n")
        print(f"{name}:")
        print(f"Improved max error from   {prev_tol}   to   {max_err}.")
        print(f"condition: {cond_in_solution}")


def _get_parser():
    import argparse
    import sys

    from .__about__ import __version__

    parser = argparse.ArgumentParser(
        description=("Optimize quadrature data."),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("infile", type=str, help="quadrature data file to optimize")

    version_text = "\n".join(
        [
            "quadpy {} [Python {}.{}.{}]".format(
                __version__,
                sys.version_info.major,
                sys.version_info.minor,
                sys.version_info.micro,
            ),
            "Copyright (c) 2016-2020 Nico Schlömer",
        ]
    )

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=version_text,
        help="display version information",
    )
    return parser
