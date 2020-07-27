import numpy
import sys


def optimize(content):
    if "domain" not in content:
        print('Missing key "domain".', file=sys.stderr)
        exit(1)

    domain = content["domain"]
    if domain.lower() == "t2":
        return _optimize_t2(content)
    # elif domain.lower() == "s2":
    #     return _optimize_t2(content)

    print(f'Don\'t know how to optimize domain "{domain}".', file=sys.stderr)
    exit(1)


def _optimize_t2(content):
    import numpy
    import orthopy
    from scipy.optimize import minimize

    from .t2._helpers import expand_symmetries_points_only

    degree = content["degree"]
    keys = list(content["data"].keys())
    num_symm = [len(item[0]) for item in content["data"].values()]

    def x_to_dict(x):
        # convert x to dictionary (without weights)
        x_split = numpy.split(x, splits)
        vals = [item.reshape(shape) for item, shape in zip(x_split, shapes)]
        d = dict(zip(keys, vals))
        return d

    def get_Ab(x):
        d = x_to_dict(x)
        points, len_symm = expand_symmetries_points_only(d)

        if numpy.any(numpy.isnan(points)):
            # return some "large" residual value
            return 1.0

        # evaluate all orthogonal polynomials up to `degree` at all points
        evaluator = orthopy.t2.Eval(points, scaling="normal")
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
        b[0] = numpy.sqrt(2)
        return A, b

    def get_w_from_x(x):
        A, b = get_Ab(x)
        w, res, rank, s = numpy.linalg.lstsq(A, b, rcond=None)
        return A, b, w, numpy.sqrt(res[0])

    def f(x):
        _, _, _, res = get_w_from_x(x)
        return res

    keys = list(content["data"].keys())
    values = list(content["data"].values())

    values_without_weights = [val[1:] for val in values]
    shapes = [numpy.array(val).shape for val in values_without_weights]
    sizes = [numpy.array(val).size for val in values_without_weights]
    splits = numpy.cumsum(sizes)[:-1]

    # collect and concatenate all point coords
    x0 = numpy.concatenate([numpy.array(val).flat for val in values_without_weights])

    out = minimize(f, x0, method="Nelder-Mead", tol=1.0e-17)

    # compute max(err)
    A, b, w, _ = get_w_from_x(out.x)
    print("cond in solution:", numpy.linalg.cond(A))
    max_err = numpy.max(numpy.abs(A @ w - b))

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
    return d, max_err


def main():
    import json

    import fjson

    args = _get_parser().parse_args()

    with open(args.infile, "r") as f:
        content = json.load(f)

    new_data, max_err = optimize(content)
    if "weight factor" in content:
        w = content["weight factor"]
        for key, item in new_data.items():
            new_data[key][0] = (numpy.array(item)[0] / w).tolist()

    name = content["name"]
    prev_tol = content["test_tolerance"]
    if max_err < content["test_tolerance"]:
        content["data"] = new_data
        content["test_tolerance"] = max_err
        with open(args.infile, "w") as f:
            fjson.dump(content, f, indent=2, float_format=".15e")
            # for POSIX compliance:
            f.write("\n")
        print(f"{name}: Improved max error from {prev_tol} to {max_err}.")
    else:
        print(f"{name}: Could not improve scheme any further.")


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
