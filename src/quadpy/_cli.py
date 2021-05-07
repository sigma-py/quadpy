import numpy as np

from ._optimize import optimize


def main():
    import json

    import fjson

    args = _get_parser().parse_args()

    with open(args.infile) as f:
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
                new_data[key][0] = (np.array(item)[0] / w).tolist()

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
        print(f"condition: {cond_in_solution:.3e}")


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
