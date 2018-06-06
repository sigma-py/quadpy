# -*- coding: utf-8 -*-
#
"""
Parse Fortran code to extract points and weight of the Xiao-Gimbutas schemes.
"""
import numpy


def _parsed_strings_to_array(strings):
    return numpy.array(
        [
            val.replace(",", "").replace("&", "").replace("/", "").replace("D", "e")
            for val in strings
        ],
        dtype=float,
    )


def _parse():
    # The Fortran file contains multiple sections like
    # ```
    #   data xs / &
    #        0.00000000000000000000000000000000D+00/
    #        [...]
    #   data ys / &
    #        0.00000000000000000000000000000000D+00/
    #        [...]
    #   data ws / &
    #        0.21934566882541541013653648363283D+00/
    #        [...]
    # ```
    # Find those and extract the data.

    data = []

    with open("symq.txt", "r") as f:
        while True:
            line = f.readline()
            if not line:
                # EOF
                break

            line = line.strip()

            # skip if not at the start of a data block
            if line[:7] != "data xs":
                continue

            # start of a data block
            xstr = []
            while line[-1] == "&":
                line = f.readline().strip()
                xstr.append(line)

            line = f.readline().strip()
            assert line[:7] == "data ys"
            ystr = []
            while line[-1] == "&":
                line = f.readline().strip()
                ystr.append(line)

            line = f.readline().strip()
            assert line[:7] == "data ws"
            wstr = []
            while line[-1] == "&":
                line = f.readline().strip()
                wstr.append(line)

            points = numpy.column_stack(
                [_parsed_strings_to_array(xstr), _parsed_strings_to_array(ystr)]
            )
            weights = _parsed_strings_to_array(wstr)
            data.append((points, weights))

    return data


def _extract_bary_data(data):
    # The points are given in terms of coordinates of a reference triangle. Convert to
    # barycentric coordinates, and check their symmetry there.
    t0 = [-1, -1 / numpy.sqrt(3)]
    t1 = [+1, -1 / numpy.sqrt(3)]
    t2 = [0, 2 / numpy.sqrt(3)]

    T = numpy.array([[t1[0] - t0[0], t2[0] - t0[0]], [t1[1] - t0[1], t2[1] - t0[1]]])

    tol = 1.0e-10

    all_dicts = []

    ref_weight = 0.21934566882541541013653648363283

    for k, item in enumerate(data):
        points, weights = item

        b = (points - t0).T
        sol = numpy.linalg.solve(T, b)
        bary = numpy.column_stack([sol[0], sol[1], 1.0 - sol[0] - sol[1]])

        d = {"s1": [], "s2": [], "s3": [], "degree": k + 1}
        for w, b in zip(weights, bary):
            if numpy.all(numpy.abs(b - 1.0 / 3.0) < tol):
                weight = w / ref_weight
                d["s3"].append([weight])
            elif abs(b[0] - b[1]) < tol:
                weight = w / ref_weight / 3
                d["s2"].append([weight, b[0]])
            elif abs(b[1] - b[2]) < tol:
                weight = w / ref_weight / 3
                d["s2"].append([weight, b[1]])
            elif abs(b[2] - b[0]) < tol:
                weight = w / ref_weight / 3
                d["s2"].append([weight, b[0]])
            else:
                srt = numpy.sort(b)
                weight = w / ref_weight / 6
                d["s1"].append([weight, srt[0], srt[1]])

        for key in ["s1", "s2", "s3"]:
            if len(d[key]) == 0:
                d.pop(key)

        all_dicts.append(d)

    return all_dicts


def _main():
    data = _parse()
    all_dicts = _extract_bary_data(data)

    # Write the json files.

    # Getting floats in scientific notation in python.json is almost impossible, so do
    # some work here. Compare with <https://stackoverflow.com/a/1733105/353337>.
    class PrettyFloat(float):
        def __repr__(self):
            return "{:.16e}".format(self)

    def pretty_floats(obj):
        if isinstance(obj, float):
            return PrettyFloat(obj)
        elif isinstance(obj, dict):
            return dict((k, pretty_floats(v)) for k, v in obj.items())
        elif isinstance(obj, (list, tuple)):
            return list(map(pretty_floats, obj))
        return obj

    for d in all_dicts:
        degree = d["degree"]
        with open("xg{:02d}.json".format(degree), "w") as f:
            string = (
                pretty_floats(d)
                .__repr__()
                .replace("'", '"')
                .replace("{", "{\n  ")
                .replace("[[", "[\n    [")
                .replace("], [", "],\n    [")
                .replace(']], "', ']\n  ],\n  "')
                .replace("}", "\n}")
            )
            f.write(string)
    return


if __name__ == "__main__":
    _main()
