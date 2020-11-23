"""
Parse Fortran code to extract points and weight of the Xiao-Gimbutas schemes.
"""
import numpy


# TODO the first two functions could go into a helper and be shared with tri
def _parsed_strings_to_array(strings):
    return numpy.array(
        "".join(strings).replace("&", "").replace("/", "").replace("D", "e").split(","),
        dtype=float,
    )


def _parse():
    # The Fortran file contains multiple sections like
    # ```
    # data xs / &
    #   -.1685037180276000D+00,0.2783799427534418D-01, &
    #   [...]
    # data ys / &
    #   0.1910914916271708D+00,-.2304932838839657D-01, &
    #   [...]
    # data zs / &
    #   -.3896267314585163D+00,0.5481350663241830D+00, &
    #   [...]
    # data ws / &
    #   0.1287213727402025D+00,0.2179034339695993D+00, &
    #   [...]
    # ```
    # (Sometimes single columns.)
    # Find those and extract the data.
    data = []

    with open("tet.txt") as f:
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
            assert line[:7] == "data zs"
            zstr = []
            while line[-1] == "&":
                line = f.readline().strip()
                zstr.append(line)

            line = f.readline().strip()
            assert line[:7] == "data ws"
            wstr = []
            while line[-1] == "&":
                line = f.readline().strip()
                wstr.append(line)

            points = numpy.column_stack(
                [
                    _parsed_strings_to_array(xstr),
                    _parsed_strings_to_array(ystr),
                    _parsed_strings_to_array(zstr),
                ]
            )
            weights = _parsed_strings_to_array(wstr)
            data.append((points, weights))

    return data


def _extract_bary_data(data):
    # The points are given in terms of coordinates of a reference tetrahedron. Convert
    # to barycentric coordinates, and check their symmetry there.
    t0 = [-1, -1 / numpy.sqrt(3), -1 / numpy.sqrt(6)]
    t1 = [+0, +2 / numpy.sqrt(3), -1 / numpy.sqrt(6)]
    t2 = [+1, -1 / numpy.sqrt(3), -1 / numpy.sqrt(6)]
    t3 = [+0, +0, 3 / numpy.sqrt(6)]

    T = numpy.array([[t1[k] - t0[k], t2[k] - t0[k], t3[k] - t0[k]] for k in range(3)])

    all_dicts = []

    ref_weight = 0.9709835434146467

    for k, item in enumerate(data):
        d = {"degree": k + 1}
        points, weights = item

        b = (points - t0).T
        sol = numpy.linalg.solve(T, b)
        bary = numpy.column_stack(
            [sol[0], sol[1], sol[2], 1.0 - sol[0] - sol[1] - sol[2]]
        )

        idx = numpy.argsort(weights)
        d["weights"] = (weights[idx] / ref_weight).tolist()
        d["bary"] = bary[idx].tolist()
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
            return f"{self:.16e}"

    def pretty_floats(obj):
        if isinstance(obj, float):
            return PrettyFloat(obj)
        elif isinstance(obj, dict):
            return {k: pretty_floats(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple, numpy.ndarray)):
            return list(map(pretty_floats, obj))
        return obj

    for d in all_dicts:
        degree = d["degree"]
        print(d)
        with open(f"xg{degree:02d}.json", "w") as f:
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
