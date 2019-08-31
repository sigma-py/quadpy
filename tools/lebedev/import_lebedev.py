"""
This little helper takes Lebedev point and weight data from [1] and produces JSON files.

[1]
https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
"""
import os
import re

import numpy


def read(filename):
    data = numpy.loadtxt(filename)
    azimuthal_polar = data[:, :2] / 180.0
    weights = data[:, 2]
    return azimuthal_polar, weights


def chunk_data(weights):
    # divide the weight vector into chunks of 6, 8, 12, 24, or 48
    chunks = []
    k = 0
    ref_weight = 0.0
    tol = 1.0e-12
    while k < len(weights):
        if len(chunks) > 0 and abs(weights[k] - ref_weight) < tol:
            chunks[-1].append(k)
        else:
            chunks.append([k])
            ref_weight = weights[k]
        k += 1
    return chunks


def sort_into_symmetry_classes(weights, azimuthal_polar):
    data = {"a1": [], "a2": [], "a3": [], "pq0": [], "llm": [], "rsw": []}
    for c in chunks:
        if len(c) == 6:
            data["a1"].append([weights[c[0]]])
        elif len(c) == 12:
            data["a2"].append([weights[c[0]]])
        elif len(c) == 8:
            data["a3"].append([weights[c[0]]])
        elif len(c) == 24:
            if any(abs(azimuthal_polar[c, 1] - 0.5) < 1.0e-12):
                # polar == pi/2   =>   X == [p, q, 0].
                # Find the smallest positive phi that's paired with `polar ==
                # pi/2`; the symmetry is fully characterized by that phi.
                k = numpy.where(abs(azimuthal_polar[c, 1] - 0.5) < 1.0e-12)[0]
                assert len(k) == 8
                k2 = numpy.where(azimuthal_polar[c, 0][k] > 0.0)[0]
                azimuthal_min = numpy.min(azimuthal_polar[c, 0][k][k2])
                data["pq0"].append([weights[c[0]], azimuthal_min])
            else:
                # X = [l, l, m].
                # In this case, there must by exactly two phi with the value
                # pi/4. Take the value of the smaller corresponding `polar`;
                # all points are characterized by it.
                k = numpy.where(abs(azimuthal_polar[c, 0] - 0.25) < 1.0e-12)[0]
                assert len(k) == 2
                k2 = numpy.where(azimuthal_polar[c, 1][k] > 0.0)[0]
                polar_min = numpy.min(azimuthal_polar[c, 1][k][k2])
                data["llm"].append([weights[c[0]], polar_min])
        else:
            assert len(c) == 48
            # This most general symmetry is characterized by two angles; one
            # could take any two here.
            # To make things easier later on, out of the 6 smallest polar
            # angle, take the one with the smallest positive phi.
            min_polar = numpy.min(azimuthal_polar[c, 1])
            k = numpy.where(abs(azimuthal_polar[c, 1] - min_polar) < 1.0e-12)[0]
            k2 = numpy.where(azimuthal_polar[c, 0][k] > 0.0)[0]
            min_azimuthal = numpy.min(azimuthal_polar[c, 0][k][k2])
            data["rsw"].append([weights[c[0]], min_azimuthal, min_polar])

    return data


def write_json(filename, d):
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

    with open(filename, "w") as f:
        string = (
            pretty_floats(d)
            .__repr__()
            .replace("'", '"')
            .replace("{", "{\n  ")
            .replace("[[", "[\n    [")
            .replace("], [", "],\n    [")
            .replace(']], "', ']\n  ],\n  "')
            .replace("}", "\n}")
            .replace("]]", "]\n  ]")
        )
        f.write(string)

    return


if __name__ == "__main__":
    directory = "data/"
    for k, file in enumerate(os.listdir(directory)):
        filename = os.fsdecode(file)
        m = re.match("lebedev_([0-9]+)\\.txt", filename)
        degree = int(m.group(1))
        azimuthal_polar, weights = read(os.path.join("data", filename))
        chunks = chunk_data(weights)
        data = sort_into_symmetry_classes(weights, azimuthal_polar)

        delete_list = []
        for key in data:
            if len(data[key]) == 0:
                delete_list.append(key)
        for key in delete_list:
            data.pop(key)
        data["degree"] = degree

        write_json("lebedev_{:03d}.json".format(degree), data)
