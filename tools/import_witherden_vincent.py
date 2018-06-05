"""
Import data from Witherden/Vincent.
"""
import numpy
import re

import import_helpers


def read_data_tri(filename):
    data = numpy.loadtxt(filename)
    if len(data.shape) == 1:
        data = numpy.array([data])
    points = data[:, :2]
    weights = data[:, 2]
    # The reference triangle is (-1, -1), (1, -1), (-1, 1). Transform the
    # points to barycentric coordinates.
    points += 1.0
    points *= 0.5
    points = numpy.array(
        [points[:, 0], points[:, 1], 1.0 - numpy.sum(points, axis=1)]
    ).T
    return points, weights * 0.5


def read_data_tet(filename):
    data = numpy.loadtxt(filename)
    if len(data.shape) == 1:
        data = numpy.array([data])
    points = data[:, :3]
    weights = data[:, 3]
    # Transform to barycentric coordinates.
    points += 1.0
    points *= 0.5
    points = numpy.array(
        [points[:, 0], points[:, 1], 1.0 - numpy.sum(points, axis=1)]
    ).T
    return points, weights * 0.75


def data_to_code(points, weights):
    # identify groups of equal weights
    tol = 1.0e-12
    count = 0
    kk = 0
    last_value = weights[0]
    for w in weights:
        if abs(last_value - w) < tol:
            count += 1
        else:
            pts = points[kk : kk + count]
            kk += count
            print(
                8 * " "
                + "(%.15e, %s),"
                % (last_value, import_helpers.get_symmetry_code_tet(pts))
            )
            last_value = w
            count = 1

    pts = points[kk : kk + count]
    print(
        8 * " "
        + "(%.15e, %s)," % (last_value, import_helpers.get_symmetry_code_tet(pts))
    )
    return


def import_triangle():
    filenames = [
        "1-1.txt",
        "2-3.txt",
        "4-6.txt",
        "5-7.txt",
        "6-12.txt",
        "7-15.txt",
        "8-16.txt",
        "9-19.txt",
        "10-25.txt",
        "11-28.txt",
        "12-33.txt",
        "13-37.txt",
        "14-42.txt",
        "15-49.txt",
        "16-55.txt",
        "17-60.txt",
        "18-67.txt",
        "19-73.txt",
        "20-79.txt",
    ]
    for k, filename in enumerate(filenames):
        out = re.match("([0-9]+)-([0-9]+)\.txt", filename)
        strength = out.group(1)
        print("elif degree == {}:".format(strength))
        print("    data = [")
        x, weights = read_data_tri(filename)
        data_to_code(x, weights)
        print(8 * " " + "]")


def import_tet():
    filenames = [
        "1-1.txt",
        "2-4.txt",
        "3-8.txt",
        "5-14.txt",
        "6-24.txt",
        "7-35.txt",
        "8-46.txt",
        "9-59.txt",
        "10-81.txt",
    ]
    for k, filename in enumerate(filenames):
        out = re.match("([0-9]+)-([0-9]+)\.txt", filename)
        strength = out.group(1)
        print("elif degree == {}:".format(strength))
        print("    data = [")
        x, weights = read_data_tet(filename)
        data_to_code(x, weights)
        print(8 * " " + "]")


if __name__ == "__main__":
    import_triangle()
    # import_tet()
