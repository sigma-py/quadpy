"""
Import data from
https://arxiv.org/src/1411.5631v2/anc/fullsymmetry.txt
"""
import re

import numpy as np


def read_data(filename, num_orbit_types):
    data = []
    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                # EOF
                break

            if line[0] == "#":
                continue

            orbit_pattern = (
                "orbits:\\s+\\["
                + ",\\s*".join(num_orbit_types * ["([0-9]+)"])
                + "\\]\\s+"
            )
            if line[:6] == "degree":
                pattern = (
                    "degree:\\s+([0-9]+)\\s+"
                    + "points:\\s+([0-9]+)\\s+"
                    + orbit_pattern
                )
                out = re.match(pattern, line)
                degree = int(out.group(1))
                orbits = [int(out.group(k)) for k in range(3, 3 + num_orbit_types)]
                num_items = sum(orbits)
                dat = np.fromfile(f, count=num_items * 4, sep=" ").reshape(
                    (num_items, 4)
                )
                d0 = dat[: orbits[0]]
                d1 = dat[orbits[0] : orbits[0] + orbits[1]]
                d2 = dat[orbits[0] + orbits[1] :]
                data.append({"degree": degree, "data": {0: d0, 1: d1, 2: d2}})

    return data


def data_to_code(data, f):
    for k, item in enumerate(data):
        print(f"elif index == {k}:")
        print("    self.degree = {}".format(item["degree"]))
        print("    data = {")

        if len(item["data"][0]) > 0:
            print(8 * " " + "'s3': [")
            for d0 in item["data"][0]:
                print(12 * " " + f"[{d0[0]:.16e}],")
            print(12 * " " + "],")

        # for d1 in item['data'][1]:
        #     # find the value that appears twice
        #     if abs(d1[0] - d1[1]) < 1.0e-12:
        #         alpha = d1[0]
        #     else:
        #         alpha = d1[2]
        #     print(8*' ' + '({:.16e}, {}({:.16e})),'.format(d1[0], f[1], alpha))

        if len(item["data"][1]) > 0:
            print(8 * " " + "'rot': [")
            for d2 in item["data"][1]:
                print(12 * " " + "[{:.16e}, {:.16e}, {:.16e}],".format(*d2))
            print(12 * " " + "],")

        print(8 * " " + "}")
    return


def data_to_json(data, f):
    for k, item in enumerate(data):
        filename = f"full{k:02d}.json"

        with open(filename, "w") as fh:
            fh.write("{\n")

            if len(item["data"][0]) > 0:
                fh.write('  "s3": [')
                for d0 in item["data"][0]:
                    fh.write(f"    [{d0[0]:+.16e}],\n")
                fh.write("  ],\n")

            if len(item["data"][1]) > 0:
                fh.write('  "s2": [\n')
                for d1 in item["data"][1]:
                    # find the value that appears twice
                    if abs(d1[0] - d1[1]) < 1.0e-12:
                        alpha = d1[0]
                    else:
                        alpha = d1[2]
                    fh.write(f"    [{d1[0]:.16e}, {alpha:.16e}],\n")
                fh.write("  ],\n")

            if len(item["data"][2]) > 0:
                fh.write('  "s1": [\n')
                for d1 in item["data"][2]:
                    fh.write("    [{:.16e}, {:.16e}, {:.16e}],\n".format(*d1[:-1]))
                fh.write("  ],\n")

            # for "rot"
            # if len(item['data'][1]) > 0:
            #     fh.write('  \"rot\": [\n')
            #     for d2 in item['data'][1]:
            #         fh.write('    [{:.16e}, {:.16e}, {:.16e}],\n'.format(*d2))
            #     fh.write('  ],\n')

            fh.write("}")
    return


if __name__ == "__main__":
    data = read_data("papanicolopulos.txt", num_orbit_types=3)
    # data_to_code(data, ('_s3', '_s21', '_s111'))
    data_to_json(data, ("_s3", "_s21", "_s111"))
    # data = read_data('rotationalsymmetry.txt', num_orbit_types=2)
    # # data_to_code(data, ('s3', 'rot'))
    # data_to_json(data, ('s3', 'rot'))
