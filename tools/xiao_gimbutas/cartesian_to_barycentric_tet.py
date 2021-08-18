"""
Generate code from cartesian coordinates.
"""
import mpmath
import numpy as np

# set precision
mpmath.mp.dps = 15

# reference tet:
t0 = [-1, -1 / mpmath.sqrt(3), -1 / mpmath.sqrt(6)]
t1 = [+0, +2 / mpmath.sqrt(3), -1 / mpmath.sqrt(6)]
t2 = [+1, -1 / mpmath.sqrt(3), -1 / mpmath.sqrt(6)]
t3 = [+0, +0, 3 / mpmath.sqrt(6)]


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# read data from file,
# e.g.,
# <https://people.sc.fsu.edu/~jburkardt/f_src/tetrahedron_arbq_rule/tetrahedron_arbq_rule.f90>
def read_data(filename, blocks="xyzw"):
    data = []
    current_block = None
    next_block = 0
    with open(filename) as f:
        for line in f:
            line = line.replace("&", "").replace("/", "").replace("D", "e").strip()

            numbers = line.split(",")
            # remove empty strings
            numbers = [number for number in numbers if number]

            if all([is_float(number) for number in numbers]):
                for number in numbers:
                    if current_block is None:
                        current_block = next_block
                        if current_block == 0:
                            data.append({})
                    if blocks[current_block] not in data[-1]:
                        data[-1][blocks[current_block]] = []
                    data[-1][blocks[current_block]].append(mpmath.mp.mpf(number))
            elif current_block is not None:
                next_block = (current_block + 1) % len(blocks)
                current_block = None
    return data


def find_duplicate_groups(lst, tol=1.0e-10):
    arr = np.array(lst)
    groups = []
    n = len(lst)
    is_handled = np.zeros(n, dtype=bool)
    for idx in range(n):
        if is_handled[idx]:
            continue
        dups = np.where(abs(arr[idx] - arr[idx:]) < tol)[0] + idx
        groups.append(dups)
        is_handled[dups] = True
    return groups


def tet_symmetries(groups):
    if len(groups) == 1:
        print(8 * " " + "_s4(),")
        multiplicity = 1
    elif len(groups) == 2:
        if len(groups[0]) == 2 and len(groups[1]) == 2:
            print(8 * " " + "_s22(%s)," % groups[0][0])
            multiplicity = 6
        else:
            assert len(groups[0]) == 3 or len(groups[1]) == 3
            if len(groups[0]) == 3:
                a = groups[0][0]
            else:
                a = groups[1][0]
            print(8 * " " + "_s31(%s)," % a)
            multiplicity = 4
    elif len(groups) == 3:
        assert len(groups[0]) == 2 or len(groups[1]) == 2 or len(groups[2]) == 2
        if len(groups[0]) == 2:
            a = groups[0][0]
            b = groups[1][0]
        elif len(groups[1]) == 2:
            a = groups[1][0]
            b = groups[2][0]
        else:
            a = groups[2][0]
            b = groups[0][0]
        print(
            (
                8 * " "
                + "_s211(\n"
                + 12 * " "
                + "%s,\n"
                + 12 * " "
                + "%s\n"
                + 12 * " "
                + "),"
            )
            % (a, b)
        )
        multiplicity = 12
    else:
        print(
            (
                8 * " "
                + "_s1111(\n"
                + 12 * " "
                + "%s,\n"
                + 12 * " "
                + "%s\n"
                + 12 * " "
                + "%s\n"
                + 12 * " "
                + "),"
            )
            % (lmbda[0], lmbda[1], lmbda[2])
        )
        multiplicity = 24

    return multiplicity


data = read_data("tet.txt")

for i, scheme_data in enumerate(data):
    X = scheme_data["x"]
    Y = scheme_data["y"]
    Z = scheme_data["z"]
    W = scheme_data["w"]

    if i == 0:
        print("if index == %d:" % (i + 1))
    else:
        print("elif index == %d:" % (i + 1))
    print("    bary = np.array([")

    # generate barycentric coordinate code
    XYZ = [[xx, yy, zz] for xx, yy, zz in zip(X, Y, Z)]
    T = mpmath.matrix([[t1[k] - t0[k], t2[k] - t0[k], t3[k] - t0[k]] for k in range(3)])
    multiplicities = []
    for k, xyz in enumerate(XYZ):
        b = [xyz[k] - t0[k] for k in range(3)]
        sol = mpmath.lu_solve(T, b)
        lmbda = [sol[0], sol[1], sol[2], 1.0 - sol[0] - sol[1] - sol[2]]
        assert abs(sum(lmbda) - 1.0) < 1.0e-10
        # print('%.15e %.15e %.15e' % (lmbda[0], lmbda[1], lmbda[2]))
        print(
            8 * " "
            + "[\n"
            + 12 * " "
            + f"{lmbda[0]}, {lmbda[1]},\n"
            + 12 * " "
            + f"{lmbda[2]}, {lmbda[3]},\n"
            + 8 * " "
            + "],"
        )
        # diffs = [
        #     abs(lmbda[k-1] - lmbda[k]) for k in range(4)
        #     ]
        # groups = find_duplicate_groups(diffs)
        # mult = tet_symmetries(groups)
        # multiplicities.append(mult)

    print("        ])")
    print("    self.weights = np.array([")
    # generate weight code
    alpha = mpmath.mp.mpf("0.9709835434146467")
    for weight in W:
        print("        %s," % (weight / alpha))
    print("        ])")

    # print('        ])')
    # print('    self.weights = np.concatenate([')
    # # generate weight code
    # alpha = mpmath.mp.mpf('0.9709835434146467')
    # for weight, m in zip(W, multiplicities):
    #     print(
    #         '        %s * np.ones(%d),'
    #         % (weight / alpha / m, m)
    #         )
    # print('        ])')
