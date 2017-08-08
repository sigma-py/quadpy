'''
Import data from
https://arxiv.org/src/1411.5631v2/anc/fullsymmetry.txt
'''
import re
import numpy


def read_data(filename):
    data = []
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                # EOF
                break

            if line[0] == '#':
                continue

            if line[:6] == 'degree':
                pattern = (
                    'degree:\s+([0-9]+)\s+'
                    'points:\s+([0-9]+)\s+'
                    'orbits:\s+\[([0-9]+),\s*([0-9]+),\s*([0-9]+)\]\s+'
                    )
                out = re.match(pattern, line)
                degree = int(out.group(1))
                orbits = [int(out.group(k)) for k in [3, 4, 5]]
                num_items = sum(orbits)
                dat = numpy.fromfile(
                        f, count=num_items*4, sep=' '
                        ).reshape((num_items, 4))
                d0 = dat[:orbits[0]]
                d1 = dat[orbits[0]:orbits[0]+orbits[1]]
                d2 = dat[orbits[0]+orbits[1]:]
                data.append({
                    'degree': degree,
                    'data': {0: d0, 1: d1, 2: d2},
                    })

    return data


def data_to_code(data):
    for k, item in enumerate(data):
        print('elif index == {}:'.format(k))
        print('    self.degree = {}'.format(item['degree']))
        print('    data = [')

        for d0 in item['data'][0]:
            print(8*' ' + '({:.16e}, _s3()),'.format(d0[0]))

        for d1 in item['data'][1]:
            # find the value that appears twice
            if abs(d1[0] - d1[1]) < 1.0e-12:
                alpha = d1[0]
            else:
                alpha = d1[2]
            print(8*' ' + '({:.16e}, _s21({:.16e})),'.format(d1[0], alpha))

        for d2 in item['data'][2]:
            print(8*' ' + '({:.16e}, _s111({:.16e}, {:.16e})),'.format(
                d2[0], d2[1], d2[2]
                ))
        print('        ]')
    return


if __name__ == '__main__':
    data = read_data('papa.txt')
    data_to_code(data)
