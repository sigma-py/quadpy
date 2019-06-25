# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

import quadpy

schemes = [
    [quadpy.triangle.BerntsenEspelid(k) for k in range(1, 5)],
    [quadpy.triangle.Centroid()],
    [quadpy.triangle.CoolsHaegemans(k) for k in [1]],
    [quadpy.triangle.Cubtri()],
    [quadpy.triangle.Dunavant(k) for k in range(1, 21)],
    [quadpy.triangle.Gatermann()],
    [quadpy.triangle.GrundmannMoeller(k) for k in range(7)],
    [quadpy.triangle.HammerMarloweStroud(k) for k in range(1, 6)],
    [quadpy.triangle.LaursenGellert(key) for key in [
      '1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
      '13', '14', '15a', '15b'
      ]],
    [quadpy.triangle.LiuVinokur(k) for k in range(1, 14)],
    [quadpy.triangle.LynessJespersen(k) for k in range(1, 22)],
    [quadpy.triangle.NewtonCotesClosed(k) for k in range(1, 6)],
    [quadpy.triangle.NewtonCotesOpen(k) for k in range(6)],
    [quadpy.triangle.SevenPoint()],
    [quadpy.triangle.Strang(k) for k in range(1, 11)],
    [quadpy.triangle.Stroud(k) for k in range(10)],
    [quadpy.triangle.TaylorWingateBos(k) for k in [1, 2, 4, 5, 8]],
    [quadpy.triangle.Triex(19), quadpy.triangle.Triex(28)],
    [quadpy.triangle.Vertex()],
    [quadpy.triangle.VioreanuRokhlin(k) for k in range(10)],
    [quadpy.triangle.Walkington(k) for k in [1, 2, 3, 5, 'p5']],
    [quadpy.triangle.WandzuraXiao(k) for k in range(1, 5)],
    [quadpy.triangle.WilliamsShunnJameson(k) for k in range(1, 9)],
    [quadpy.triangle.XiaoGimbutas(k) for k in range(1, 20)],
    [quadpy.triangle.ZhangCuiLiu(k) for k in [1, 2, 3]],
    ]

for scheme in schemes:
    # filter out schemes with negative weights
    p_scheme = [s for s in scheme if numpy.all(s.weights > 0.0)]
    # filter out schemes with points outside the triangle weights
    pi_scheme = [s for s in p_scheme if numpy.all(s.bary > 0.0)]
    print(scheme[0].name, len(scheme), len(pi_scheme))
    if pi_scheme:
        # get degrees and number of points
        degree = [s.degree for s in pi_scheme]
        num_points = [len(s.points) for s in pi_scheme]
        # plot
        plt.plot(
                degree, num_points,
                'o',
                label=pi_scheme[0].name
                )

plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
plt.xlabel('degree')
plt.ylabel('number of points')
plt.grid('on')
plt.show()
