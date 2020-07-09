from matplotlib import pyplot as plt

import quadpy

scheme = quadpy.u2.krylov(30)
scheme.savefig("circle-krylov-30.svg")
plt.close()

scheme = quadpy.s2.hammer_stroud_20()
scheme.savefig("disk-hammer-stroud-20.svg")
plt.close()

scheme = quadpy.e1r.gauss_laguerre(3)
scheme.savefig("e1r-gauss-laguerre-3.svg")
plt.close()

scheme = quadpy.e1r2.gauss_hermite(8)
scheme.savefig("e1r2-gauss-hermite-8.svg")
plt.close()

scheme = quadpy.e2r.rabinowitz_richter_5()
scheme.savefig("e2r-rabinowitz-richter-5.svg")
plt.close()

scheme = quadpy.e2r2.rabinowitz_richter_3()
scheme.savefig("e2r2-rabinowitz-richter-3.svg")
plt.close()

scheme = quadpy.c1.gauss_legendre(20)
scheme.savefig("line-segment-gauss-legendre-20.svg")
plt.close()

scheme = quadpy.c2.maxwell()
scheme.savefig("quad-maxwell.svg")
plt.close()

scheme = quadpy.t2.dunavant_15()
scheme.savefig("triangle-dunavant-15.svg")
plt.close()
