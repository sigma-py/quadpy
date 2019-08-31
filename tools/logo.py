import matplotlib.pyplot as plt

import quadpy

quadpy.quadrilateral.maxwell().plot()

# plt.xlim(-0.0, 1.0)
# plt.ylim(-0.0, 1.0)
# plt.show()
plt.savefig("out.svg", transparent=True, bbox_inches="tight")
