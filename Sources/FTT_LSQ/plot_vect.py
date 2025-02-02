import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("nvec.000", delimiter=" ", dtype=float)

plt.quiver(data[:, 0], data[:, 1], data[:, 2], data[:, 3])
plt.show()
