import numpy as np
from matplotlib import pyplot as plt

depth = 6

data = np.loadtxt("./Examples/DATA/fgrd.000", usecols=[0, 1, 4])

print(data[0])

uniform_data = np.zeros((2**depth, 2**depth))
X = np.zeros((2**depth, 2**depth))
Y = np.zeros((2**depth, 2**depth))

for i in range(data.shape[0]):
    x = data[i, 0]
    y =  data[i, 1]
    kappa = data[i, 2]
    if 0 <= x <= 1:
         if 0 <= y <= 1:
              X[int(x*(2**depth)), int(y*(2**depth))] = x
              Y[int(x*(2**depth)), int(y*(2**depth))] = y
              uniform_data[int(x*(2**depth)), int(y*(2**depth))] = kappa

ax = plt.gcf().add_subplot(projection="3d")
ax.plot_surface(X, Y, uniform_data)
plt.show()
              

