from matplotlib import pyplot as plt
import numpy as np

# mxp = -1
# mzp = -0.7
# alpha = 1

mxp = -1.000000;      mzp = 0.511609;         alpha=0.428170
mxp = 0.478753;       mzp = -1.000000;        alpha=0.017880

mx = abs(mxp)
mz = abs(mzp)

plt.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], "k")

x0 = np.clip((alpha)/mx, a_max=1, a_min=0) # x eval at y=0
x1 = np.clip((alpha-mz)/mx, a_max=1, a_min=0) # x eval at y=1

y0 = np.clip((alpha)/mz, a_max=1, a_min=0) # x eval at x=0
y1 = np.clip((alpha-mx)/mz, a_max=1, a_min=0) # x eval at x=1

plt.plot([x1, x0], [y0, y1], 'b:')

if mxp < 0 :
    x0 = 1- x0
    x1 = 1-x1

if mzp < 0 :
    y0 = 1-y0
    y1 = 1-y1

plt.plot([x1, x0], [y0, y1], 'r')
plt.show()

