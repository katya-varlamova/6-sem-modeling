
from numpy import e
import matplotlib.pyplot as plt
import numpy as np
from math import fabs
from mpl_toolkits.mplot3d import Axes3D


file = open("result.txt", "r")

x = file.readline()
x = (list(map(float, x.split())))
te = file.readline()
te = (list(map(float, te.split())))

file_line = True

arrs = []

while file_line:
    file_line = file.readline()

    if file_line:
        arrs.append(list(map(float, file_line.split())))


file.close()


arrs = np.array(arrs)


# 3д график
xgrid, ygrid = np.meshgrid(x, te)
xgrid = np.array(xgrid)
ygrid = np.array(ygrid)
zgrid = np.array(arrs)

fig = plt.figure(figsize=(7, 4))
ax_3d = Axes3D(fig)
ax_3d.plot_surface(ygrid, xgrid, zgrid, cmap='GnBu')
ax_3d.set_ylabel("x, cm")
ax_3d.set_xlabel("t, c")
ax_3d.set_zlabel("T, K")
plt.show()
