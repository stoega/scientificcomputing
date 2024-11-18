import numpy as np
import matplotlib.pyplot as plt
import sys

## depending on platform/os
## e.g. build/poisson_solution.txt
## or   build/Release/poisson_solution.txt
print("pass filename: ", sys.argv[1])
filename = sys.argv[1]

data = np.loadtxt(filename)

nx = data.shape[0]
print("nx = ", nx)

x = np.linspace(0,1,nx)
y = np.linspace(0,1,nx)

meshx, meshy = np.meshgrid(x, y)


## simple plot
# plt.pcolormesh(meshx, meshy, data) 
# plt.show()


## 3D plot
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
surf = ax.plot_surface(meshx, meshy, data, cmap="jet")

plt.show()
