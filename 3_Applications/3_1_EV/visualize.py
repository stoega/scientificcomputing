import numpy as np
import matplotlib.pyplot as plt
import sys

print("pass filename: ", sys.argv[1])
filename = sys.argv[1]

data = np.loadtxt(filename, skiprows=1)
lam = np.loadtxt(filename, max_rows=1)

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
ax.set_title(f"lambda = {lam}")

plt.show()
