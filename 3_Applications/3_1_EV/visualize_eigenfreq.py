import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("build/eigenfrequencies_n10.txt")
data2 = np.loadtxt("build/eigenfrequencies_exact.txt")

plt.plot(data1[:,0],np.zeros(data1.shape[0]), "x")
plt.plot(data2[:,0],np.zeros(data2.shape[0]), "x")

plt.show()

input("..")