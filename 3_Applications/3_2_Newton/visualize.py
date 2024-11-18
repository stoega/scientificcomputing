#%%


import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("data_20C.txt", skiprows=1)

fit = np.loadtxt("computed.txt")


plt.loglog(data[:,0], data[:,1], "xb", label="E'")
plt.loglog(fit[:,0], fit[:,1], "--b", label="E' fit")
plt.legend()
plt.show()
plt.loglog(data[:,0], data[:,2], "xr", label="E''")
plt.loglog(fit[:,0], fit[:,2], "--r", label="E'' fit")
plt.legend()
plt.show()



#%%