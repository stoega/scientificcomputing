import matplotlib.pyplot as plt
import numpy as np

xvals = np.zeros(5)
yvals = np.zeros(5)
xvals[0] = 0
xvals[1] = 1
xvals[2] = 2
xvals[3] = 3
xvals[4] = 5.1
yvals[0] = 2
yvals[1] = 2
yvals[2] = 6
yvals[3] = 7
yvals[4] = 7
xvals_fine = np.linspace(np.min(xvals), np.max(xvals), 100)

def f_reg(x):
	return 1.2639* x**0 + 2.63411* x**1 + -0.28888* x**2

plt.plot(xvals, yvals, 'x')
plt.plot(xvals_fine, f_reg(xvals_fine))
plt.show()
