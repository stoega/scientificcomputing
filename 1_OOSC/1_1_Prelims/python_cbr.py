#%%

import numpy as np

#%%

def IncreaseByOne(a):
    a = 1
    print(f"increased, a = {a}")

#%%
n = np.zeros(3)
print(f"n = {n}")
IncreaseByOne(n)
print(f"after function call: n = {n}")

#%%