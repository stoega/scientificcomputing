import numpy as np

def Increase(a): 
    a += 1
    return

# simple type -> call by value -> a is not modified outside function body
a = 4
# list -> += operation not defined -> error
# a = [4]
# np.array -> call by reference -> a is modified
# a = np.array([4])


Increase(a)
print("a = ", a)