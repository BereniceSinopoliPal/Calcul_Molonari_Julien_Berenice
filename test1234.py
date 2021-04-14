import numpy as np

a = [(1,2,3), (45, 33, 22), (35, 45, 2)]

b = np.asarray(list(zip(*a)))

print(b[1])




