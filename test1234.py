import numpy as np
 
a = [(1, 10, 100), (4, 40, 400), (5, 50, 500)]

print(np.quantile(a, [0.05, 0.95], axis=0))