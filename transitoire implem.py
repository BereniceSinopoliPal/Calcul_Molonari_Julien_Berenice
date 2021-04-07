### transitoire implem
import numpy as np
import matplotlib.pyplot as plt
def run_modele_direct():
    pass

p_init = np.linspace(1,2,100).reshape((100,1))
print(p_init)
def run_k_solve(K, Ss ,dt, dz, Zi, z0, zn, alpha = 0.7):
    A=np.zeros((100,100))
    A[0][0]=1
    A[99][99]=1

    B=np.zeros((100,100))
    B[0][0]=0
    B[99][99]=0
    
    for i in range(1,99):
        A[i][i-1]=K*alpha/dz**2
        A[i][i]=-(2*K*alpha/dz**2) - Ss/dt
        A[i][i+1]=K*alpha/dz**2

        B[i][i-1]=-K*(1-alpha)/dz**2
        B[i][i]=(2*K*(1-alpha)/dz**2) - Ss/dt
        B[i][i+1]=-K*(1-alpha)/dz**2

    C = np.matmul(B,Zi)
    C[0]= z0
    C[99]=zn

    res = np.linalg.solve(A,C)
    return res


for i in range(5):
    plt.scatter(p_init,range(len(p_init)),label='%s' %i)
    p_init = run_k_solve(1e-6,1e-5,30*60,1,p_init,1.5,1.7)
plt.legend()
plt.show()

