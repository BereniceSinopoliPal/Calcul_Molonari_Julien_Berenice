### transitoire implem
import numpy as np
import matplotlib.pyplot as plt
def run_modele_direct():
    pass

p_init = np.linspace(1,2,100).reshape((100,1))
print(p_init)
def run_k_solve(K, Ss ,dt, dz, Zi, z0, zn, alpha = 0.7):

    n= Zi.shape[0]
    
    A=np.zeros((n,n))
    B=np.zeros((n,n))
    A[0][0]= 1
    A[1][0]= 2*alpha *K*4/(3*(dz)**2) ##terme limite
    A[1][1]= -alpha*K*4/(dz)**2
    A[1][2]= alpha *K*4/(3*(dz)**2)

    A[n-1][n-1]=1
    A[n-2][n-3]= (alpha)*K*4/(3*(dz)**2)
    A[n-2][n-2]= -(alpha)*K*4/(dz)**2
    A[n-2][n-1] = 2*(alpha)*K*4/(3*(dz)**2)

    
    B[0][0]=0
    B[1][0]= 2*(1-alpha)*K*4/(3*(dz)**2)
    B[1][1]= -(1-alpha)*K*4/(dz)**2
    B[1][2] = (1-alpha)*K*4/(3*(dz)**2)

    B[n-1][n-1]=0
    B[n-2][n-3]= (1-alpha)*K*4/(3*(dz)**2)
    B[n-2][n-2]= -(1-alpha)*K*4/(dz)**2
    B[n-2][n-1] = 2*(1-alpha)*K*4/(3*(dz)**2)

   
    
    for i in range(2,n-2):
        A[i][i-1]=K*alpha/dz**2
        A[i][i]=-(2*K*alpha/dz**2) - Ss/dt
        A[i][i+1]=K*alpha/dz**2

        B[i][i-1]=-K*(1-alpha)/dz**2
        B[i][i]=(2*K*(1-alpha)/dz**2) - Ss/dt
        B[i][i+1]=-K*(1-alpha)/dz**2

    C = np.matmul(B,Zi)
    C[0]= z0
    C[n-1]=zn

    res = np.linalg.solve(A,C)
    return res


for i in range(5):
    plt.scatter(p_init,range(len(p_init)),label='%s' %i)
    p_init = run_k_solve(1e-6,1e-5,30*60,1,p_init,1.5,1.7)
plt.legend()
plt.show()

