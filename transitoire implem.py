### transitoire implem
import numpy as np
import matplotlib.pyplot as plt

p_init = np.linspace(1,2,100).reshape((100,1))
#print(p_init) #réfléchir à passer la porosité plutot de Ss
def run_k_solve(K, nporo ,dt, dz, Zi, z0, zn, alpha = 0.7):
    
    n= Zi.shape[0]
    Ss = nporo/n

    A=np.zeros((n,n))
    B=np.zeros((n,n))

    A[0][0]= 1
    A[1][0]= 8*alpha*K/(3*((dz)**2)) ##terme limite
    A[1][1]= -alpha*K*4/((dz)**2) - Ss/dt
    A[1][2]= 4*alpha*K/(3*((dz)**2))

    A[n-1][n-1]=1
    A[n-2][n-3]= (alpha)*K*4/(3*(dz)**2)
    A[n-2][n-2]= -(alpha)*K*4/((dz)**2) - Ss/dt
    A[n-2][n-1] = 8*(alpha)*K/(3*(dz)**2)

    
    B[0][0]=1
    B[1][0]= -2*(1-alpha)*K*4/(3*(dz)**2)
    B[1][1]= (1-alpha)*K*4/(dz)**2 - Ss/dt
    B[1][2] = -(1-alpha)*K*4/(3*(dz)**2)

    B[n-1][n-1]=1
    B[n-2][n-3]= -(1-alpha)*K*4/(3*(dz)**2)
    B[n-2][n-2]= (1-alpha)*K*4/(dz)**2 - Ss/dt
    B[n-2][n-1] = -8*(1-alpha)*K/(3*(dz)**2)

   
    
    for i in range(2,n-2):
        A[i][i-1]=K*alpha/dz**2
        A[i][i]=-(2*K*alpha/dz**2) - Ss/dt
        A[i][i+1]=K*alpha/dz**2

        B[i][i-1]=-K*(1-alpha)/dz**2
        B[i][i]=(2*K*(1-alpha)/dz**2) - Ss/dt
        B[i][i+1]=-K*(1-alpha)/dz**2

    C = B @ Zi
    C[0]= z0
    C[n-1]=zn

    res = np.linalg.solve(A,C)
    return res

for i in range(5):
    #plt.scatter(p_init,range(len(p_init)),label='%s' %i)
    p_init = run_k_solve(1e-6,1e-3,15*60,1,p_init,1,1.1)
#plt.legend()
#plt.show()

def run_modele_solve(Mat_q, Temp_prev,lbdm, pmcm, dz, dt, T0, Tn, alpha=1):
    delta_H = np.zeros((100,1))
    for i in range(len(Mat_q)-1):
        #print((Mat_q[i+1]-Mat_q[i]))
        delta_H[i]=(Mat_q[i+1]-Mat_q[i])/dz
    #print(delta_H)
    n= Temp_prev.shape[0]
    A=np.zeros((n,n))
    B=np.zeros((n,n))

    A[0][0]= 1
    A[n-1][n-1]=1

    B[0][0]= 1
    B[n-1][n-1]=1
    ke = lbdm/pmcm
    K = 1e-6
    ae = K*(1000*4185)/pmcm


    for i in range(1,n-1):
        A[i][i-1]=alpha*(ke/dz**2 - ae*delta_H[i]/(2*dz))
        A[i][i]=alpha*(-2*ke/dz**2) - 1/dt
        A[i][i+1]=alpha*(ke/dz**2 + ae*delta_H[i]/(2*dz))

        B[i][i-1]=-(1-alpha)*(ke/dz**2 - ae*delta_H[i]/(2*dz))
        B[i][i]=-(1-alpha)*(-2*ke/dz**2) - 1/dt
        B[i][i+1]=-(1-alpha)*(ke/dz**2 + ae*delta_H[i]/(2*dz))


    C = B @ Temp_prev
    C[0]= T0
    C[n-1]=Tn

    res = np.linalg.solve(A,C)
    #print(res)
    return res

print(np.random.normal(0,1,100))
T_init = np.linspace(10+273,30+273,100).reshape((100,1)) 
for i in range(len(T_init)):
    T_init[i] += np.random.normal(0,1,1)[0]
print(T_init)
plt.close('all')
p_init = np.linspace(1,1.1,100)
#print(p_init)
for i in range(30):
    plt.plot(T_init,range(len(T_init)),label='%s' %i)
    T_init = run_modele_solve(p_init, T_init, 3, 4e5, 0.1, 24*60*60, 10+273+1*(-i), 30+273, alpha=0.7)
plt.legend()
plt.show() 
# fichier d'initialisation du projet à faire clean
# date time
# capacité calorifique
# travailler avec les paramètres équivalents ou non