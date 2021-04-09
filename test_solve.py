import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
class Column:

    @classmethod
    def from_dict(cls, col_dict):
        return cls(**col_dict)

    def __init__(self, z_mesure, t_mesure, delta_z, p_mesure, temp_mesure, sigma_p, sigma_temp):
        self._dH = p_mesure
        self._T_mesure = temp_mesure
        self._h = z_mesure[3]-z_mesure[0]
        self._profondeur_mesure = z_mesure
        self._dh = delta_z
        self._rhom_cm = 4e6
        self._t_mesure = t_mesure
        self._rho_w = 1000
        self._c_w = 4180
        self.distribution = None  # [(k,lambda_s,n)]

    def solve_transi(self, param: tuple, nb_cel: int, alpha=0.7):
        K= 10**(-param[0])
        lbdm = param[1]
        n = param[2]
        dz = self._h/nb_cel
        Ss = n/self._h

        list_P = [[] for i in range(len(self._t_mesure))]
        list_P[0] = np.linspace(*self._dH[0], nb_cel)
        
        for j in range(1, len(self._t_mesure)):
            dt = self._t_mesure[j] - self._t_mesure[j-1]
            A = np.zeros((nb_cel, nb_cel))
            B = np.zeros((nb_cel, nb_cel))
            A[0][0] = 1
            A[1][0] = 8*alpha*K/(3*((dz)**2))  # terme limite
            A[1][1] = -alpha*K*4/((dz)**2) - Ss/dt
            A[1][2] = 4*alpha*K/(3*((dz)**2))
            A[-1][nb_cel-1] = 1
            A[nb_cel-2][nb_cel-3] = (alpha)*K*4/(3*(dz)**2)
            A[nb_cel-2][nb_cel-2] = -(alpha)*K*4/((dz)**2) - Ss/dt
            A[nb_cel-2][nb_cel-1] = 8*(alpha)*K/(3*(dz)**2)

            B[1][0] = -2*(1-alpha)*K*4/(3*(dz)**2)
            B[1][1] = (1-alpha)*K*4/(dz)**2 - Ss/dt
            B[1][2] = -(1-alpha)*K*4/(3*(dz)**2)
            B[nb_cel-2][nb_cel-3] = -(1-alpha)*K*4/(3*(dz)**2)
            B[nb_cel-2][nb_cel-2] = (1-alpha)*K*4/(dz)**2 - Ss/dt
            B[nb_cel-2][nb_cel-1] = -8*(1-alpha)*K/(3*(dz)**2)

            for i in range(2, nb_cel-2):
                A[i][i-1] = K*alpha/dz**2
                A[i][i] = -(2*K*alpha/dz**2) - Ss/dt
                A[i][i+1] = K*alpha/dz**2

                B[i][i-1] = -K*(1-alpha)/dz**2
                B[i][i] = (2*K*(1-alpha)/dz**2) - Ss/dt
                B[i][i+1] = -K*(1-alpha)/dz**2
            #print(list_P[i-1])
            C = B @ list_P[j-1]
            C[0], C[nb_cel-1] = self._dH[j]

            res = np.linalg.solve(A, C)
            list_P[j] = res

        list_temp = [[] for i in range(len(self._t_mesure))]

        coef = lagrange(self._profondeur_mesure,self._T_mesure[0])
        profondeur = np.linspace(0.1,0.4,nb_cel)
        profondeur_inter = coef(profondeur)
        print(profondeur_inter)
        list_temp[0] = profondeur_inter

        ke = lbdm/self._rhom_cm##lbm/pmcm
        ae = K # K *pwcw/pmcm

        for j in range(1, len(self._t_mesure)):
            delta_H=[[] for p in range(len(list_P[j]))]
            for p in range(len(list_P[j])-1):
                delta_H[p] =  (list_P[j][p+1]-list_P[j][p])/dz
            A=np.zeros((nb_cel,nb_cel))
            B=np.zeros((nb_cel,nb_cel))

            A[0][0]= 1
            A[nb_cel-1][nb_cel-1]=1

            B[0][0]= 1
            B[nb_cel-1][nb_cel-1]=1

            for i in range(1,nb_cel-1):
                A[i][i-1]=alpha*(ke/dz**2 - ae*delta_H[i]/(2*dz))
                A[i][i]=alpha*(-2*ke/dz**2) - 1/dt
                A[i][i+1]=alpha*(ke/dz**2 + ae*delta_H[i]/(2*dz))

                B[i][i-1]=-(1-alpha)*(ke/dz**2 - ae*delta_H[i]/(2*dz))
                B[i][i]=-(1-alpha)*(-2*ke/dz**2) - 1/dt
                B[i][i+1]=-(1-alpha)*(ke/dz**2 + ae*delta_H[i]/(2*dz))
            C = B @ list_temp[j-1]
            C[0],C[nb_cel-1]= self._T_mesure[j][0],self._T_mesure[j][-1]
            res = np.linalg.solve(A,C)
            list_temp[j]=res
            
        return np.array(list_temp),list_P



col_dict = {
    "z_mesure": [.1, .2, .3, .4], # Profondeur ou sont les capteurs
	"t_mesure": np.arange(0,60*60,15*60), #temps des mesures
    "delta_z": .05, # Decalage de profondeur des capteurs
    "p_mesure": np.array([(0,0.1+np.random.normal(0,0.01,1)[0]) for i in range(4)]), # shape (N,) Chaque pression au temps t
    "temp_mesure": [[275.15,278.25,280.15,282.13],
[275.45,279.05,281.15,282.23],
[275.75,279.25,281.35,282.53],
[276.15,280.25,281.75,283.03]], # shape (N,4) Chaque 4-uplets de mesure de temperature au temps t
    "sigma_p": .4, #incertitude sur la pression
    "sigma_temp" : [3., 2., 4., 3.]
}



column = Column.from_dict(col_dict)
param=(4,2,0.15)
temp, pression = column.solve_transi(param, 100, alpha=0.7)
#print(res)
for i,j in enumerate(temp):
    print(j)
    plt.plot(j, np.linspace(0,0.4,100),label=f'{i}')
plt.legend()
plt.show()
