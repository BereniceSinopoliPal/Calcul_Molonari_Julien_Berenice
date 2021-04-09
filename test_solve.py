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
            
        return list_temp,list_P

    
    def mcmc(self, priors: dict, nb_iter: int, nb_cel: int, alpha: float):
        def pi(T_mesure, T_calcul, sigma_obs):
            T_mesure = np.array(T_mesure)
            T_calcul = np.array(T_calcul)
            return ((1/sigma_obs**4)*np.exp((-0.5/(sigma_obs**2))*np.linalg.norm(T_mesure - T_calcul)**2))

        def compute_energy(T_mesure, T_calcul, sigma_obs):
            T_mesure = np.array(T_mesure)
            T_calcul = np.array(T_calcul)
            return ((1/sigma_obs**2)*(-0.5/(sigma_obs**2))*np.linalg.norm(T_mesure - T_calcul)**2)

        def perturbation(borne_inf, borne_sup, previous_value, sigma):
            new_value = np.random.normal(previous_value, sigma)
            while new_value - borne_sup >0:
                new_value = borne_inf + (new_value - borne_sup)
            while borne_inf - new_value >0:
                new_value = borne_sup - (borne_inf - new_value)
            return new_value


        #Initialisation des paramètres selon le prior et calcul des valeurs initiales
        k_0 = np.random.uniform(priors['moinslog10K'][0][0], priors['moinslog10K'][0][1])
        lambda_s_0 = np.random.uniform(priors['lambda_s'][0][0], priors['lambda_s'][0][1])
        n_0 = np.random.uniform(priors['n'][0][0], priors['n'][0][1])
        sigma_obs_0 = np.random.uniform(priors['sigma_obs'][0][0], priors['sigma_obs'][0][1])
        
        T_mesure_0,*reste = self.solve_transi((k_0, lambda_s_0, n_0), nb_cel, alpha)
        energie_init = compute_energy(self._T_mesure, T_mesure_0, sigma_obs_0)


        #Initialisation des tableaux de valeurs 
        
        params = [(k_0, lambda_s_0, n_0)] #Distribution a posteriori des paramètres (k, lambda_s, n)
        sigma_obs_distrib = [sigma_obs_0]
        energie = [energie_init]
        profils_temp = [T_mesure_0] #Profils de température
        proba_acceptation = [] #Probabilité acceptation à chaque itération
        moy_acceptation = [] #Moyenne des probabilités d'acceptation à chaque itération
            
        for i in range(nb_iter):
            #Génération d'un état candidat

            moinslog10K_new = perturbation(priors['moinslog10K'][0][0], priors['moinslog10K'][0][1],params[-1][0], priors['moinslog10K'][1])
            lambda_s_new = perturbation(priors['lambda_s'][0][0], priors['lambda_s'][0][1],params[-1][1], priors['lambda_s'][1])
            n_new = perturbation(priors['n'][0][0], priors['n'][0][1], params[-1][2], priors['n'][1])
            sigma_obs_new  = perturbation(priors['sigma_obs'][0][0], priors['sigma_obs'][0][1], sigma_obs_distrib[-1], priors['sigma_obs'][1])

            T_res = self.solve_transi((moinslog10K_new, lambda_s_new, n_new), nb_cel, alpha) #verifier qu'on a bien un array en sortie

            #Calcul de la probabilité d'acceptation
            piX = pi(self._T_mesure, profils_temp[-1][i], sigma_obs_distrib[-1])
            piY = pi(self._T_mesure, T_res, sigma_obs_new)

            if piX > 0:
                alpha = min(1, piY/piX)
            else :
                alpha = 1

            #Acceptation ou non
            if np.random.uniform(0,1) < alpha: #si le candidat est accepté
                params.append((moinslog10K_new, lambda_s_new, n_new))
                sigma_obs_distrib.append(sigma_obs_new)
                profils_temp.append(T_res)
                energie.append(compute_energy(self._T_mesure, T_res, sigma_obs_new))
                proba_acceptation.append(alpha)
                moy_acceptation.append(np.mean([proba_acceptation[k] for k in range(i+1)]))

            else: #si le candidat n'est pas accepté, on reprend les valeurs précédentes dans les tableaux
                params.append(params[-1])
                profils_temp.append(profils_temp[-1])
                energie.append(energie[-1])
                proba_acceptation.append(alpha)
                moy_acceptation.append(np.mean([proba_acceptation[k] for k in range(i+1)]))

        self.distribution = params

        k_param = [params[i][0] for i in range(len(params))]
        plt.hist(k_param, bins=20)
        plt.show()

        return None
    



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



priors = {
    "moinslog10K": ((3, 10), 1), # (intervalle, sigma)
    "n": ((0.01, 0.25), .06),
    "lambda_s": ((1,5), .5),
    "sigma_obs":((0.1, 2), 0.05),
    #"rho_m_cm":...
}


column = Column.from_dict(col_dict)
param=(4,2,0.15)
temp, pression = column.solve_transi(param, 100, alpha=0.7)
#print(res)
for i,j in enumerate(temp):
    #print(j)
    plt.plot(j, np.linspace(0,0.4,100),label=f'{i}')
plt.legend()
plt.show()

column.mcmc(priors, 100, 100, 0.7)
