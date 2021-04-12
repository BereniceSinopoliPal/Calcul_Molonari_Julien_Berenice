import numpy as np
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt

col_dict = {
    "z_mesure": [.1, .2, .3, .4], # Profondeur ou sont les capteurs
	"t_mesure": list, #temps des mesures
    "delta_z": .05, # Decalage de profondeur des capteurs
    "p_mesure": np.array, # shape (N,) Chaque pression au temps t
    "temp_mesure": np.array, # shape (N,4) Chaque 4-uplets de mesure de temperature au temps t
    "sigma_p": .4, #incertitude sur la pression
    "sigma_temp" : [3., 2., 4., 3.]
}

priors = {
    "moinslog10K": ((3, 10), 1), # (intervalle, sigma)
    "n": ...,
    "lambda_s": ...,
    "sigma_obs":...,
    "rho_m_cm":...
}

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
        self._rhom_cm = 4e6 ###Provisoire à modifier plus tard avec le MCMC
        self._t_mesure = t_mesure
        self._rho_w = 1000
        self._c_w = 4180
        self.distribution = None # [(k,lambda_s,n)]

    def solve_hydro(self, param: tuple, nb_cel: int, alpha=0.7):
        K= param[0]
        n = param[1]
        dz = self._h/nb_cel
        Ss = n/self._h

        list_P = [[] for i in range(len(self._t_mesure))]
        list_P[0] = np.linspace(self._dH[0],0,nb_cel)

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
            
            C = B @ list_P[j-1]
            C[0], C[nb_cel-1] = self._dH[j],0

            res = np.linalg.solve(A, C)
            list_P[j] = res
        
        delta_H=[[] for p in range(len(self._t_mesure))]
        for j in range(len(self._t_mesure)):    
            for p in range(len(list_P[j])-1):
                delta_H[j].append((list_P[j][p+1]-list_P[j][p])/dz)   
        return np.asarray(delta_H)

    def solve_thermique(self, param: tuple, nb_cel: int, grad_h, alpha=0.7):
        K= param[0]
        lbds = param[1]
        n = param[2] ##normal ?
        pmcm = param[3]

        dz = self._h/nb_cel
        lbdm = (n*np.sqrt(0.6)+(1-n)*np.sqrt(lbds))**2


        list_temp = [[] for i in range(len(self._t_mesure))]

        coef = lagrange(self._profondeur_mesure,self._T_mesure[0])
        profondeur = np.linspace(0.1,0.4,nb_cel)
        profondeur_inter = coef(profondeur)
        list_temp[0] = profondeur_inter

        ke = lbdm/pmcm##lbm/pmcm
        ae = K*(self._c_w*self._rho_w)/pmcm # K *pwcw/pmcm

        for j in range(1, len(self._t_mesure)):
            dt = self._t_mesure[j] - self._t_mesure[j-1]
            delta_H= grad_h[j]
            A=np.zeros((nb_cel,nb_cel))
            B=np.zeros((nb_cel,nb_cel))

            A[0][0]= 1
            A[nb_cel-1][nb_cel-1]=1
            A[1][0]=alpha*(2*ke/dz**2 - ae*delta_H[1]/(2*dz))
            A[1][1]=alpha*(-2*ke/dz**2)*(3/2) - 1/dt
            A[1][2]=alpha*(ke/dz**2 + ae*delta_H[1]/(2*dz))
            A[nb_cel-2][nb_cel-3]=alpha*(ke/dz**2 + ae*delta_H[nb_cel-2]/(2*dz))
            A[nb_cel-2][nb_cel-2]=alpha*(-2*ke/dz**2)*(3/2) - 1/dt
            A[nb_cel-2][nb_cel-1]=alpha*(2*ke/dz**2 - ae*delta_H[nb_cel-2]/(2*dz))


            B[0][0]= 1
            B[nb_cel-1][nb_cel-1]=1

            B[1][0]=-(1-alpha)*(2*ke/dz**2 - ae*delta_H[1]/(2*dz))
            B[1][1]=-(1-alpha)*(-2*ke/dz**2)*(3/2) - 1/dt
            B[1][2]=-(1-alpha)*(ke/dz**2 + ae*delta_H[1]/(2*dz))
            B[nb_cel-2][nb_cel-3]=-(1-alpha)*(ke/dz**2 + ae*delta_H[nb_cel-2]/(2*dz))
            B[nb_cel-2][nb_cel-2]=-(1-alpha)*(-2*ke/dz**2)*(3/2) - 1/dt
            B[nb_cel-2][nb_cel-1]=-(1-alpha)*(2*ke/dz**2 - ae*delta_H[nb_cel-2]/(2*dz))


            for i in range(2,nb_cel-2):
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
        list_temp=np.asarray(list_temp)
            
        return list_temp

    def solve_transi(self, param: dict, alpha=0.7):
        K = 10**(-param['K'])
        lbds = param['lambda_S']
        n = param['n']
        pmcm = param['Rho_m_C_m']
        nb_cel = param['nb_cel']

        delta_H = self.solve_hydro((K,n),nb_cel)

        res_temp= self.solve_thermique((K,lbds,n,pmcm),nb_cel,delta_H)

        return res_temp,delta_H

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
        energie_init = compute_energy(self._T_mesure, [T_mesure_0[:,i] for i in [0,33,66,99]], sigma_obs_0)


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

            T_res,*reste = self.solve_transi((moinslog10K_new, lambda_s_new, n_new), nb_cel, alpha) #verifier qu'on a bien un array en sortie

            #Calcul de la probabilité d'acceptation
            piX = pi(self._T_mesure, [profils_temp[-1][:,i] for i in [0,33,66,99]], sigma_obs_distrib[-1])
            piY = pi(self._T_mesure, [T_res[:,i] for i in [0,33,66,99]], sigma_obs_new)

            if piX > 0:
                alpha = min(1, piY/piX)
            else :
                alpha = 1

            #Acceptation ou non
            if np.random.uniform(0,1) < alpha: #si le candidat est accepté
                params.append((moinslog10K_new, lambda_s_new, n_new))
                sigma_obs_distrib.append(sigma_obs_new)
                profils_temp.append(T_res)
                energie.append(compute_energy(self._T_mesure, [T_res[:,i] for i in [0,33,66,99]], sigma_obs_new))
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
    
    #Ici la list les méthodes (non exhaustives)
    #pour recup les choses liées à la mcmc
    #la liste est vouée à évoluer.
"""
    #@mcmc_needed
    def sample_param(self):
    #Tire un tuple de param au hasard parmis
    # #ceux retenus avec la mcmc
        return param

    #@mcmc_needed
    def get_best_params(self):
        return best_param

    #@mcmc_needed
    def get_all_params(self):
        return params

    #@mcmc_needed
    def get_all_moinslog10K(self):
        return moinslog10K_list 

    #@mcmc_needed
    def get_all_n(self):
        return n_list

    #@mcmc_needed
    def get_all_lambda_s(self):
        return lambda_s_list

    #@mcmc_needed
    def get_all_energy(self):
        return energy_list

    #@mcmc_needed
    def get_all_acceptance_ratio(self):
        return acceptance_ratio_list"""