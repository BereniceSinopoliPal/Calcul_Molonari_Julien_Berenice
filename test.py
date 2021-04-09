import numpy as np
from typing import Sequence
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange


class Params():
    def __init__(self, params_range, sigma = 0):
        if isinstance(params_range, Sequence):
            self._down = params_range[0]
            self._up = params_range[1]
            self._sigma = sigma
        else :
            self._down=params_range
            self._up=params_range
            self._sigma = 0

        self._value = None

    def generate(self):
        self._value = np.random.uniform(self._down,self._up) 
        return self._value

    def perturb(self):
        self._value +=  np.random.randn(1)[0]*self._sigma
        ecart = self._up -self._down
        while self._value > self._up:
            self._value += -ecart
        while self._value < self._down:
            self._value += ecart
        return self._value

class Column():
    def __init__(self, col_dict, dt, rho_w=1000, c_w=4180):
        self._h = col_dict['Column_Height']
        self._dH = col_dict['P_measures']
        self._T_mesure = col_dict['T_measures']
        self._profondeur_mesure = col_dict['Depth_Sensors']
        self._dh = col_dict['Delta_h']
        self._rhom_cm = col_dict['Rhom']
        self._nb_cellule = self._h/self._dh
        self._dt = dt
        self._rho_w=rho_w
        self._c_w=c_w
        self.distribution = None # [(k,lambda_s,n)]


    def run_MCMC(self, N, sigma_obs_param, k_param, lambda_s_param, n_param):

        def pi(T_mesure, moinslog10K, lambda_s, n, sigma_obs):
            res = self.run_modele_direct(moinslog10K, lambda_s, n)
            FY = np.array(res) #attention, vérifier ce que run_modele_direct renvoie en sortie
            Z = np.array(T_mesure)
            return np.exp((-0.5/(sigma_obs**2))*np.linalg.norm(FY-Z)**2)

        def energie_calcul(T_mesure, moinslog10K, lambda_s, n, sigma_obs):
            return -np.log(pi(T_mesure, moinslog10K, lambda_s, n, sigma_obs))


        #Initialisation des paramètres selon le prior et calcul des valeurs initiales
        k_0 = k_param.generate()
        lambda_s_0 = lambda_s_param.generate()
        n_0 = n_param.generate()
        sigma_obs = sigma_obs_param.generate()

        modele_direct_init = self.run_modele_direct(k_0, lambda_s_0, n_0)
        energie_init = energie_calcul(self._T_mesure, k_0, lambda_s_0, n_0, sigma_obs)

        #Initialisation des tableaux de sortie

        params = []
        distribution_a_posteriori = [[k_0, lambda_s_0, n_0]]
        energie = [energie_init]
        profils_temp = [modele_direct_init] #Profils de température
        proba_acceptation = [] #Probabilité acceptation à chaque itération
        moy_acceptation = [] #Moyenne des probabilitéd d'acceptation à chaque itération
            
        #Chaine de calcul

        for i in range(N):
            #Génération d'un état candidat
            moinslog10K_new = k_param.perturb()
            lambda_s_new = lambda_s_param.perturb()
            n_new = n_param.perturb()

            #res = modele_direct

            #Calcul de la probabilité d'acceptation
            piX = pi(self._T_mesure, distribution_a_posteriori[i][0], distribution_a_posteriori[i][1], distribution_a_posteriori[i][2], sigma_obs)
            piY = pi(self._T_mesure, moinslog10K_new, lambda_s_new, n_new, sigma_obs)

            if piX > 0:
                alpha = min(1, piY/piX)
            else :
                alpha = 1

            #Acceptation ou non
            if np.random.uniform(0,1) < alpha: #si le candidat est accepté
                params.append([moinslog10K_new, lambda_s_new, n_new])
                modele_direct_i = self.run_modele_direct(moinslog10K_new, lambda_s_new, n_new)
                profils_temp.append(modele_direct_i)
                energie.append(energie_calcul(self._T_mesure, moinslog10K_new, lambda_s_new, n_new, sigma_obs))
                proba_acceptation.append(alpha)
                moy_acceptation.append(np.mean([proba_acceptation[k] for k in range(i+1)]))

            else: #si le candidat n'est pas accepté, on reprend les valeurs précédentes dans les tableaux
                params.append(params[-1])
                profils_temp.append(profils_temp[-1])
                energie.append(energie[-1])
                proba_acceptation.append(alpha)
                moy_acceptation.append(np.mean([proba_acceptation[k] for k in range(i+1)]))


        return(distribution_a_posteriori, energie, profils_temp, proba_acceptation, moy_acceptation)


    
    def run_modele_direct(self,k, lambda_s, n):
        return True


#### Exemple utilisation

#colums = Columns(8, 100, 15*60, *)

k_param = Params([1,10], sigma = 1)
print(k_param.generate())
print(k_param.perturb())


