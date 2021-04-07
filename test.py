import numpy as np
from typing import Sequence


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
    def __init__(self,h, nb_cellule, dt, data_profondeur, data_temperature, data_pression, rhom_cm=4e6, rho_w=1000, c_w=4180):
        self._h = h
        self._nb_cellule = nb_cellule
        self._dh = h/nb_cellule
        self._dt = dt
        self._rhom_cm=rhom_cm
        self._rho_w=rho_w
        self._c_w=c_w
        self._dH = data_pression
        self._T_mesure = data_temperature
        self._profondeur_mesure = data_profondeur

        self.distribution = None # [(k,lambda_s,n)]

    def run_MCMC(self, N, sigma_obs_param, k_param, lambda_s_param, n_param):

        def pi(T_mesure, k, lambda_s, n, sigma_obs):
            FY = np.array(modele_direct(moinslog10K, lambda_s, n)[0])
            Z = np.array(T_mesure)
            return np.exp((-0.5/(sigma_obs**2))*np.linalg.norm(FY-Z)**2)

        def energie(T_mesure, T_compute):
            return T_mesure


        #Initialisation des paramètres selon le prior et calcul des valeurs initiales
        k_0 = k_param.generate()
        lambda_s_0 = lambda_s_param.generate()
        n_0 = n_param.generate()

        modele_direct_init = self.run_modele_direct(k_0, lambda_s_0, n_0)
        energie_init = energie( modele_direct_init, self._T_mesure )

        #Initialisation des tableaux de sortie


        distribution_a_posteriori = [[k_0, lambda_s_0, n_0]]
        energie = [energie_init]
        profils_temp = [modele_direct_init] #Profils de température
        proba_acceptation = [] #Probabilité acceptation

        for _ in range(N):
            new_k = k_param.perturb()
            new_lbd = lambda_s_param.perturb()
            new_n = n_param.perturb()

            res = modele_direct


    
    def run_modele_direct(self,k, lambda_s, n):
        return True
#### Exemple utilisation

#colums = Columns(8, 100, 15*60, *)

k_param = Params([1,10], sigma = 1)
print(k_param.generate())
print(k_param.perturb())


