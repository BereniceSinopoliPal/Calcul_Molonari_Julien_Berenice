import numpy as np
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
import random as rd

col_dict = {
    "river_bed": 1, ##hauteur de la rivière en m
    "offset" : 0.05,#décalage
    "depth_sensores": [.1, .2, .3, .4], # Profondeur ou sont les capteurs
    "dH_mesures": list, # Decalage de profondeur des capteurs
    "temp_mesure": np.array, # shape (N,4,2) Chaque 4-uplets de mesure de temperature au temps t
    "sigma_meas_P": .4, #incertitude sur la pression
    "sigma_meas_T" : [3., 2., 4., 3.]
}

priors = {
    "moinslog10K": ((3, 10), 1), # (intervalle, sigma)
    "n": ...,
    "lambda_s": ...,
    "rhos_cs":...,
    "c_s":...
}

class Column:
    @classmethod
    def from_dict(cls, col_dict):
        return cls(**col_dict)
    def __init__(self, river_bed, offset, depth_sensors,dH_measures, temp_measure, sigma_meas_P, sigma_meas_T):
        self._dH = dH_measures
        self._T_mesure = temp_measure
        self._h = depth_sensors[-1]
        self._profondeur_mesure = depth_sensors
        self._dh = offset
        self._sigma_p = sigma_meas_P
        self._sigma_temp = sigma_meas_T
        self._rho_w = 1000
        self._c_w = 4180
        self.grad_H = []
        self.res_T = []
        self.run_mcmc = False
        self._t_mesure = []

        for i in range(len(self._dH)):
            self._t_mesure.append(self._dH[i][0])
        
        self.distrib_a_posteriori = None
        self.energie = None
        self.moy_acceptation = None
        self.profils_temp = None
        self.run_mcmc = False

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
        self.grad_H.append(np.asarray(delta_H))
        return np.asarray(delta_H)

    def solve_thermique(self, param: tuple, nb_cel: int, grad_h, alpha=0.7):
        K= param[0]
        lbds = param[1]
        n = param[2] ##normal ?
        pscs = param[3]

        dz = self._h/nb_cel
        lbdm = (n*np.sqrt(0.6)+(1-n)*np.sqrt(lbds))**2

        pmcm = n*self._rho_w*self._c_w + (1-n)*pscs


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
        K = 10**(-param['moinslog10K'])
        lbds = param['lambda_s']
        n = param['n']
        pscs = param['rhos_cs']
        nb_cel = param['nb_cel']

        delta_H = self.solve_hydro((K,n),nb_cel)

        res_temp= self.solve_thermique((K,lbds,n,pscs),nb_cel,delta_H)
        self.res_T.append(res_temp)

        return res_temp,delta_H

    def mcmc(self, priors: dict, nb_iter: int, nb_cel: int):

        self.run_mcmc = True 

        def pi(T_mesure, T_calcul, sigma_obs):
            T_mesure = np.array(T_mesure)
            T_calcul = np.array(T_calcul)
            return ((1/sigma_obs**6)*np.exp((-0.5/(sigma_obs**2))*np.linalg.norm(T_mesure - T_calcul)**2))

        def compute_energy(T_mesure, T_calcul, sigma_obs):
            T_mesure = np.array(T_mesure)
            T_calcul = np.array(T_calcul)
            return (-np.log((1/sigma_obs**6))*(-0.5/(sigma_obs**2))*np.linalg.norm(T_mesure - T_calcul)**2)

        def perturbation(borne_inf, borne_sup, previous_value, sigma):
            new_value = np.random.normal(previous_value, sigma)
            while new_value - borne_sup >0:
                new_value = borne_inf + (new_value - borne_sup)
            while borne_inf - new_value >0:
                new_value = borne_sup - (borne_inf - new_value)
            return new_value

       # def densite_rhos_cs(x, cs1=priors['c_s'][0][0], cs2=priors['c_s'][0][1], rho1=priors['rho_s'][0][0], rho2=priors['rho_s'][0][1]):
       #     if x < rho1*cs1 or x > rho2*cs2:
       #         return 0
       #     else:
       #         return (np.log(rho2*cs2/(rho1*cs1)) - abs(np.log(rho2*cs1/x)) - abs(np.log(rho1*cs2/x)))/(2*(rho2-rho1)*(cs2-cs1))


        #Calcul des indices de cellule correspondant à la profondeur des capteurs (on ne conserve pas ceux aux extrémités car ils servent pour les CL)

        indice_capteurs = np.rint(self._profondeur_mesure*nb_cel/self._h)
        indice_capteurs_interieur = indice_capteurs[1:4]


        #Initialisation des paramètres selon le prior et calcul des valeurs initiales
        moinslog10K_0 = np.random.uniform(priors['moinslog10K'][0][0], priors['moinslog10K'][0][1])
        lambda_s_0 = np.random.uniform(priors['lambda_s'][0][0], priors['lambda_s'][0][1])
        n_0 = np.random.uniform(priors['n'][0][0], priors['n'][0][1])
        rhos_cs_0 = np.random.uniform(priors['rhos_cs'][0][0], priors['rhos_cs'][0][1])
        param_0 = (moinslog10K_0, lambda_s_0, n_0, rhos_cs_0)

        dict_params_0 = {
            "moinslog10K": moinslog10K_0,
            "lambda_s": lambda_s_0,
            "n": n_0,
            "rhos_cs": rhos_cs_0,
            "nb_cel": nb_cel
        }
        
        T_mesure_0,*reste = self.solve_transi(dict_params_0)
        energie_init = compute_energy(self._T_mesure, [T_mesure_0[:,i] for i in indice_capteurs_interieur], self._sigma_temp)


        #Initialisation des tableaux de valeurs 
        
        params = [param_0] #Distribution a posteriori des paramètres (k, lambda_s, n, rho_s, c_s)
        all_dict_params = [dict_params_0]
        energie = [energie_init]
        profils_temp = [T_mesure_0] #Profils de température
        proba_acceptation = [] #Probabilité acceptation à chaque itération
        moy_acceptation = [] #Moyenne des probabilités d'acceptation à chaque itération
            
        for i in range(nb_iter):
            #Génération d'un état candidat

            moinslog10K_new = perturbation(priors['moinslog10K'][0][0], priors['moinslog10K'][0][1],params[-1][0], priors['moinslog10K'][1])
            lambda_s_new = perturbation(priors['lambda_s'][0][0], priors['lambda_s'][0][1],params[-1][1], priors['lambda_s'][1])
            n_new = perturbation(priors['n'][0][0], priors['n'][0][1], params[-1][2], priors['n'][1])
            rhos_cs_new = perturbation(priors['rhos_cs'][0][0], priors['rhos_cs'][0][1], params[-1][3], priors['rhos_cs'][1])
            param_new = (moinslog10K_new, lambda_s_new, n_new, rhos_cs_new)

            #Résolution du régime transitoire
            
            dict_params_new = {
                "moinslog10K": moinslog10K_new,
                "lambda_s": lambda_s_new,
                "n": n_new,
                "rhos_cs": rhos_cs_new,
                "nb_cel": nb_cel
            }
                

            T_res,*reste = self.solve_transi(dict_params_new) #verifier qu'on a bien un array en sortie

            #Calcul de la probabilité d'acceptation
            piX = pi(self._T_mesure, [profils_temp[-1][:,i] for i in indice_capteurs_interieur], self._sigma_temp)
            piY = pi(self._T_mesure, [T_res[:,i] for i in indice_capteurs_interieur], self._sigma_temp)
            

            if piX > 0:
                alpha = min(1, piY/piX)
            else :
                alpha = 1

            #Acceptation ou non
            if np.random.uniform(0,1) < alpha: #si le candidat est accepté
                params.append(param_new)
                all_dict_params.append(dict_params_new)
                profils_temp.append(T_res)
                energie.append(compute_energy(self._T_mesure, [T_res[:,i] for i in indice_capteurs_interieur], self._sigma_temp))
                proba_acceptation.append(alpha)
                moy_acceptation.append(np.mean([proba_acceptation[k] for k in range(i+1)]))

            else: #si le candidat n'est pas accepté, on reprend les valeurs précédentes dans les tableaux
                params.append(params[-1])
                all_dict_params.append(all_dict_params[-1])
                profils_temp.append(profils_temp[-1])
                energie.append(energie[-1])
                proba_acceptation.append(alpha)
                moy_acceptation.append(np.mean([proba_acceptation[k] for k in range(i+1)]))

        self.distrib_a_posteriori = params
        self.energie = energie
        self.moy_acceptation = moy_acceptation
        self.profils_temp = profils_temp

        """
        k_param = [params[i][0] for i in range(len(params))]
        plt.hist(k_param, bins=20)
        plt.show()"""
    
    #Ici la list les méthodes (non exhaustives)
    #pour recup les choses liées à la mcmc
    #la liste est vouée à évoluer.

    #@mcmc_needed
    def sample_param(self):
        a = np.random.randint(0, len(self.distrib_a_posteriori))
        sampled_param = self.distrib_a_posteriori[a] #vérifier la forme de distrib
        return sampled_param

    def get_all_params(self):
        return np.asarray(self.distrib_a_posteriori)

    def get_best_params(self):
        argmin_energie = np.argmin(self.energie)
        return self.distrib_a_posteriori[argmin_energie]

    def get_all_moinslog10K(self):
        return np.asarray(list(zip(*self.distrib_a_posteriori)))[0]

    def get_all_lambda_s(self):
        return np.asarray(list(zip(*self.distrib_a_posteriori)))[1]
    
    def get_all_lambda_n(self):
        return np.asarray(list(zip(*self.distrib_a_posteriori)))[2]

    def get_all_rhos_cs(self):
        return np.asarray(list(zip(*self.distrib_a_posteriori)))[3]

    def get_all_acceptance_ratio(self):
        return np.asarray(self.moy_acceptation)

    def get_temps_quantile(self, quantile): 
        temp_zip = list(zip(*self.profils_temp))
        profil_quantile = np.quantile(temp_zip, quantile, axis=1)
        return profil_quantile 

    def get_quantile_param(self, quantile):
        return np.quantile(self.distrib_a_posteriori, quantile, axis=0)

    #def get_flow_quantile(self, quantile):


        
