import numpy as np


class Params():
    def __init__(self, params_range, sigma = 0):
        if isinstance(params_range, list):
            self.down = params_range[0]
            self.up = params_range[1]
            self.sigma = sigma
        else :
            self.down=params_range
            self.up=params_range
            self.sigma = 0
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

    def pertub(self):
        self._value +=  np.random.randn(1)[0]*self._sigma
        ecart = self._up -self._down
        while self._value > self._up:
            self._value += -ecart
        while self._value < self._down:
            self._value += ecart
        return self._value

class Column():
    def __init__(h, nb_cellule, dt, data_profondeur, data_temperature, data_pression, rhom_cm=4e6, rho_w=1000, c_w=4180):
        self._h = h
        self._nb_cellule = nb_cellule
        self._dh = h/nb_cellule
        self._dt = dt
        self.rhom_cm=rhom_cm
        self.rho_w=rho_w
        self.c_w=c_w
        self.dH = data_pression
        self.T_mesure = data_temperature
        self.profondeur_mesure = data_profondeur

        self.distribution = None # [(k,lambda_s,n)]

    def run_MCMC(self, N, sigma_obs, k_param, lambda_s_param, n_param):

        def pi(T_mesure, k, lambda_s, n, sigma_obs):
            FY = np.array(modele_direct(moinslog10K, lambda_s, n)[0]) #attention ici à bien prendre les valeurs du modèle pour les points de mesure
            Z = np.array(T_mesure)
            return np.exp((-0.5/(sigma_obs**2))*np.linalg.norm(FY-Z)**2)

        def energie(T_mesure, k, lambda_s, n, sigma_obs):
            return -np.log(pi(T_mesure, moinslog10K, lambda_s, n, sigma_obs))


        #Initialisation des paramètres selon le prior

        k_0 = k_param.generate()
        lambda_s_0 = lambda_s_param.generate()
        n_0 = n_param.generate()

        # Application du modèle direct

        modele_direct_O = run_modele_direct(k_0, lambda_s_0, n_0)

        # Calcul de l'énergie initiale

        energie_0 = energie(T_mesure, k_0, lambda_s_0, n_0 )

        #Initialisation des tableaux de sortie

        distribution_a_posteriori = [[] for i in range(N+1)]
        distribution_a_posteriori[0] = [k_0, lambda_s_0, n_0]

        energie = [0 for i in range(N+1)]
        energie[0] = energie_O

        profils_temp = [[] for i in range(N+1)] #Profils de température
        profils_temp[0] = modele_direct_0

        proba_acceptation = [0 for i in range(N+1)]
        proba_acceptation[0] = 1

        moy_acceptation = [0 for i in range(N+1)]
        moy_acceptation[0] = 1

        #Chaine de calcul

        for i in range(N):
            #Génération d'un état candidat
            moinslog10K_new = k_param.perturb()
            lambda_s_new = lambda_s_param.perturb()
            n_new = n_param.perturb()

            #Calcul de la probabilité d'acceptation
            piX = pi(T_mesure, distribution_a_posteriori[i][0], distribution_a_posteriori[i][1], distribution_a_posteriori[i][2], sigma_obs)
            piY = pi(T_mesure, moinslog10K_new, lambda_s_new, n_new, sigma_obs)

            if piX > 0:
                alpha = min(1, piY/piX)
            else :
                alpha = 1

            #Acceptation ou non
            if np.random.uniform(0,1) < alpha: #si le candidat est accepté
                params[i+1] = [moinslog10K_new, lambda_s_new, n_new]
                modele_direct_i = run_modele_direct(moinslog10K_new, lambda_s_new, n_new)
                profils_temp[i+1] = modele_direct_i[1]
                energy[i+1] = energie(T_mesure, moinslog10K_new, lambda_s_new, n_new, sigma_obs)
                proba_acceptation[i+1] = alpha
                moy_acceptation[i+1] = np.mean([proba_acceptation[k] for k in range(i+1)])

            else: #si le candidat n'est pas accepté, on reprend les valeurs précédentes dans les tableaux
                params[i+1] = params[i]
                profils_temp[i+1] = profils_temp[i]
                energy[i+1] = energy[i]
                proba_acceptation[i+1] = alpha
                moy_acceptation[i+1] = np.mean([proba_acceptation[k] for k in range(i+1)])


        return(distribution_a_posteriori, energie, profils_temp, proba_acceptation, moy_acceptation)



    def run_modele_direct(k, lambda_s, n):


#### Exemple utilisation

colums = Columns(8, 100, 15*60, *)

k_param = Params(range, sigma = 0)
n_param =..
lambda_s_param =...
sigma_obs=....

colums.run_MCMC(nbsim,k,....)

columns.run_modele_direct()

