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

    def pertub


##pw cw lambda


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

    def run_MCMC(N, sigma_obs, k_param, lambda_s_param, n_param):

        def pi(T_mesure, k, lambda_s, n, sigma_obs):
            FY = np.array(modele_direct(moinslog10K, lambda_s, n)[0])
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

        proba_acceptation = [0 for i in range(N+1)] #Probabilité d'

        moy_acceptation = [0 for i in range(N+1)]



    def run_modele_direct(k, lambda_s, n):



#### Exemple utilisation

colums = Columns(8, 100, 15*60, *)

k_param = Params(range, sigma = 0)
n_param =..
lambda_s_param =...
sigma_obs=....

colums.run_MCMC(nbsim,k,....)

columns.run_modele_direct()

