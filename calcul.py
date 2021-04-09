col_dict = {
    "z_mesure": [.1, .2, .3, .4], # Profondeur ou sont les capteurs
	"t_mesure": list, #temps des mesures
    "delta_z": .05, # Decalage de profondeur des capteurs
    "p_mesure": np.array, # shape (N,) Chaque pression au temps t
    "temp_mesure": np.array, # shape (N,4) Chaque 4-uplets de mesure de temperature au temps t
    "sigma_p": .4, #incertitude sur la pression
    "sigma_temp" : [3., 2., 4., 3.]
    "rho_w": 3 #mettre la vraie valeur
    "c_w":2 #idem
}

priors = {
    "moinslog10K": ((3, 10), 1), # (intervalle, sigma)
    "n": ...,
    "lambda_s": ...,
	“rho_m_cm”: ...
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
        self._rhom_cm = prior['rho_m_cm']
        self._t_mesure = t_mesure
        self._rho_w = rho_w
        self._c_w = c_w
        self.distribution = None # [(k,lambda_s,n)]

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
        print(self._h)
        profondeur = np.linspace(0.1,0.4,nb_cel)
        profondeur_inter = coef(profondeur)
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

	def compute_mcmc(self, priors: dict, nb_iter: int) -> NoReturn:
		…
	
	
	#Ici la list les méthodes (non exhaustives)
	#pour recup les choses liées à la mcmc
	#la liste est vouée à évoluer.

	@mcmc_needed
	def sample_param(self):
	  #Tire un tuple de param au hasard parmis
		#ceux retenus avec la mcmc
		return param

	@mcmc_needed
	def get_best_params(self):
		return best_param

	@mcmc_needed
	def get_all_params(self):
		return params
	
	@mcmc_needed
	def get_all_moinslog10K(self):
		return moinslog10K_list

	@mcmc_needed
	def get_all_n(self):
		return n_list

	@mcmc_needed
	def get_all_lambda_s(self):
		return lambda_s_list

	@mcmc_needed
	def get_all_energy(self):
		return energy_list

	@mcmc_needed
	def get_all_acceptance_ratio(self):
		return acceptance_ratio_list