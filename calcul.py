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
	 “rho_m_cm”: ...
}

class Column:
    	
	@classmethod
	def from_dict(cls, col_dict):
		…
		return column

	def solve_transi(self, param: tuple, nb_cel: int):
		…
		return (z_array, t_array), (temps, debit)

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