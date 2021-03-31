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
    def __init__(h, nb_cellule, dt, rhom_cm=4e6, rho_w=1000, c_w=4180):
        self._h = h
        self._nb_cellule = nb_cellule
        self._dh = h/nb_cellule
        self._dt = dt
        self.rhom_cm=rhom_cm
        self.rho_w=rho_w
        self.c_w=c_w

        self.distribution = None # [(k,lambda_s,n)]

    def run_MCMC():
        self.distribtuio = list

    def run_modele_direct():
        


#### Exemple utilisation

colums = Columns(8, 100, 15*60, *)

k = Params(range, sigma = 0)
n=..
lambda_s=...
sigma_obs=....

colums.run_MCMC(nbsim,k,....)

columns.run_modele_direct()

   