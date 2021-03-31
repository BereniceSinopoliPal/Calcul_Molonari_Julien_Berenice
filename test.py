import numpy as np


class Params():
    def __init__(self, params_range, sigma = 0):
        if isinstance(params_range, list):
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

    #def run_modele_direct():
#### Exemple utilisation