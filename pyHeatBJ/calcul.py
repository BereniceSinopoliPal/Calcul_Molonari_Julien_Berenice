import numpy as np
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
import random as rd
import tqdm.auto as tqn
from numba import njit

##pylint


@njit
def solve_hydro(K, n, h, dH, len_t_mesure, t_mesure, dz, nb_cel, alpha=0.7):
    Ss = n / h

    H_array = np.zeros((len_t_mesure, nb_cel), dtype=np.float64)

    H_array[0] = np.linspace(dH[0], 0, nb_cel).astype(np.float64)

    for j in range(1, len_t_mesure):
        dt = t_mesure[j - 1]
        A = np.zeros((nb_cel, nb_cel), dtype=np.float64)
        B = np.zeros((nb_cel, nb_cel), dtype=np.float64)
        A[0][0] = 1
        A[1][0] = 8 * alpha * K / (3 * ((dz) ** 2))  # terme limite
        A[1][1] = -alpha * K * 4 / ((dz) ** 2) - Ss / dt
        A[1][2] = 4 * alpha * K / (3 * ((dz) ** 2))
        A[-1][nb_cel - 1] = 1
        A[nb_cel - 2][nb_cel - 3] = (alpha) * K * 4 / (3 * (dz) ** 2)
        A[nb_cel - 2][nb_cel - 2] = -(alpha) * K * 4 / ((dz) ** 2) - Ss / dt
        A[nb_cel - 2][nb_cel - 1] = 8 * (alpha) * K / (3 * (dz) ** 2)

        B[1][0] = -2 * (1 - alpha) * K * 4 / (3 * (dz) ** 2)
        B[1][1] = (1 - alpha) * K * 4 / (dz) ** 2 - Ss / dt
        B[1][2] = -(1 - alpha) * K * 4 / (3 * (dz) ** 2)
        B[nb_cel - 2][nb_cel - 3] = -(1 - alpha) * K * 4 / (3 * (dz) ** 2)
        B[nb_cel - 2][nb_cel - 2] = (1 - alpha) * K * 4 / (dz) ** 2 - Ss / dt
        B[nb_cel - 2][nb_cel - 1] = -8 * (1 - alpha) * K / (3 * (dz) ** 2)

        for i in range(2, nb_cel - 2):
            A[i][i - 1] = K * alpha / dz ** 2
            A[i][i] = -(2 * K * alpha / dz ** 2) - Ss / dt
            A[i][i + 1] = K * alpha / dz ** 2

            B[i][i - 1] = -K * (1 - alpha) / dz ** 2
            B[i][i] = (2 * K * (1 - alpha) / dz ** 2) - Ss / dt
            B[i][i + 1] = -K * (1 - alpha) / dz ** 2

        C = np.dot(B, H_array[j - 1])
        C[0], C[nb_cel - 1] = dH[j], 0

        res = np.linalg.solve(A, C)
        H_array[j] = res

    delta_H = np.zeros((len_t_mesure, nb_cel - 1), dtype=np.float64)

    for j in range(len_t_mesure):
        for p in range(len(H_array[j]) - 1):
            delta_H[j] = (H_array[j][p + 1] - H_array[j][p]) / dz
    return np.asarray(delta_H)


@njit
def solve_thermique(
    K,
    lbds,
    n,
    pscs,
    h,
    len_temp,
    t_mesure,
    nb_cel,
    grad_h,
    rho_w,
    c_w,
    profondeur_init,
    dH,
    T_mesure,
    alpha=0.7,
):

    dz = h / nb_cel
    lbdm = (n * np.sqrt(0.6) + (1 - n) * np.sqrt(lbds)) ** 2
    pmcm = n * rho_w * c_w + (1 - n) * pscs
    list_temp = np.zeros((len_temp, nb_cel), dtype=np.float64)
    list_temp[0] = profondeur_init.astype(np.float64)

    ke = lbdm / pmcm  ##lbm/pmcm
    ae = K * c_w * rho_w / pmcm  # K *pwcw/pmcm

    for j in range(1, len_temp):
        dt = t_mesure[j - 1]
        delta_H = grad_h[j]
        A = np.zeros((nb_cel, nb_cel), dtype=np.float64)
        B = np.zeros((nb_cel, nb_cel), dtype=np.float64)

        A[0][0] = 1
        A[nb_cel - 1][nb_cel - 1] = 1
        A[1][0] = alpha * (2 * ke / dz ** 2 - ae * delta_H[1] / (2 * dz))
        A[1][1] = alpha * (-2 * ke / dz ** 2) * (3 / 2) - 1 / dt
        A[1][2] = alpha * (ke / dz ** 2 + ae * delta_H[1] / (2 * dz))
        A[nb_cel - 2][nb_cel - 3] = alpha * (
            ke / dz ** 2 + ae * delta_H[nb_cel - 2] / (2 * dz)
        )
        A[nb_cel - 2][nb_cel - 2] = alpha * (-2 * ke / dz ** 2) * (3 / 2) - 1 / dt
        A[nb_cel - 2][nb_cel - 1] = alpha * (
            2 * ke / dz ** 2 - ae * delta_H[nb_cel - 2] / (2 * dz)
        )
        B[0][0] = 1
        B[nb_cel - 1][nb_cel - 1] = 1
        B[1][0] = -(1 - alpha) * (2 * ke / dz ** 2 - ae * delta_H[1] / (2 * dz))
        B[1][1] = -(1 - alpha) * (-2 * ke / dz ** 2) * (3 / 2) - 1 / dt
        B[1][2] = -(1 - alpha) * (ke / dz ** 2 + ae * delta_H[1] / (2 * dz))
        B[nb_cel - 2][nb_cel - 3] = -(1 - alpha) * (
            ke / dz ** 2 + ae * delta_H[nb_cel - 2] / (2 * dz)
        )
        B[nb_cel - 2][nb_cel - 2] = (
            -(1 - alpha) * (-2 * ke / dz ** 2) * (3 / 2) - 1 / dt
        )
        B[nb_cel - 2][nb_cel - 1] = -(1 - alpha) * (
            2 * ke / dz ** 2 - ae * delta_H[nb_cel - 2] / (2 * dz)
        )

        for i in range(2, nb_cel - 2):
            A[i][i - 1] = alpha * (ke / dz ** 2 - ae * delta_H[i] / (2 * dz))
            A[i][i] = alpha * (-2 * ke / dz ** 2) - 1 / dt
            A[i][i + 1] = alpha * (ke / dz ** 2 + ae * delta_H[i] / (2 * dz))

            B[i][i - 1] = -(1 - alpha) * (ke / dz ** 2 - ae * delta_H[i] / (2 * dz))
            B[i][i] = -(1 - alpha) * (-2 * ke / dz ** 2) - 1 / dt
            B[i][i + 1] = -(1 - alpha) * (ke / dz ** 2 + ae * delta_H[i] / (2 * dz))
        C = np.dot(B, list_temp[j - 1])
        C[0], C[nb_cel - 1] = dH[j], T_mesure[j][-1]
        res = np.linalg.solve(A, C)
        list_temp[j] = res

    return np.asarray(list_temp)


class Column:
    @classmethod
    def from_dict(cls, col_dict):
        return cls(**col_dict)

    def __init__(
        self,
        river_bed,
        offset,
        depth_sensors,
        dH_measures,
        T_measures,
        sigma_meas_P,
        sigma_meas_T,
    ):
        self._dH = dH_measures
        self._T_mesure = T_measures
        self._h = depth_sensors[-1]
        self._profondeur_mesure = depth_sensors
        self._dh = offset
        self._sigma_p = sigma_meas_P
        self._sigma_temp = sigma_meas_T
        self._rho_w = 1000
        self._c_w = 4180
        self._t_mesure = [i[0] for i in self._dH]
        self._T_mesure_int = [i[1][0:3] for i in T_measures]

        self.grad_H = []
        self.res_T = []
        self.debit = []
        self.distrib_a_posteriori = None
        self.energie = None
        self.moy_acceptation = None
        self.run_mcmc = False
        self.profil_temp_quantile = None
        self.param_quantile = None
        self.advec_flows = None
        self.conduc_flows = None
        self.debit_quantile = None
        self.flux_adv_quantile = None
        self.flux_cond_quantile = None

    def solve_transi(self, param: dict, alpha=0.7):
        K = 10 ** (-param["moinslog10K"])
        lbds = param["lambda_s"]
        n = param["n"]
        pscs = param["rhos_cs"]
        nb_cel = param["nb_cel"]
        lbdm = (n * np.sqrt(0.6) + (1 - n) * np.sqrt(lbds)) ** 2

        array_times = np.zeros((len(self._t_mesure)), dtype=np.float64)
        for j in range(1, len(self._t_mesure)):
            a = float((self._t_mesure[j] - self._t_mesure[j - 1]).total_seconds())
            array_times[j - 1] = a
        H_use = np.array([i[1][0] for i in self._dH])
        T_up_use = np.array([i[1][1] for i in self._dH])
        len_time = len(self._t_mesure)
        delta_H = solve_hydro(
            K, n, self._h, H_use, len_time, array_times, self._h / float(nb_cel), nb_cel
        )
        self.grad_H.append(np.asarray(delta_H))
        self.debit = [-K * i[0] for i in self.grad_H[-1]]

        Temp_use = np.asarray([i[1] for i in self._T_mesure], dtype=np.float64)

        coef = lagrange(
            [0] + self._profondeur_mesure, [self._dH[0][1][1]] + self._T_mesure[0][1]
        )
        profondeur = np.linspace([0], self._profondeur_mesure[-1], nb_cel)
        profondeur_inter = coef(profondeur)
        profondeur_inter = np.array([i[0] for i in profondeur_inter], dtype=np.float64)

        res_temp = solve_thermique(K,lbds,n,pscs,self._h,len_time,array_times,nb_cel,delta_H,self._rho_w,self._c_w,profondeur_inter,T_up_use,Temp_use)
        dz = self._h / nb_cel

        self.advec_flows = self._rho_w * self._c_w * delta_H * res_temp[:, :-1]
        self.conduc_flows = lbdm * np.gradient(
            res_temp, np.linspace(0, self._h, nb_cel), axis=-1
        )

        self.res_T.append(res_temp)
        return res_temp, delta_H

    def mcmc(self, priors: dict, nb_iter: int, nb_cel: int, quantile):

        self.run_mcmc = True

        def pi(T_mesure, T_calcul, sigma_obs, norm_init=1):
            T_mesure = np.array(T_mesure)
            T_calcul = np.array(T_calcul).transpose()
            return np.exp(
                (-0.5 / (sigma_obs ** 2)) * np.linalg.norm(T_mesure - T_calcul) ** 2
            )
        def compute_energy(T_mesure, T_calcul, sigma_obs):
            T_mesure = np.array(T_mesure)
            T_calcul = np.array(T_calcul).transpose()
            return (0.5 / (sigma_obs ** 2)) * np.linalg.norm(T_mesure - T_calcul) ** 2

        def perturbation(borne_inf, borne_sup, previous_value, sigma):
            new_value = np.random.normal(previous_value, sigma)
            while new_value - borne_sup > 0:
                new_value = borne_inf + (new_value - borne_sup)
            while borne_inf - new_value > 0:
                new_value = borne_sup - (borne_inf - new_value)
            return new_value

        # def densite_rhos_cs(x, cs1=priors['c_s'][0][0], cs2=priors['c_s'][0][1], rho1=priors['rho_s'][0][0], rho2=priors['rho_s'][0][1]):
        #     if x < rho1*cs1 or x > rho2*cs2:
        #         return 0
        #     else:
        #         return (np.log(rho2*cs2/(rho1*cs1)) - abs(np.log(rho2*cs1/x)) - abs(np.log(rho1*cs2/x)))/(2*(rho2-rho1)*(cs2-cs1))

        # Calcul des indices de cellule correspondant ?? la profondeur des capteurs (on ne conserve pas ceux aux extr??mit??s car ils servent pour les CL)

        indice_capteurs_interieur = np.rint(self._profondeur_mesure * nb_cel / self._h)[
            0:3
        ]

        # Initialisation des param??tres selon le prior et calcul des valeurs initiales
        moinslog10K_0 = np.random.uniform(
            priors["moinslog10K"][0][0], priors["moinslog10K"][0][1]
        )
        lambda_s_0 = np.random.uniform(
            priors["lambda_s"][0][0], priors["lambda_s"][0][1]
        )
        n_0 = np.random.uniform(priors["n"][0][0], priors["n"][0][1])
        rhos_cs_0 = np.random.uniform(priors["rhos_cs"][0][0], priors["rhos_cs"][0][1])
        param_0 = (moinslog10K_0, lambda_s_0, n_0, rhos_cs_0)

        dict_params_0 = {
            "moinslog10K": moinslog10K_0,
            "lambda_s": lambda_s_0,
            "n": n_0,
            "rhos_cs": rhos_cs_0,
            "nb_cel": nb_cel,
        }

        indice_capteurs_interieur = [int(i) for i in indice_capteurs_interieur]
        T_mesure_0, *reste = self.solve_transi(dict_params_0)
        debit_0 = self.get_flows_solve()
        advec_flow_0 = self.get_advec_flows_solve()
        cond_flow_0 = self.get_conduc_flows_solve()
        energie_init = compute_energy(
            self._T_mesure_int,
            [T_mesure_0[:, i] for i in indice_capteurs_interieur],
            self._sigma_temp,
        )
        # Initialisation des tableaux de valeurs

        params = [
            param_0
        ]  # Distribution a posteriori des param??tres (k, lambda_s, n, rho_s, c_s)
        all_dict_params = [dict_params_0]
        energie = [energie_init]
        profils_temp = [T_mesure_0]  # Profils de temp??rature
        debits = [debit_0]
        flux_adv = [advec_flow_0]
        flux_cond = [cond_flow_0]
        proba_acceptation = []  # Probabilit?? acceptation ?? chaque it??ration
        moy_acceptation = (
            []
        )  # Moyenne des probabilit??s d'acceptation ?? chaque it??ration

        for i in tqn.tqdm(range(nb_iter)):
            # for i in range(nb_iter):
            # G??n??ration d'un ??tat candidat

            moinslog10K_new = perturbation(
                priors["moinslog10K"][0][0],
                priors["moinslog10K"][0][1],
                params[-1][0],
                priors["moinslog10K"][1],
            )
            lambda_s_new = perturbation(
                priors["lambda_s"][0][0],
                priors["lambda_s"][0][1],
                params[-1][1],
                priors["lambda_s"][1],
            )
            n_new = perturbation(
                priors["n"][0][0], priors["n"][0][1], params[-1][2], priors["n"][1]
            )
            rhos_cs_new = perturbation(
                priors["rhos_cs"][0][0],
                priors["rhos_cs"][0][1],
                params[-1][3],
                priors["rhos_cs"][1],
            )
            param_new = (moinslog10K_new, lambda_s_new, n_new, rhos_cs_new)

            # R??solution du r??gime transitoire

            dict_params_new = {
                "moinslog10K": moinslog10K_new,
                "lambda_s": lambda_s_new,
                "n": n_new,
                "rhos_cs": rhos_cs_new,
                "nb_cel": nb_cel,
            }

            T_res, *reste = self.solve_transi(
                dict_params_new
            )  # verifier qu'on a bien un array en sortie

            # Calcul de la probabilit?? d'acceptation
            enX = compute_energy(
                self._T_mesure_int,
                [profils_temp[-1][:, i] for i in indice_capteurs_interieur],
                self._sigma_temp,
            )
            enY = compute_energy(
                self._T_mesure_int,
                [T_res[:, i] for i in indice_capteurs_interieur],
                self._sigma_temp,
            )
            rapport = np.exp(enX - enY)

            if True:
                alpha_accept = min(1, rapport)
            else:
                alpha_accept = 1

            # Acceptation ou non
            if np.random.uniform(0, 1) < alpha_accept:  # si le candidat est accept??
                params.append(param_new)
                all_dict_params.append(dict_params_new)
                profils_temp.append(T_res)
                energie.append(
                    compute_energy(
                        self._T_mesure_int,
                        [T_res[:, i] for i in indice_capteurs_interieur],
                        self._sigma_temp,
                    )
                )
                proba_acceptation.append(alpha_accept)
                moy_acceptation.append(
                    np.mean([proba_acceptation[k] for k in range(i + 1)])
                )
                debits.append(self.get_flows_solve())
                flux_adv.append(self.get_advec_flows_solve())
                flux_cond.append(self.get_conduc_flows_solve())

            else:  # si le candidat n'est pas accept??, on reprend les valeurs pr??c??dentes dans les tableaux
                params.append(params[-1])
                all_dict_params.append(all_dict_params[-1])
                profils_temp.append(profils_temp[-1])
                energie.append(energie[-1])
                proba_acceptation.append(alpha_accept)
                moy_acceptation.append(
                    np.mean([proba_acceptation[k] for k in range(i + 1)])
                )

                debits.append(debits[-1])
                flux_adv.append(flux_adv[-1])
                flux_cond.append(flux_cond[-1])

        self.distrib_a_posteriori = params
        self.energie = energie
        self.moy_acceptation = moy_acceptation
        self.proba_acceptation = proba_acceptation

        # Calcul des quantiles pour les param??tres
        self.param_quantile = np.quantile(self.distrib_a_posteriori, quantile, axis=0)

        # Calcul des quantiles pour les temp??ratures
        self.profil_temp_quantile = np.quantile(profils_temp, quantile, axis=0)

        # Calcul des quantiles pour le d??bit
        self.debit_quantile = np.quantile(debits, quantile)

        # Calcul des quantiles pour les flux thermiques
        self.flux_adv_quantile = np.quantile(flux_adv, quantile, axis=0)
        self.flux_cond_quantile = np.quantile(flux_cond, quantile, axis=0)

        # On r??initialise le tableau des profils de temp??rature, du d??bit et des fluxs pour ne pas les stocker en m??moire
        del profils_temp
        del debits
        del flux_adv
        del flux_cond

    def sample_param(self):
        a = np.random.randint(0, len(self.distrib_a_posteriori))
        sampled_param = self.distrib_a_posteriori[a]  # v??rifier la forme de distrib
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

    def get_temps_quantile(self):
        return self.profil_temp_quantile

    def get_moinslog10K_quantile(self):
        return self.param_quantile[0]

    def get_lambda_quantile(self):
        return self.param_quantile[1]

    def get_n_quantile(self):
        return self.param_quantile[2]

    def get_rhoscs_quantile(self):
        return self.param_quantile[3]

    def get_flows_solve(self):
        return self.debit

    def get_advec_flows_solve(self):
        return self.advec_flows

    def get_conduc_flows_solve(self):
        return self.conduc_flows

    def get_times_mcmc(self):
        return self._T_mesure[:, 1]

    def get_depths_mcmc(self):
        return self._profondeur_mesure
