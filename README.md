# Calcul MOLONARI project Berenice and Julien
You can find in this repository the code of our calcul part of the MOLONARI project !\
You will find our code, corresponding to the [API](https://docs.google.com/document/d/1MTkTiIvkihPglCRGPvlOzKta0icA32QY/edit).\
The different functions are presented in the table of contents.
If you need more details, you should read the [API](https://docs.google.com/document/d/1MTkTiIvkihPglCRGPvlOzKta0icA32QY/edit).

Don't forget to install all packages from requirements.txt by writing in your terminal : pip install -r requirements.txt


Enjoy :)

## Table of contents

- [Create your column with geometrical parameters](#Instance-a-column)
- [Find temperatures and flows](#Solve-transitoire)
- [Compute MCMC](#MCMC)
- [Get results and parameters from simulations](#Get_params)

## Instance a column
Before running a 'solve transitoire' or a 'mcmc', you have to create a column. To do so, create an instance of the Column class.
You can create a column by using the 'from_dict' method and create a dictionnary with the following keys :
col_dict = {\
    "river_bed": 1, ##hauteur de la rivière en m\
    "offset" : 0.05,#décalage\
    "depth_sensors": np.array([.2, .4, .6, .8]), # Profondeur ou sont les capteurs\
    "dH_measures": dH_test,#[datetime,[P,T]]\
    "T_measures": T_test, # shape (N,4,2) Chaque 4-uplets de mesure de temperature au temps t\
    "sigma_meas_P": .4, #incertitude sur la pression\
    "sigma_meas_T" : 0.1\
}

Than you can use instance.mcmc(params) or instance.solve_transi(params)
we will explain below how to use these two methods.

## Solve transitoire
You have already created your instance of the column class.
You want to solve your problem, given boundaries conditions. 
To use this methods, you also have to create a dictionnary with the following keys :

params ={\
'moinslog10K':-np.log(8e-4),\
'lambda_s': 1.08,\ 
'n':0.15,\
'rhos_cs':4e6,\
'nb_cel':100\
}

Then call :\
instance.solve_transi(params)


## MCMC
You have already created your instance of the column class.
You want to solve the inverse problem, and find the parameters, given measures. 
To use this methods, you also have to create a dictionnary with the following keys :

priors = {\
    "moinslog10K": ((3, 10), 0.2), # (intervalle, sigma)\
    "n": ((0.01,0.2),0.01),\
    "lambda_s": ((1,5),0.1),\
    "rhos_cs":((700*650,2800*1080),100000),\
}\
The priors are used to explore the space of the differents parameters and find states with a low energy.

Then call :\
instance.mcmc(priors)

## Get params ans results
The previous methods return nothing. It means that you have to use our custom methods to access any attributs. To access flows and temperatures profils, you have to use :\
instance.method()

The list of methods is :
sample_param()\
get_all_params()\
get_best_params()\
get_all_moinslog10K()\
get_all_lambda_s()\
get_all_lambda_n()\
get_all_rhos_cs()\
get_all_acceptance_ratio()\
get_temps_quantile()\
get_moinslog10K_quantile()\
get_lambda_quantile()\
get_n_quantile()\
get_rhoscs_quantile()\
get_flows_solve()\
get_advec_flows_solve()\
get_conduc_flows_solve()\
get_times_mcmc()\
get_depths_mcmc()
