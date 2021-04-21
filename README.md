# Calcul MOLONARI project Berenice and Julien
You can find in this repository the code of our calcul part of the MOLONARI project ! 
You will find our code, corresponding to the API (insert link to the API doc). The different functions are presented in the table of contents.
If you need more details, you should read the API

Don't forget to install all packages from requirements.txt by writing in your terminal : pip install -r requirements.txt


Enjoy :)

## Table of contents

- [Create your column with geometrical parameters](#Instance-a-column)
- [Find temperatures and flows](#Solve-transitoire)
- [Compute MCMC](#MCMC)
- [Get results and parameters from simulations](#Get_params)

## Instance a column
Before running a 'solve transitoire' or a 'mcmc', you have to create a column. To do so, create an instance of the Column class.
You can create a column by using the 'from_dict' method and create a dictionnary with the following index :
col_dict = {
    "river_bed": 1, ##hauteur de la rivière en m
    "offset" : 0.05,#décalage
    "depth_sensors": np.array([.2, .4, .6, .8]), # Profondeur ou sont les capteurs
    "dH_measures": dH_test,#[datetime,[P,T]]
    "T_measures": T_test, # shape (N,4,2) Chaque 4-uplets de mesure de temperature au temps t
    "sigma_meas_P": .4, #incertitude sur la pression
    "sigma_meas_T" : 0.1
}


## Solve transitoire
In process

## MCMC
In process

## Get params ans results
In process
