# Calcul MOLONARI project Berenice and Julien
You can find in this repository the code of our project related to the DPE database ! 
You will find our code, corresponding to each step of our framework. The different steps are presented in the table of contents.

Don't forget to install all packages from requirements.txt by writing in your terminal : pip install -r requirements.txt


Enjoy :)

## Table of contents

- [Extract data from DPE databases and create base 0, base 1, base 2 and base 3](#Extract-data)
- [Bayesian Network for missing values inference](#Bayesian-Network-for-missing-values-inference)

## Extract data
You will find the code in the folder 'data_exploitation'.
We decided to split the code for the differents steps (extracting data, creating base 0, creating base 1...).
The main reason is to limit the use of the RAM (we could not open all the databases at the same time). 
If you want to run the code, donwload the dpe database at this link : and run the code in the same directory.
You have to run the code in the right order : from 1 to 9. 

The document "Notice DPE" explains our work.

ps: Depending on your computer, it can be long. On our laptop (Intel i5 - Quad core - 2.5 Ghz, 8GB RAM), it takes approximatively 24 hours. 


## Bayesian Network for missing values inference
You will find the code in the folder 'Bayesian network'.
There a two files in this folder, one file named 'Bayesian network implementation' that contains the code for the bayesian network applied to the database using structure learning, and one file named 'graphs' that contains the code for the visualization (graphs, heatmaps...).
