# Replication code for "Equilibrium Trade in Automobiles"
### by Gillingham, Iskhakov, Munk-Nielsen, Rust & Schjerning
### Journal of Political Economy, 2022 
### [https://doi.org/10.1086/720463](https://doi.org/10.1086/720463)


### Research team:
- [Kenneth Gillingham](https://resources.environment.yale.edu/gillingham/) (Yale University)
- [Anders Munk-Nielsen](https://sites.google.com/view/munk-nielsen/) (University of Copenhagen)
- [John Rust](https://editorialexpress.com/jrust/) (Georgetown University) 
- [Fedor Iskhakov](https://fedor.iskh.me) (Australian National University)
- [Bertel Schjerning](https://bschjerning.com)  (University of Copenhagen)

## MATLAB scripts used for replication
The MATLAB script [run_all.m](https://github.com/fediskhakov/EQB/blob/master/run_all.m) reads in the data and replicates all the main results in the paper by calling the following MATLAB scripts in sequence

1. [run_illustrations.m](https://github.com/fediskhakov/EQB/blob/master/run_illustrations.m) 
    - MATLAB script that produces numerical illustrations of equilibrium in a stylized economy with 2 cars and 2 consumers 
    - Replicates Figure 2 and 3 in the theory section
    - Results stored in the folder results/illustration
1. [run_estimation.m](https://github.com/fediskhakov/EQB/blob/master/run_estimation.m): 
    - MATLAB script that reads in the tabulated data and runs DNFXP to obtain parameter estimates used below.
    - Results stored in the matfile results/estimation/mle_step2.mat
    - Before running the scripts below, make sure converged parameter estimates are stored as results/estimation/mle_converged.mat. In [run_all.m](https://github.com/fediskhakov/EQB/blob/master/run_all.m) this is done automatically by copying mle_step2.mat to mle_converged.mat.
1. [run_param_tables.m](https://github.com/fediskhakov/EQB/blob/master/run_param_tables.m)
    - MATLAB script that produces tex tables with estimation output.
    - Replicates tables with parameter estimates and standard errors in the online appendix.
    - Results stored in the folder results/tables
    - Note: This code also produces the summary statistics for all car types found in the oline appendix. To obtain the table with summary statistics for each of the 8 demographic groups, you must also execute the Jupyter Notebook [/data/sumstats.ipynb](https://github.com/fediskhakov/EQB/blob/master/data/sumstats.ipynb). 
1. [run_model_fit.m](https://github.com/fediskhakov/EQB/blob/master/run_model_fit.m) 
    -   MATLAB script that produces graphs with model fit 
    - Replicates Figure 4-6 in the empirical section of the paper
     - Results stored in the folder results/model_fit
1. [run_laffer_3d.m](https://github.com/fediskhakov/EQB/blob/master/run_laffer_3d.m)
    -   MATLAB script that produces contour plots and 3d-graphs with Tax-revenue (Laffer curve), CO2-emissions, and social surplus when varying the fuel and registration tax rates 
    - Replicates Figure 8 (contour plots) in empirical section of the paper and Figure 1 (3d-curves) in the online appendix. 
     - Results stored in the folder results/laffer
1. [run_iruc_CF.m](https://github.com/fediskhakov/EQB/blob/master/run_iruc_CF.m)
   - MATLAB script that produces counterfactual policy simulation results. This code run simulations assuming i) perfect pass-through (as in the paper) and ii) imperfect pass-through (not published)
   - Replicates Figure 7 (prices) and Table 1 (summary of policy simulation) in the empirical section of the paper.
   - Results (with perfect pass-through) stored in the folder results/iruc_CF.


