%% Replicate results for "Equilibrium Trade in Automobiles" 
% by 
% 	Kenneth Gillingham <kenneth.gillingham@yale.edu>, 
% 	Anders Munk-Nielsen <amn@econ.ku.dk> 
% 	John Rust <jrust@editorialexpress.com>
% 	Fedor Iskhakov <fediskhakov@gmail.com>
% 	Bertel Schjerning <bertel.schjerning@econ.ku.dk>
%
% This version: Revision for JPE on Dec, 2021

%% Illustrations for theory section
run_illustrations;  

%% Estimate model: 
run_estimation; 

% copy parameter estimates:  
copyfile('./results/estimation/mle_step2.mat','./results/estimation/mle_converged.mat'); 

% create tables for paper w. estimation output
run_param_tables; 

% run scrip that poroduces graphs with modelfit
run_model_fit; 

% counterfactuals
run_laffer_1d; 
run_laffer_3d; 
run_iruc_CF; 












