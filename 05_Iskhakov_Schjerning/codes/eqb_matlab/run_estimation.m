%% Estimation with new data and endogenous scrappage 
% 2021-9-11

% clc; 
close all; 
clear all; 
addpath('matlabinclude');
addpath('autotrade');

years = 1996:2008;   % sample period: max 1996:2009

%% STEP 0:  Read data 
% Read in data from .xlsx files or prev. stored .mat files
% Type "help data.setup" for a list of available datasets and instructions on how to add new data sets
datadir='data/8x4/'; 

% set up parameters and indexes consistent with loaded data
LOADMATFILES = true; % set to false first time you run this program to generate matfiles 
[mp, s, dta, car, demo, fuelp, est]=data.setup(datadir, years, LOADMATFILES); 

%  Estimate model outcomes (ccps, holding distributions ect. ) based on tabulated choice/state/type data
[sol_d] = estim.outcomes(mp, s, dta); 

save('results/estimation/estimation_data.mat', 'mp', 's', 'dta', 'car', 'demo', 'fuelp', 'est', 'sol_d')

data.sumstat_table_households(dta, demo, mp, './results/tables/sumstats_demo.tex'); 

%% STEP 1: Estimate model
% If mle_step1.mat is previously saved you can uncomment this block start straight from the first stage estimates. 
mp = trmodel.setparams(mp); 

mp.pnames = { 'u_0', 'u_a', 'u_even','mum', 'psych_transcost', 'psych_transcost_nocar', 'acc_0', 'acc_a', 'sigma_s'}; 
mp.pnames = {mp.pnames{:}, 'tc_sale', 'tc_sale_even'};

mp.p0=getfields(mp, mp.pnames); % set starting values to parameters in mp

mp.p0.u_0     = repcell({3.5}, mp.ntypes, mp.ncartypes);
mp.p0.u_a     = repcell({-0.5}, mp.ntypes, mp.ncartypes);
mp.p0.u_a_sq  = repcell({0}, 1, mp.ncartypes);
mp.p0.u_even  = repcell({0}, mp.ntypes, mp.ncartypes);
mp.p0.acc_0   = repcell({-5}, 1, mp.ncartypes); % logit transformed parameter 
mp.p0.acc_a   = repcell({0}, 1, mp.ncartypes); % logit transformed parameter 
mp.p0.mum     = repcell({0.1}, mp.ntypes,1);
mp.p0.psych_transcost = repcell({5}, mp.ntypes, 1); 
mp.p0.tc_sale =0;
mp.p0.sigma_s =1;

% scrappage prices
% pscrap = pnew * 0.87^25: the DAF suggested used-car price for a 
% 25 y/o car of the corresponding car type. (In practice, there is 
% virtually no variation in beta across cartypes/ages.)
mp.pscrap_notax=num2cell([mp.pnew{:}].*car.beta(car.age==0)'.^ [mp.abar_j0{:}]);  % ... same formula, but based on tax inclusive price 

% update parameters
mp = trmodel.setparams(mp);

%% estimate model
[mp_mle, theta_mle, Avar_mle, sol_mle, pvec0_mle] = estim.mle(mp, s, dta, {'bhhh','bhhh','bhhh'}, [20,20,10]);
save('results/estimation/mle_step1.mat', 'mp_mle', 'sol_mle', 'Avar_mle')

%% STEP 2: Run again using different algorithms with estimates from STEP 1 as starting values
% We only do this to see if there are any futher improvements in likelihood. 
% Step 2 takes some additional time, and only gives a minor imporve in the mean likelihood at the 4 digit. 
load('results/estimation/mle_step1.mat')
mp = mp_mle; 
mp.p0=getfields(mp, mp.pnames); % set starting values to parameters in mp
[mp_mle, theta_mle, Avar_mle, sol_mle, pvec0_mle] = estim.mle(mp, s, dta, {'bhhh','bfgs','bhhh','bhhh','bhhh'}, [10,20,10,10,10]); % original 
save('results/estimation/mle_step2.mat', 'mp_mle', 'sol_mle', 'Avar_mle')


