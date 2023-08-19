clear all
% ************************************
% Model solution used for DGP
% ************************************
% Switches:

nMC=2; 	  	% number of MC samples
N=50;			% Number of busses to simulate 
T=119;			% Number of time periods to simulate 

% Read default parameters in to struct mp 
mp0.pnames_P={};   % set mp0.pnames_P={} to do PMLE (FULL MLE is not implemented for MPEC)
mp0.bellman_type='iv';
mp0=zurcher.setup(mp0);	


% ************************************
% START MONTE CARLO HERE
% ************************************

fprintf('Begin Monte Carlo, with n=%d replications\n', nMC);
rand('seed',300);

% solve model (used for simulation)
P0 = zurcher.statetransition(mp0);	% Transition matrix for mileage
u0 = zurcher.statetransition(mp0);	% Payoff function
bellman= @(ev) zurcher.bellman(ev, mp0);
[~, pk0]=dpsolver.poly(bellman, 0, mp0.ap, mp0.beta);	

mpec_results=struct;
nfxp_results=struct;

% starting values
mp=mp0;
% mp.RC=0;
% mp.c=0;


%main loop
for i_mc = 1:nMC;
	% ************************************
	% STEP 1: SIMULATE DATA 
	% ************************************
	timetosimulate=tic;
	data = zurcher.simdata(N, T, mp, P0, pk0);
	fprintf('i_mc=%d, Time to simulate data : %1.5g seconds\n', i_mc, toc(timetosimulate));

  	% ************************************
	% STEP 2a: ESTIMATE parameters using NFXP
	% ************************************
	[nfxp_result_i, theta_hat, Avar]=nfxp.estim(data, mp);
	nfxp_results=output.savemc(nfxp_result_i, nfxp_results, i_mc);

	% ************************************
 	% STEP 2b: ESTIMATE parameters using MPEC
 	% ************************************
	[mpec_result_i, theta_hat]=mpec.estim(data, mp);
	mpec_results=output.savemc(mpec_result_i, mpec_results, i_mc);

end  % End Monte Carlo

nfxp_results.title = 'Nested Fixed Point Algorithm, NFXP';
mpec_results.title = 'Mathematical Programming with Equilibrium Constraints (MPEC)';

results = {nfxp_results,mpec_results};

output.table_mc(mp.pnames_u, mp0, results); 
output.table_np(results); 
