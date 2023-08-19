% run_mc_nfxp: Monte Carlo experiement to asses the performance of NFXP on Rust's engine repplacement model
clear all
% Read default parameters in to struct mp 
mp0.bellman_type='ev';

%% Adjustments to parameters 
mp0.pnames_P={'p'};   % estimate transition parameters  
%mp0.pnames_P={};   % skip estimation of transition parameters  

mp=zurcher.setup(mp0);		

% local constants:
nMC=2; 		% number of MC samples
N=50;			% Number of busses to simulate 
T=119;			% Number of time periods to simulate 

% Model solution used for DGP
bellman= @(ev) zurcher.bellman(ev, mp);
[~, pk0]=dpsolver.poly(bellman, 0, mp.ap, mp.beta);	

% Transition matrix for mileage
[P0] = zurcher.statetransition(mp);
[u0] = zurcher.u(mp);

fprintf('Begin Monte Carlo, with n=%d replications\n', nMC);

rand('seed',301);
nfxp_results=struct; 
for i_mc=1:nMC;
	% ************************************
	% SIMULATE DATA 
	% ************************************

	timetosimulate=tic;
	data = zurcher.simdata(N, T, mp, P0, pk0);
	timetosimulate=toc(timetosimulate);
	fprintf('i_mc=%d, Time to simulate data : %1.5g seconds\n', i_mc, timetosimulate);

	% ************************************
	% Estimate parameters and collect results
	% ************************************
	result_i=nfxp.estim(data, mp);
	nfxp_results=output.savemc(result_i, nfxp_results, i_mc);
end  % End Monte Carlo

nfxp_results.title = 'Monte Carlo results, NFXP';
output.table_mc([mp.pnames_u mp.pnames_P], mp, {nfxp_results}); 
output.table_np({nfxp_results}); 

