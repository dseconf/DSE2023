close all
clc
clear all;

mp.n=10;  % fixed point dimension

% Read default parameters in to struct mp 
mp=zurcher.setup(mp);

n_u = numel(struct2vec(mp,mp.pnames_u)); % n_u: number of preference parameters to be estimated (e.g. mp.RC and mp.c)
n_P = numel(struct2vec(mp,mp.pnames_P)); % n_P: number of free transition matrix parameters to be estimated 
N =mp.n;  			% N: number of EV parameters (i.e. number of grid points) 
M= n_P+1; 			% M: number transition matrix parameters. 

% ************************************************************************
% INSPECT SPARSITY PATTERNS 
% Jacobian of constraints and Hessian of Lagrangian for likelihood
% ************************************************************************
[J_pattern, H_pattern] = mpec.sparsity_pattern(n_u,n_P,N, M);

fprintf('Number of cost function parameters to be estimated:           n_u=%g\n', n_u);
fprintf('Number of free transition matrix parameters to be estimated   n_P=%g\n', n_P);
fprintf('Number of EV parameters (i.e. number of gridpoints)           N =%g\n', N);
fprintf('Number transition matrix parameters                           M =%g\n', M);
fprintf('Total number variables in constrained optimization problem:   n_u+n_P+N =%g\n', n_u+n_P+N);
fprintf('Number of non-zero elements in Jacobian =%g\n', sum(sum(J_pattern)));
fprintf('Number of non-zero elements in Hessian =%g\n', sum(sum(H_pattern)));

figure(1)
subplot(1,2,1), spy(J_pattern); 
set(gca,'FontSize',14)
title('Jacobian of constraints')

subplot(1,2,2), spy(H_pattern);
title('Hessian of likelihood');
set(gca,'FontSize',14)

% ************************************************************************
% INSPECT SPARSITY PATTERNS 
% Transition matrix for mileage
% ************************************************************************

% Cost function
cost=0.001*mp.c*mp.grid;				


%solve for fixed point in integrated value function space
mp.bellman_type='iv';
[V, pk, dV, iterinfo]=dpsolver.poly(@(ev) zurcher.bellman(ev, mp), 0, mp.ap, mp.beta);	

%solve for fixed point in expected value function space
mp.bellman_type='ev';
[ev, pk, dev, iterinfo]=dpsolver.poly(@(ev) zurcher.bellman(ev, mp), 0, mp.ap, mp.beta);	

% Transition matrix for mileage
P = zurcher.statetransition(mp);

figure(2)
subplot(1,2,1), spy(P{1})
set(gca,'FontSize',14)
title(sprintf('Transition matrix for mileage \n (keep)'));

subplot(1,2,2), spy(P{2})
set(gca,'FontSize',14)
title(sprintf('Transition matrix for mileage \n (replacement)'));

% dev is the Frechet derivative of bellman operator (which we need for the NK iteration)
% dev is equal to the discount factor time the unconditional (i.e. choice probability weighted) 
% transition matrix, Pu, for mileage, i.e. dev = mp.beta*Pu

figure(3)
subplot(1,2,1), spy(dev)
set(gca,'FontSize',14)
title(sprintf('Frechet derivative of ev=bellman_{ev}(ev) \n'));
subplot(1,2,2), spy(dV)
set(gca,'FontSize',14)
title(sprintf('Frechet derivative of V=bellman_{V}(V) \n'));

