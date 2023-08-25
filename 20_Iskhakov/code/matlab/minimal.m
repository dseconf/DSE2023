% This script is the 13 lines implementation of 
% The ENDOGENOUS GRIDPOINT METHOD
% to solve the stochastic dynamic model of 
% consumption and savings with credit constraint
% Written by Fedor Iskhakov, University of New South Wales, 2015

addpath('utils')

% Parameters
EXPN=10         ;% Number of quadrature points to calculate expectation
MMAX=10         ;% Maximum wealth
NM=100          ;% Number of grid points
TBAR=25         ;% Number of time periods
SIGMA=0.25      ;% Sigma parameter in logNormal distribution
Y=1             ;% Wage income
R=0.05          ;% Interest rate
DF=0.95         ;% Discount factor

% 13 lines EGM implementation
[quadp quadw]=quadpoints(EXPN,0,1);     % create quadrature notes and weights
quadstnorm=norminv(quadp,0,1);          % prepare quadrature points for calculation of expectations of Normal
savingsgrid=linspace(0,MMAX,NM);        % post-decision grid on savings
policy{TBAR}.w=[0 MMAX];                % terminal period wealth
policy{TBAR}.c=[0 MMAX];                % terminal period optimal consumption
for it=TBAR-1:-1:1                      % main backwards induction loop
 w1=Y+exp(quadstnorm*SIGMA)*(1+R)*savingsgrid;  % next period wealth (budget equation), matrix for all savings and all shocks
 c1=interp1(policy{it+1}.w,policy{it+1}.c,w1,'linear','extrap'); %next period optimal consumption
 rhs=quadw'*(1./c1);                    % RHS of the Euler equation (with log utility)
 policy{it}.c=[0 1./(DF*(1+R)*rhs)];    % current period optimal consumption rule
 policy{it}.w=[0 savingsgrid+policy{it}.c(2:end)]; % current period endogenous grid on wealth
end

% Plot the optimal policy functions
for it=TBAR:-1:1
	plot(policy{it}.w,policy{it}.c)
	hold all
end
set(gca,'XLim',[0 MMAX])
xlabel(gca,'Wealth')
ylabel(gca,'Optimal consumption')
title(gca,'Optimal consumption rules by age')


