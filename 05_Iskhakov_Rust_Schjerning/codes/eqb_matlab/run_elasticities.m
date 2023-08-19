%% Calculate elasticities 
addpath('matlabinclude');
addpath('autotrade');

close all; 
clear all; 
% colormap(summer); 


% ****************************************************************************
% Parameters  
% ****************************************************************************
% mphat = setparams_denmark_jpe_submission_0; 
loaded = load('results/estimation/mle_converged.mat');

% parameters
mp0=loaded.mp_mle; % baseline 
sol0=loaded.sol_mle; 
s=trmodel.index(mp0); 

market_shares = @(sol) sum(sol.marketshares, 1); % sum over HH types 

%% evaluate baseline demands 
% Overall goal: two elasticities 
% - in equilibrium 
% - out of equilibrium 
%
% to compute in equilibrium, simply pick one of the car prices and shock it
% to see how the q vector (equilibrium holdings) responds 

if false 
    
    % extract starting values for prices
    global p0
    
    sol0 = equilibrium.solve(mp0, s, p0);
    q0 = market_shares(sol0);
    q = nan(mp0.ncartypes, mp0.ncartypes+1); % [car1, car2, ..., carJ, nocar]
    
    tic;
    fprintf('Evaluating price elasticities: ');
    for j = 1:mp0.ncartypes % the car to change
        fprintf('%d/%d ', j, mp0.ncartypes);
        mp1 = mp0;
        mp1.pnew_notax{j} = mp0.pnew_notax{j} * 1.01; % 1% increase
        sol1 = equilibrium.solve(mp1, s, p0);
        q(j, :) = market_shares(sol1);
    end
    fprintf(' (%5.2fs)\n', toc);
    
    %%
    fprintf('Price elasticities\n');
    e = (q - q0) ./ q0
    
    labsr = sprintfc('price of car %d', 1:mp0.ncartypes);
    labsc = [sprintfc('car %d', 1:mp0.ncartypes), 'no car'];
    array2table(e, 'VariableNames', labsc, 'RowNames', labsr)
    
    
end


%% partial equilibrium elasticities

global p0 % load starting values 
        equilibrium.solve(mp0, s, p0); % to get starting values close 
p0 = sol0.p; 
sol0  = equilibrium.solve(mp0, s, p0);
qnew0 = new_car_sales(sol0,mp0,s,sol0.q_tau); 

qnew = nan(mp0.ncartypes, mp0.ncartypes); 

tic; 
fprintf('Evaluating price elasticities: '); 
for j = 1:mp0.ncartypes % the car to change 
    fprintf('%d/%d ', j, mp0.ncartypes); 
    mp1 = mp0; 
    mp1.fixprices = 1; 
    p0 = sol0.p; 
    mp1.pnew_notax{j} = mp0.pnew_notax{j} * 1.01; % 1% increase in *pre*-tax 
    sol1 = equilibrium.solve(mp1, s, p0);
    qnew(j, :) = new_car_sales(sol1,mp1,s,sol0.q_tau);
end 
fprintf(' (%5.2fs)\n', toc);

% Compute elasticities 
fprintf('Partial EQ Price Elasticities\n'); 
ep = (qnew - qnew0) ./ qnew0 

% format as table 
labsr = sprintfc('price of car %d', 1:mp0.ncartypes);
labsc = [sprintfc('car %d', 1:mp0.ncartypes)];%, 'no car'];
T = array2table(round(ep,3), 'VariableNames', labsc, 'RowNames', labsr); 
T



%% Helper functions 

function qnew = new_car_sales(sol, mp, s, q_tau)
% OUTPUT: 
%   qnew: (1*ncartypes)
qnew = zeros(1, mp.ncartypes); 
for j=1:mp.ncartypes 
    for tau=1:mp.ntypes 
        Eccp = q_tau{tau}' * sol.ccp_tau{tau}; % 1*102 (num decisions) 
        ii = [s.id.trade_new{:}]; % rows of new car sales
        qnew = qnew + Eccp(ii); % 1*ncartypes 
    end
end
end
