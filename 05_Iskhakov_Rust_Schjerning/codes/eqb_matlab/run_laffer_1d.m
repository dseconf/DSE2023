%% Laffer curve
 
% This script creates the laffer curve with respect to the fuel tax,
% holding fixed car taxes at the counterfactual level (50% of the
% baseline).  %

%% SETUP 

close all; 
clear all; 
addpath('matlabinclude');
addpath('autotrade');

% set up parameters and indexes consistent with loaded data
load('results/estimation/estimation_data.mat');

% load parameter estimates
load('results/estimation/mle_converged.mat');

% main switch
DOSAVE = 1;
figOutPath = './results/laffer'; 
if ~exist(figOutPath, 'dir') && DOSAVE
    mkdir(figOutPath); 
    fprintf('Created dir, "%s"\n', figOutPath); 
end


%% STEP 0: baseline 
mp0=mp_mle;   % use mle estimates
mp0.modeltype = 'structuralform';
p0=sol_mle.p; % price vector at mle estimates

[sol0]=equilibrium.solve(mp0, s, p0); % solve model in baseline
[outcomes0] = stats.compute_outcomes(mp0, s, sol0); % comptue market outcomes 

%% STEP 1: Laffer iteration
% change in the fuel price in DKK/liter. Gets recalculated to the implied tax. Baseline = 0. 
num_points = 22; 
delta_fuelpList = linspace(-0.5*mp0.p_fuel,1.0*mp0.p_fuel, num_points)'; 

% multiplicative change in the new car tax rates (both high and low)
cartax_fraction = 1.0; 

h = waitbar(0, 'Evaluating'); 
outcomeList = cell(numel(delta_fuelpList), 1); % allocate space to save results 
for iDelta = 1:numel(delta_fuelpList)
    waitbar(iDelta/numel(delta_fuelpList), h, sprintf('Evaluating, %d/%d',iDelta,numel(delta_fuelpList)));
        
    mp1=mp0;
    delta_fuelp = delta_fuelpList(iDelta); 

    % set counterfactual taxes and compute implied prices
    mp1.tax_fuel    = (delta_fuelp + (mp0.p_fuel-mp0.p_fuel_notax)) / mp0.p_fuel_notax; % <--- implied change necessary in the fuel tax to achieve delta_fuelp 
    mp1.cartax_hi=mp0.cartax_hi*cartax_fraction;
    mp1.cartax_lo=mp0.cartax_lo*cartax_fraction;
    mp1=trmodel.update_mp(mp1); 

    [sol1]=equilibrium.solve(mp1, s, p0);
    p0=sol1.p; % update starting values for next iteration 

    [outcomes1] = stats.compute_outcomes(mp1, s, sol1);

    outcomeList{iDelta} = outcomes1; 
end

close(h) 

%% Car and fuel tax revenue separately 
xx = delta_fuelpList + mp0.p_fuel; 
% convert to fuel tax in pct. 
xx = 100.0 * (xx - mp0.p_fuel_notax) / mp0.p_fuel_notax; 

% loop and fill out tax revenues for each fuel price that we have solved
% above.
yy0 = nan(size(xx)); 
yy1 = nan(size(xx)); 
yy2 = nan(size(xx)); 
for i=1:numel(xx)
    yy0(i) = outcomeList{i}.total_revenue; 
    yy1(i) = outcomeList{i}.revenue_car; 
    yy2(i) = outcomeList{i}.revenue_fuel; 
end

figure(1); 
plot(xx,yy1,'-o', xx,yy2,'-x', xx,yy0,'-d');
xlabel('Fuel tax (%)'); ylabel('Total tax revenue'); 
xline(100, '--'); 
d=1.02; % offset 
text(100*d, outcomes0.total_revenue*0.95, 'Baseline tax'); 
%text(1, outcomes0.total_revenue*d, 'Baseline revenue'); 
set(gca, 'fontsize', 14);
leg = legend('Registration', 'Fuel', 'Total', 'Location', 'SouthEast'); 
title(leg, 'Tax'); 
if DOSAVE
    name_ = sprintf('%s/laffer.eps', figOutPath); 
    saveas(gcf, name_, 'epsc'); 
    fprintf('Figure saved as <a href="%s">%s</a>\n', figOutPath, name_);
end
