%% IRUC Counterfactual: exploring the role of the secondary market 
% 
% 2021-07-14
%
% Total runtime ~15min 
%
% STEP 0: baseline 
% STEP 1: new taxes + approximated price equilibrium => revenue neutrality
% STEP 2: new taxes + solve for equilibrium @ prices from STEP 2 
% STEP 3: new taxes + solve for equilibrium => revenue neutrality 

addpath('matlabinclude');
addpath('autotrade');

assert(isdir('results'), 'Folder "./results/" must be in present working directory.');
assert(isfile('results/estimation/mle_converged.mat'), 'Saved estimates file, "results/estimation/mle_converged.mat", not found.')

if ~isdir('results/iruc_CF')
    mkdir iruc_CF
end
if ~isdir('results/iruc_CF_imperfect')
    mkdir iruc_CF_imperfect
end

close all; 
clear all; 
colormap(summer); 

this_t = tic; 

% IMPERFECTPASSTHROUGH = true;
for  IMPERFECTPASSTHROUGH = 0:1 

fprintf('--- Running counterfactuals for IMPERFECTPASSTHROUGH = %d --- \n', IMPERFECTPASSTHROUGH);

% choose the decision variable and the outcome for the policy maker: the
% code will then choose the decision that makes the outcome identical in
% the baseline and the counterfactual
policy_outcome = 'total_revenue'; % 'total_revenue', 'total_co2', 'consumer_surplus', or any other field name in outcomes (assuming that equivalence is attainable using the policyvar) 
policyvar = 'tax_fuel'; % the decision variable for the policy maker 

if IMPERFECTPASSTHROUGH
    passthrough_rate = 0.9; 
    out_dir = 'results/iruc_CF_imperfect'; 
else
    passthrough_rate = 1.0; 
    out_dir = 'results/iruc_CF'; 
end

% paths (leave empty to drop saving)
outputfile      = sprintf('%s/counterfactual_%s_equivalence', out_dir, policy_outcome);
figpath         = out_dir;

ftol = 1e-5; % in Bln DKK or tonn Co2

% ****************************************************************************
% Parameters  
% ****************************************************************************
% mphat = setparams_denmark_jpe_submission_0; 
loaded = load('results/estimation/mle_converged.mat');

% parameters
mp0=loaded.mp_mle; % baseline 
sol0=loaded.sol_mle; 
s=trmodel.index(mp0); 

% set to structural form: 
mp0.modeltype   = 'structuralform'; 

% XXX: if VATFIRST changed wrt. estimated parameters, we have to update the
% prices ex tax. 
[a,b] = trmodel.price_notax(mp0); 

[~, mp0.pnew_notax, ~] = trmodel.price_notax(mp0); 

mp0 = trmodel.update_mp(mp0); 

% Verify integrity of pre-tax car prices
% if we have changed 
[a,b,c] = trmodel.price_notax(mp0); 
for j=1:mp0.ncartypes
    assert(b{j} == mp0.pnew_notax{j}, 'internal inconsistency in pre/post tax values! Did you change the tax system between estimation and counterfactuals?');
end


% counterfactual change multiplicative change in the new car tax rates (both high and low)
mp_cf=loaded.mp_mle; % cf = counterfactual
mp_cf.modeltype = 'structuralform'; 
mp_cf.cartax_hi=mp0.cartax_hi*.5; %0.5;
mp_cf.cartax_lo=mp0.cartax_lo*.5; %0.5;

if IMPERFECTPASSTHROUGH
    mp_cf.passthrough = trmodel.set_up_passthrough(mp0, passthrough_rate); 
end

mp_cf=trmodel.update_mp(mp_cf);

% ****************************************************************************
% STEP 0: baseline 
% ****************************************************************************

sol0 = equilibrium.solve(mp0, s, sol0.p); % solve model in baseline
outcomes0 = stats.compute_outcomes(mp0, s, sol0); % comptue market outcomes 

p0 = sol0.p; 

% ****************************************************************************
% STEP 1: Naive, expected: 
% what a policy maker would be lead to do if decisions were based on a
% model with proportional passthrough
% ****************************************************************************

% set parameters equal to baseline counter_factual
mp1              = mp_cf; 
mp1.fixprices    = 1; % circumvents the equilibrium solver

% rescale used-car prices proportionally to equilibrium prices in baseline
p1 = sol0.p; 
for j=1:mp0.ncartypes
    p1(s.ip{j}) = p1(s.ip{j}) * mp1.pnew{j} / mp0.pnew{j}; 
end

% solve for policy and update parameters
policy_objective1 = @(tax) dktax.policy_objective(mp1, p1, tax, outcomes0, policyvar, policy_outcome); 
mp1.tax_fuel = bisection(policy_objective1, 0.5, 4.0, ftol); 
mp1=trmodel.update_mp(mp1);
assert(not(isnan(mp1.tax_fuel)), 'Bisection failure!');

%% Solve model (holding price change proportional) and compute outcoms
sol1 = equilibrium.solve(mp1, s, p1); 
outcomes1 = stats.compute_outcomes(mp1, s, sol1); 






%% ****************************************************************************
% STEP 2: Naive, realized
% What the actual equilibrium looks like at the policy value chosen by the
% naive 
% ****************************************************************************

mp2 = mp_cf; 
mp2.tax_fuel = mp1.tax_fuel; % we use the policy choice from STEP 1
mp2 = trmodel.update_mp(mp2); 

% 3. solve model solving for equilibrium prices
mp2.fixprices=0;

% solve model and compute coutcomes 
sol2 = equilibrium.solve(mp2, s, p1); 
outcomes2 = stats.compute_outcomes(mp2, s, sol2); 


%% ****************************************************************************
% STEP 3: Sophisticated (accounting for equilibrium dynamics) 
% ****************************************************************************
% 1: set counterfactual taxes and compute implied prices
mp3 = mp_cf; 

% solve for policy and update parameters
policy_objective3 = @(tax) dktax.policy_objective(mp3, sol2.p, tax, outcomes0, policyvar, policy_outcome); 
mp3.tax_fuel = bisection(policy_objective3, .5, 4.0, ftol);
mp3=trmodel.update_mp(mp3);

% solve equilirbrium model and compute outcomes
sol3 = equilibrium.solve(mp3, s, p1); 
outcomes3 = stats.compute_outcomes(mp3, s, sol3); 



%% ****************************************************************************
% Print outcome comparison
% ****************************************************************************

titles = {'Baseline', 'Naive, expected', 'Naive, realized', 'Sophisticated'}; 
for DOSHORT=1
    for DOLATEX=0:1
        dktax.print_outcomes_comparison({outcomes0, outcomes1, outcomes2, outcomes3}, {mp0, mp1, mp2, mp3}, sprintf('RESULTS - Policy variable: %s, Outcome variable: %s', policyvar, policy_outcome), titles, outputfile, DOLATEX, DOSHORT); 
    end
    dktax.print_outcomes_comparison({outcomes0, outcomes1, outcomes2, outcomes3}, {mp0, mp1, mp2, mp3}, sprintf('RESULTS - Policy variable: %s, Outcome variable: %s', policyvar, policy_outcome), titles, [], false, DOSHORT); 
end

fprintf('--- Welfare changes from %s to %s --- \n', titles{1}, titles{4}); 
delta_welfare_tau = outcomes3.consumer_surplus_tau - outcomes0.consumer_surplus_tau; 
for tau=1:mp0.ntypes
    fprintf('%30s: %8.4f\n', mp0.lbl_types{tau}, delta_welfare_tau(tau)); 
end

% ONLY FOR DEBUGGING:
%% plot car prices
graphs.myfigure(); 
tiledlayout(2,2, 'TileSpacing', 'compact')
for j=1:mp0.ncartypes
    %subplot(3,2,j); 
    nexttile
    plot(...
     s.id.age(s.id.trade{j}), [mp0.pnew{j}; sol0.p(s.ip{j})], '-d', ... % baseline,   STEP 0
     s.id.age(s.id.trade{j}), [mp1.pnew{j}; sol1.p(s.ip{j})], '-o', ... % non-EQ,     STEP 1
     s.id.age(s.id.trade{j}), [mp2.pnew{j}; sol2.p(s.ip{j})], '-x', ... % EQ,         STEP 2 
     s.id.age(s.id.trade{j}), [mp3.pnew{j}; sol3.p(s.ip{j})], '-s');    % EQ, neutral STEP 3 
    ylabel('Price'); 
    title(sprintf('Car %d: %s', j, mp0.lbl_cartypes{j})); 
    xlabel('Car age'); 
    set(gca, 'fontsize', 14); set(gcf,'Color',[1 1 1]); set(gca, 'box', 'off', 'ygrid', 'on', 'ticklength', [0,0]); axis('tight'); 
    %graphs.set_fig_layout_post(gcf); 
end
lg = legend(titles, 'Location', 'southoutside', 'numcolumns', 4); 
lg.Position = [0.0875    0.0143    0.8277    0.0452]; 


if ~isempty(figpath)   
    name_ = sprintf('%s/prices_car%d_%s.eps', figpath, j, policy_outcome);
    saveas(gcf, name_, 'epsc');
    fprintf('Figure saved as <a href="%s">%s</a>\n', figpath, name_);
end

%% Uncomment to plot equilibrium
%plots={'prices', 'holdings', 'keep', 'taxes'}; % leave empty to do all plots
% plots = {'agg_holdings'};
% close all; 
% graphs.outcomes(mp0, s, plots, sol0, [], {'Model'})
% graphs.outcomes(mp1, s, plots, sol1, [], {'Model'})
% graphs.outcomes(mp2, s, plots, sol2, [], {'Model'})
% graphs.outcomes(mp3, s, plots, sol3, [], {'Model'})

% graphs.show(plots, mp0, s0, price_j0, q_tau0, ev_tau0, ccp_tau0); 
% graphs.show(plots, mp1, s1, price_j1, q_tau1, ev_tau1, ccp_tau1); 
% graphs.show(plots, mp2, s2, price_j2, q_tau2, ev_tau2, ccp_tau2); 
% graphs.show(plots, mp3, s3, price_j3, q_tau3, ev_tau3, ccp_tau3); 




%% --- Laffer Curves --- 
% Naive vs. Sophisticated Laffer curve
% this plots the Naive and Sophisticated Laffer curves, with =0 indicating
% "equal to the baseline." 

% 1. computations 
xx = linspace(0.0, 3., 12); 
yy1 = nan(size(xx)); 
yy3 = nan(size(xx)); 
h = waitbar(0, 'Computing Laffer curve points'); 
for i=1:numel(xx)
    waitbar((i-1)/numel(xx), h, 'Computing Laffer curve'); 
    yy1(i) = policy_objective1(xx(i)); 
    yy3(i) = policy_objective3(xx(i)); 
end
close(h); 

%% 2. Plot Naive vs. sophisticated Laffer curves 
figure
f=plot(xx, yy1, '-x', xx,yy3,'-or', 'LineWidth', 2.0);  xline(1.0, ':'); yline(0.0, ':'); 
legend('Naive', 'Equilibrium');
graphs.set_fig_layout_post(f); title('Naive, realized'); 
ylabel([policy_outcome ': CF - baseline'], 'Interpreter', 'none'); 
xlabel(policyvar, 'Interpreter', 'none'); 

if ~isempty(figpath)
    name_ = sprintf('%s/laffer_naive_vs_sophisticated.eps', figpath);
    saveas(gcf, name_, 'epsc');
    fprintf('Figure saved as <a href="%s">%s</a>\n', figpath, name_);
end


%% --- Outcomes ---

policy_outcome = 'total_revenue'; 
policyvar = 'tax_fuel'; 

tolx = 1e-3; 

num_points = 6; 
rates = linspace(0,1,num_points); 
tt = nan(size(rates)); 
outs = cell(numel(rates), 1);  

bounds0 = [1.0, 3.5]; 
opts = optimset('display', 'iter', 'tolx', tolx); 

h = waitbar(0, 'Computing'); 
for ir=1:numel(rates)
    waitbar((ir-1)/numel(rates),h,sprintf('Computing, %d/%d', ir, numel(rates))); 
    rate = rates(ir); 
    
    % 1: set counterfactual taxes and compute implied prices
    mp_cf=loaded.mp_mle; % cf = counterfactual
    mp_cf.modeltype = 'structuralform';
    mp_cf.cartax_hi=mp0.cartax_hi*rate;
    mp_cf.cartax_lo=mp0.cartax_lo*rate;
    mp_=trmodel.update_mp(mp_cf);
    
    % solve for policy and update parameters
    policy_objective_ = @(tax) dktax.policy_objective(mp_, sol2.p, tax, outcomes0, policyvar, policy_outcome);
    x0 = fzero(policy_objective_, bounds0, opts);
    mp_.tax_fuel = x0; 
    mp_=trmodel.update_mp(mp_);
    
    sol_ = equilibrium.solve(mp_, s, p1);
    outcomes_ = stats.compute_outcomes(mp_, s, sol_);

    outs{ir} = outcomes_; 
    tt(ir) = mp_.tax_fuel; 
end
close(h); 

%% 
f = graphs.myfigure(); 
plot(rates, tt, '-o', 'linewidth', 2); 
xlabel('Registration tax (relative to baseline)');
ylabel('Fuel tax (relative to baseline)', 'interpreter', 'none'); 
graphs.set_fig_layout_post();

if ~isempty(figpath)
    name_ = sprintf('%s/revenue_level_curve_tax_rates.eps', figpath);
    saveas(gcf, name_, 'epsc');
    fprintf('Figure saved as <a href="%s">%s</a>\n', figpath, name_);
end

%% Plot outcomes over the car tax rate 

% select outcomes to be plotted 
vars = {'social_surplus_total', 'social_surplus_ex_co2', 'total_co2', 'total_revenue', 'consumer_surplus'}; 

% convert from array of structs to matrix-form 
vv = nan(numel(rates), numel(vars)); 
for j=1:numel(vars) 
    v = vars{j};
    for ir=1:numel(rates)
        vv(ir,j) = outs{ir}.(v); 
    end
end

% plot a separate graph for each 
for j=1:numel(vars)
    f=graphs.myfigure(); 
    plot(rates, vv(:,j), '-o', 'linewidth', 2); 
    xlabel('Registration tax (relative to baseline)'); ylabel(stats.get_outcome_name(vars{j}), 'interpreter', 'none'); 
    graphs.set_fig_layout_post(f); 
    
    % save graph
    if ~isempty(figpath)
        name_ = sprintf('%s/revenue_level_curve_%s.eps', figpath, vars{j});
        saveas(gcf, name_, 'epsc');
        fprintf('Figure saved as <a href="%s">%s</a>\n', figpath, name_);
    end
end

end % for IMPERFECTPASSTHROUGH = 0:1 

fprintf('run_iruc_CF: total runtime = %5.2f min\n', toc(this_t)/60); 