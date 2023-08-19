% This script creates the figures used in Section 4

addpath('matlabinclude');
addpath('autotrade');
addpath('..');

clear
clc

%% 2: Persistently heterogeneous consumer economy
close all;

% common
mp = setparams.default(); % parameters used for illustration
mp.ntypes = 2;
mp.ncartypes=1; % switch to 1 car type
mp.lbl_cartypes = {' '}; % no label for the only car
s = trmodel.index(mp);

% 1. Low transaction cost
% 1.a set parameters
mp.transcost = 0;
mp = trmodel.update_mp(mp); % update parameters

% 1.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline

% 1.c graphs
f = graphs.agg_holdings(mp,s,sol);
ylim = get(gca, 'ylim'); 
legend('off');
saveas(f, 'results/illustration/example_holdings_tc_0.eps', 'epsc');
f1 = graphs.prices_with_hom(mp,s,sol);
legend('off');
saveas(f1, 'results/illustration/example_prices_tc_0.eps', 'epsc');

% 2. With high transaction costs

% 2.a set parameters
mp.transcost=10;
mp.psych_transcost = {0}; 
mp=trmodel.update_mp(mp); % standard model parameters

% 2.b solve model
[sol]=equilibrium.solve(mp, s); % solve model in baseline

% 2.c graphs
f = graphs.agg_holdings(mp,s,sol);
set(gca, 'ylim', ylim);
legend('location', 'north'); 
saveas(f, 'results/illustration/example_holdings_tc_10.eps', 'epsc');
f2 = graphs.prices_with_hom(mp,s,sol);
saveas(f2, 'results/illustration/example_prices_tc_10.eps', 'epsc');



