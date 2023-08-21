%% Print parameter estimates 
clear all; close all; 

%datadir='../../../../data/8x4/'; 
addpath('matlabinclude');
addpath('autotrade');

% set up parameters and indexes consistent with loaded data
load('results/estimation/estimation_data.mat');

% load parameter estimates
load('results/estimation/mle_converged.mat');


%% Print parameter esimates

mp=mp_mle;

%% Market shares
% ----------------------------------------------------
mp=mp_mle;
mp.lbl_cartypes={'LB', 'LG', 'HB', 'HG'};
set(groot, 'defaultaxesfontsize', graphs.fontsize);
graphs.fit_market_shares(mp, s, sol_mle, sol_d);
set(gcf, 'units','normalized','outerposition',[0 0 1 1])
saveas(gcf, 'results/model_fit/fit_marketshares.eps', 'epsc'); 

graphs.fit_market_shares_agg(mp, s, sol_mle, sol_d);
axis square                                                        
title('')
saveas(gcf, 'results/model_fit/fit_marketshares_agg.eps', 'epsc'); 

graphs.fit_market_shares_45(mp, s, sol_mle, sol_d);
saveas(gcf, 'results/model_fit/fit_marketshares_45.eps', 'epsc'); 

%% Purchases
% ----------------------------------------------------
%graphs.fit_purchases(mp_mle,s, sol_mle, sol_d, 'by car type');
set(groot, 'defaultaxesfontsize', graphs.fontsize); 
graphs.fit_purchases(mp_mle,s, sol_mle, sol_d, 'all cars'); 
saveas(gcf, 'results/model_fit/fit_purchases.eps', 'epsc'); 

%% Keep
% ----------------------------------------------------
%graphs.fit_keep(mp_mle,s, sol_mle, sol_d, 'all cars'); 
set(groot, 'defaultaxesfontsize', graphs.fontsize); 
graphs.fit_keep(mp_mle,s, sol_mle, sol_d, 'by car type');
set(gcf, 'units','normalized','outerposition',[0 0 .6 .8])
saveas(gcf, 'results/model_fit/fit_keep.eps', 'epsc'); 


%% holdings - normalized
% ----------------------------------------------------
mp=mp_mle;
mp.lbl_types={'LWD, C, Poor' ,'LWD, C, Rich' ,'LWD, S, Poor' ,'LWD, S, Rich' ,...
               'HWD, C, Poor','HWD, C, Rich','HWD, S, Poor','HWD, S, Rich'};

set(groot, 'defaultaxesfontsize', graphs.fontsize);
graphs.outcomes(mp, s, {'holdings'}, sol_mle, sol_d, {'Data'})
set(gcf, 'units','normalized','outerposition',[0 0 1 1])
saveas(gcf, 'results/model_fit/holdings_data.eps', 'epsc'); 
graphs.outcomes(mp, s, {'holdings'}, sol_mle, sol_d, {'Model'})
set(gcf, 'units','normalized','outerposition',[0 0 1 1])
saveas(gcf, 'results/model_fit/holdings_model.eps', 'epsc'); 

%% post trade distribution (quantity post trade)
% ----------------------------------------------------
mp=mp_mle;
mp.lbl_cartypes={'LB', 'LG', 'HB', 'HG'};
mp.lbl_types={'LWD, C, Poor' ,'LWD, C, Rich' ,'LWD, S, Poor' ,'LWD, S, Rich' ,...
               'HWD, C, Poor','HWD, C, Rich','HWD, S, Poor','HWD, S, Rich'};

set(groot, 'defaultaxesfontsize', graphs.fontsize);
graphs.outcomes(mp, s, {'post_trade_dist'}, sol_mle, sol_d, {'Data'})
ylim([0,1.1]);
set(gcf, 'units','normalized','outerposition',[0 0 .5 .5])
saveas(gcf, 'results/model_fit/pt_dist_data.eps', 'epsc');
set(groot, 'defaultaxesfontsize', graphs.fontsize); 
graphs.outcomes(mp, s, {'post_trade_dist'}, sol_mle, sol_d, {'Model'})
ylim([0,1.1]);
set(gcf, 'units','normalized','outerposition',[0 0 .5 .5])
saveas(gcf, 'results/model_fit/pt_dist_model.eps', 'epsc'); 


%% keep shares - by household
% ----------------------------------------------------
set(groot, 'defaultaxesfontsize', graphs.fontsize);
graphs.outcomes(mp, s, {'keep'}, sol_mle, sol_d, {'Data'})
set(gcf, 'units','normalized','outerposition',[0 0 1 1])
saveas(gcf, 'results/model_fit/keep_data_tau.eps', 'epsc');
set(groot, 'defaultaxesfontsize', graphs.fontsize); 
graphs.outcomes(mp, s, {'keep'}, sol_mle, sol_d, {'Model'})
set(gcf, 'units','normalized','outerposition',[0 0 1 1])
saveas(gcf, 'results/model_fit/keep_model_tau.eps', 'epsc'); 

%% scrap probability
% ----------------------------------------------------
mp=mp_mle;
set(groot, 'defaultaxesfontsize', graphs.fontsize);
% graphs.ccp_scrap_compare(mp, s, sol_d, sol_mle);
graphs.fit_scrap(mp_mle,s, sol_mle, sol_d, 'by car type'); 
set(gcf, 'units','normalized','outerposition',[0 0 .4 .6])
saveas(gcf, 'results/model_fit/fit_scrap.eps', 'epsc'); 

%% prices
% ----------------------------------------------------
set(groot, 'defaultaxesfontsize', graphs.fontsize);
graphs.outcomes(mp, s, {'prices'}, sol_mle, sol_d, {'Model'})
set(gcf, 'units','normalized','outerposition',[0 0 .3 .6])
set(gca, 'fontname', graphs.fontname, 'Ygrid', 'on', 'box', 'off'); 
ylim([0, 300])
saveas(gcf, 'results/model_fit/prices_model.eps', 'epsc'); 


