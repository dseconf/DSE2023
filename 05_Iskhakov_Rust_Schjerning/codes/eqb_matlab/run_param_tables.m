% run_param_tables

%% Print parameter estimates 
clear all; close all; 
addpath('matlabinclude');
addpath('autotrade');

% set up parameters and indexes consistent with loaded data
load('results/estimation/estimation_data.mat');

% load parameter estimates
load('results/estimation/mle_converged.mat');


DOSAVE = 1;
DOLATEX = 1; 

%% Parameter names 

pn = struct(); % Parameter Names 
pn.acc_0 = 'Intercept'; 
pn.acc_a = 'Age slope'; 
pn.sigma_s = '$\sigma_s$: Scrap utility error variance'; 
pn.tc_sale = 'Intercept: selling (baseline is scrapping)'; 
pn.tc_sale_even = 'Selling in inspection years'; 

pn.mum = '$\mu_\tau$: marginal utility of money'; 

pn.psych_transcost = 'Intercept'; 
pn.psych_transcost_nocar = 'No car'; 

pn.u_0 = '$u_{\\tau,j,0} : $ intercept in indirect utility for car ownership';
pn.u_a = '$u_{\\tau,j,1} : $ coefficient on age in indirect utility for car ownership';


mp=mp_mle;

%% summary stats for cars
if DOSAVE 
    nam = 'results/tables/sumstats_cars.tex';
else
    nam = '';
end
data.sumstat_table_cars(dta, car, mp, nam);

%% Print parameter esimates
% 1st stage driving estimates 
if DOSAVE 
    nam = 'results/tables/tab_est_driving.tex';
else
    nam = '';
end
dktax.print_driving_parameters(mp,est,DOLATEX,nam);

se_mle=getfields(estim.update_mp(sqrt(diag(Avar_mle)), mp_mle), mp.pnames); 

% ntypes * ncartypes 
vv = {'u_0', 'u_a', 'u_a_sq'}; 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [mp.ntypes, mp.ncartypes]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    if ~DOSAVE
        fprintf('\n--- %s ---\n\n', v); 
    end
    tab = cell2table(mp_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', mp_mle.lbl_types); 
    tab_se = cell2table(se_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', mp_mle.lbl_types); 
    dktax.table2latex_se(tab, tab_se, fname(v, DOSAVE), pn.(v)); 
end


% ncartypes
vv = {'acc_0','acc_a'};
tab = cell(1, numel(vv)); 
tab_se = cell(1, numel(vv)); 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [1, mp.ncartypes]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    
    tab{i} = cell2table(mp_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', {pn.(v)}); 
	tab_se{i} = cell2table(se_mle.(v), 'VariableNames', mp_mle.lbl_cartypes, 'RowNames', {pn.(v)});     
end
if ~DOSAVE
    fprintf('\n--- %s ---\n\n', 'Accidents'); 
end
dktax.table2latex_se(vertcat(tab{:}), vertcat(tab_se{:}), fname('acc', DOSAVE)); 



% ntypes 
vv = {'mum'};
tab = cell(numel(vv),1); 
tab_se = cell(numel(vv),1); 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [mp.ntypes, 1]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    if ~DOSAVE
        fprintf('\n--- %s ---\n\n', v); 
    end
    tab{i} = cell2table(mp_mle.(v), 'VariableNames', {pn.(v)}, 'RowNames', mp_mle.lbl_types); 
    tab_se{i} = cell2table(se_mle.(v), 'VariableNames', {pn.(v)}, 'RowNames', mp_mle.lbl_types); 
end
dktax.table2latex_se([tab{:}], [tab_se{:}], fname('mum', DOSAVE)); 

% ntypes 
vv = {'psych_transcost', 'psych_transcost_nocar'};
title_str = 'Utility cost of transacting'; 
tab = cell(numel(vv),1); 
tab_se = cell(numel(vv),1); 
for i=1:numel(vv)
    v = vv{i}; 
    assert(all(size(mp_mle.(v)) == [mp.ntypes, 1]), 'Unexpected dimensions for %s', v);
    if all(cell2mat(mp_mle.(v)) == 0, 'all')
        fprintf('Skipping %s: all zeros\n', v); 
        continue
    end
    
    tab{i} = cell2table(mp_mle.(v), 'VariableNames', {pn.(v)}, 'RowNames', mp_mle.lbl_types); 
    tab_se{i} = cell2table(se_mle.(v), 'VariableNames', {pn.(v)}, 'RowNames', mp_mle.lbl_types); 
end
if ~DOSAVE
    fprintf('\n--- %s ---\n\n', 'Psychological transactions costs'); 
end
dktax.table2latex_se([tab{:}], [tab_se{:}], fname('psych_transcost', DOSAVE), title_str); 


% scrappage-related
if ~DOSAVE
    fprintf('\n--- scrappage ---\n\n'); 
end
vv = {'sigma_s', 'tc_sale', 'tc_sale_even'}; 
vv_lab = cell(size(vv)); 
x = []; 
x_se = []; 
for i=1:numel(vv) 
    if ismember(vv{i},mp_mle.pnames)
        vv_lab{i} = pn.(vv{i}); 
	    x_ = mp_mle.(vv{i}); 
	    x_se_ = se_mle.(vv{i}); 
	    if ~isscalar(x_)
	        warning('Parameter %s appears to be non-scalar: move it elsewhere', vv{i}); 
	        x_ = nan; 
	    end
	    x = [x; x_];
	    x_se = [x_se; x_se_];
	end
end
tab = array2table(x, 'rownames', vv_lab, 'variablenames', {'Estimate'}); 
tab_se = array2table(x_se, 'rownames', vv_lab, 'variablenames', {'Estimate'}); 
dktax.table2latex_se(tab, tab_se, fname('scrap', DOSAVE)); 

%% Helper function
function filename = fname(tabnam, DOSAVE)
    if DOSAVE 
        filename = sprintf('./results/tables/tab_mle_%s.tex', tabnam);
    else 
        filename = ''; 
    end
end
