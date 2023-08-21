classdef dktax
    properties (Constant)
        scaling=1;
        params={'vat','cartax_lo','cartax_hi','tax_fuel','K_cartax_hi','dk_population','fe','p_fuel'};
    end
    
    methods (Static)
        
        function diff = policy_objective(mp0, car_prices, tax, outcomes0, tax_variable_name, target_outcome_name)
            % policy_objective: the difference in the outcome,
            % "outcomes0.(outname)", and the outcome in the counterfactual
            % implied by setting the value "tax". The function is intended
            % to be used with a bisection search. 
            %
            % INPUTS: 
            %   mp0: baseline parameter values 
            %   car_prices: vector of car prices at which to solve the
            %   model 
            %   tax: (scalar) value of the tax lever 
            %   outcomes0: (struct) outcomes computed in the baseline.  
            %   tax_variable)name: (cell array of length 1) Name of the policy
            %    variable (must be a field in mp0.)
            %   target_outcome_name: (string) Name of the outcome 
            %    (must be a field in outcomes0.)
            %
            % OUTPUT: 
            %   diff: (scalar) target1 - target0 (NB! not in absolute
            %    value) 
            
            assert(ischar(tax_variable_name) && isfield(mp0, tax_variable_name), 'tax_variable_name must be a char specifying a field in mp0.'); 
            assert(ischar(target_outcome_name) && isfield(outcomes0, target_outcome_name), 'target_outcome_name, "%s" must be a string specifying a field in outcomes0.', target_outcome_name); 
            
            % update parameter vector, setting mp1.tax_variable{1} = tax. 
            mp1 = vec2struct(tax, {tax_variable_name}, mp0);
            mp1 = trmodel.update_mp(mp1);
            
            s = trmodel.index(mp1); 
            sol1 = equilibrium.solve(mp1, s, car_prices); 
            outcomes1 = stats.compute_outcomes(mp1, s, sol1); 

            % objective function value 
            diff=outcomes1.(target_outcome_name) - outcomes0.(target_outcome_name);
        end

        function print_driving_parameters(mp, est, DOLATEX, fname)

            if nargin>3 
                fid = fopen(fname,'w');
            else
                fid=1;
            end
                
            % 2021-07-16: print estimates from the first-stage driving
            % regression.
            
            % assign cell separator and line break 
            if DOLATEX
                sepa = ' & ';        
                endline = '\\\\ \n'; 
            else
                sepa = ' ';
                endline = '\n';
            end
            
            symbols = struct(); 
            names   = struct(); 
            
            v = matlab.lang.makeValidName('(Intercept)'); 
            symbols.(v) = '$\gamma_0$'; 
            names.(v) = 'Intercept'; 
            
            symbols.pt_car_age = '$\hat{\gamma}^a_1 / \phi_\tau$'; 
            names.pt_car_age = 'Car age'; 
            
            symbols.pt_car_age_sq = '$\hat{\gamma}^a_2 / \phi_\tau$'; 
            names.pt_car_age_sq = 'Car age squared';             
            
            symbols.ppk = '$\mu / \phi$'; 
            names.ppk   = 'Price (common)'; 
            for tau=1:mp.ntypes 
                v = sprintf('tau_cat_%d', tau); 
                v = matlab.lang.makeValidName(v); 
                symbols.(v) = '$\hat{\gamma}_{\tau} / \phi_\tau$'; 
                names.(v) = sprintf('Intercept, %s', mp.lbl_types{tau}); 
                
                v = sprintf('ppk:tau_cat_%d', tau); 
                v = matlab.lang.makeValidName(v); 
                symbols.(v) = '$\hat{\mu}_{\tau} / \phi_\tau$'; 
                names.(v) = sprintf('Price, %s', mp.lbl_types{tau}); 
            end
            
            for j=1:mp.ncartypes
                v = sprintf('pt_car_type_cat_%d', j); 
                v = matlab.lang.makeValidName(v); 
                symbols.(v) = '$\hat{\gamma}_j / \phi_\tau$';
                names.(v) = sprintf('Car dummy: %s', mp.lbl_cartypes{j}); 
            end
            
            if DOLATEX
                fprintf(fid, '\\begin{tabular}{llrr} \n');
            end
            
            print_str = sprintf('%%-30s %5s %%12.4g %5s (%%1.2f) %s', sepa, sepa, endline);
            print_str2 = '%40s %s %40s';
            print_header = @(txt) fprintf(fid, sprintf('\\\\toprule %s{%s} %s\\\\midrule \n', '\\multicolumn{4}{c}', txt, endline));
                        
            print_header('Dependent variable: thousands of kilometers driven per year');
            
            print_row = @(symbol, name, estimate, stderr) fprintf(fid, print_str, sprintf(print_str2, symbol, sepa, name), estimate, stderr);
                       
            vv_ = est.Properties.RowNames; 
            
            % change order of ppk and car age squared 
            vv = vv_; 
            vv(3) = vv_(15); 
            vv(15) = vv_(3); 
            vv(11) = []; % delete the reference car category 
            for i=1:numel(vv)
                v = vv{i};
                v_ = matlab.lang.makeValidName(v); 
                assert(isfield(symbols, v_), 'Field "%s% not added to manual list above.', v);                     
                print_row(symbols.(v_), names.(v_), est(v,:).Estimate, est(v,:).SE); 
            end

            % no. observations 
            nobs_str = '19,635,940';
            if DOLATEX
                fprintf(fid, '\\midrule \n');  
                fprintf(fid, 'N & Driving periods & \\multicolumn{2}{r}{%s}  \\\\ \n', nobs_str)
            else
                frpintf(fid, 'N (driving periods): %s', nobs_str)
            end
            
            % end of table 
            if DOLATEX
                fprintf(fid, '\\bottomrule\n\\end{tabular} \n');
            end

            if nargin>3 
                fclose(fid);
            end
        end % print_driving_parameters 
        
        function fid = get_fid(fname)
            if isempty(fname)
                fid = 1;
            else
                fid = fopen(fname, 'w'); 
            end
        end

        function close_fid(fid)
            if fid~=1
                fclose(fid);
            end
        end

        function table2latex(tab, fname)
            
            fid = dktax.get_fid(fname);
            
            assert(not(isempty(tab.Properties.RowNames)), 'Only implemented with RowNames for the table');
            
            sepa = ' & ';
            endline = '\\\\ \n';
            
            vars = tab.Properties.VariableNames; 
            rows = tab.Properties.RowNames; 
            
            nvars = numel(vars); 
            ncol = nvars + 1; 
            
            fprintf(fid, ['\\begin{tabular}{l' repmat('r', 1, numel(vars)) '} \n']);
            
            table_row_str  = ['%30s ' repmat([sepa '%12.4f  '], 1, nvars) endline]; % %30s & %12.4f & %12.4f ... %12.ff \\
            table_head_str = ['%30s ' repmat([sepa '%12s  '  ], 1, nvars) endline]; 
            
            fprintf(fid, table_head_str, '', vars{:}); 
            
            for i=1:numel(rows) 
                row = rows{i}; 
                fprintf(fid, table_row_str, row, tab(i,:).Variables); 
            end
            
            fprintf(fid, '\\bottomrule \\end{tabular} \n');

            if fid~=1
                fprintf('File written: %s\n', fname);
                dktax.close_fid(fid);
            end            
            
        end

        function table2latex_se(tab,tab_se,fname,title_str)
            assert(not(isempty(tab.Properties.RowNames)), 'Only implemented with RowNames for the table');
            
            if nargin < 4 || isempty(title_str) 
                DOTITLE = false; 
            else
                DOTITLE = true; 
            end
            
            fid = dktax.get_fid(fname);

            sepa = ' & ';
            endline = '\\\\ \n';
            
            vars = tab.Properties.VariableNames; 
            rows = tab.Properties.RowNames; 
            
            nvars = numel(vars); 
            ncol = nvars + 1; 
            
            fprintf(fid, ['\\begin{tabular}{l' repmat('r', 1, numel(vars)) '} \n']);
            
            table_row_str  = ['%30s ' repmat([sepa '%12.4f  '], 1, nvars) endline]; % %30s & %12.4f & %12.4f ... %12.ff \\
            table_head_str = ['%30s ' repmat([sepa '%12s  '  ], 1, nvars) endline]; 
            table_row_str_se  = ['%30s ' repmat([sepa '(%1.4f)  '], 1, nvars) endline]; % %30s & %12.4f & %12.4f ... %12.ff \\
            
            fprintf(fid, '\\toprule \n'); 
            
            if DOTITLE
                mystr = ['\\multicolumn{' num2str(numel(vars)+1) '}{c}{' title_str '} \\\\ \n\\midrule\n'];
                fprintf(fid, mystr);
            end
            
            fprintf(fid, table_head_str, '', vars{:}); 
            
            fprintf(fid, '\\midrule \n'); 
            
            for i=1:numel(rows) 
                row = rows{i}; 
                fprintf(fid, table_row_str, row, tab(i,:).Variables); 
                fprintf(fid, table_row_str_se, '', tab_se(i,:).Variables); 
            end
            
            fprintf(fid, '\\bottomrule \\end{tabular} \n');

            if fid~=1
                fprintf('File written: %s\n', fname);
                dktax.close_fid(fid);
            end        
        end
        
        % FIXME: update to new model setup 
        function print_model_parameters(mp, DOLATEX, fid)
            db = mp.db; % for short refencing below

            % assign cell separator and line break 
            if DOLATEX
                sepa = ' & ';        
                endline = '\\\\ \n'; 
            else
                sepa = ' ';
                endline = '\n';
            end
            
            if DOLATEX
                fprintf(fid, '\\begin{tabular}{llr} \n');
            end
            
            print_str = sprintf('%%-30s %s %%12.4g %s', sepa, endline);
            print_str2 = '%30s %s %25s';
            print_header = @(txt) fprintf(fid, sprintf('\\\\midrule %s{%s} %s\\\\midrule \n', '\\multicolumn{3}{c}', txt, endline));
            
            types = {'rich', 'poor'};
            cars = {'normal', 'luxury'};
            
            print_header('First-stage: linear driving equation');
            
            % intercept
            for tp=1:2
                type = types{tp};
                var = sprintf('%s %s', 'Intercept', type);
                fprintf(fid, print_str, sprintf(print_str2, '$\hat{\gamma}_\tau / \phi$', sepa, var), db.tau{tp});
            end
            % price per km
            for tp=1:2
                type = types{tp};
                var = sprintf('%s %s', 'Price per km, ', type);
                fprintf(fid, print_str, sprintf(print_str2, '$\hat{\mu}_\tau / \phi$', sepa, var), db.pkm{tp});
            end
            % car type dummies
            for car=1:2
                this_car = cars{car};
                var = sprintf('%s %s', 'Car type dummy, ', this_car);
                fprintf(fid, print_str, sprintf(print_str2, '$\hat{\gamma}_j / \phi$', sepa, var), db.car{car});
            end
            % car age
            fprintf(fid, print_str, sprintf(print_str2, '$\hat{\gamma}^{a}_1 / \phi$', sepa, 'Car age'), db.a1{1});
            fprintf(fid, print_str, sprintf(print_str2, '$\hat{\gamma}^{a}_2 / \phi$', sepa, 'Car age squared'), db.a2{1});
            
            assert(db.a1{1} == db.a1{2}, 'Expected identical coefficients across HH types for the age polynomial');
            
            
            print_header('Second-stage: Method of moments');
            
            f_ = @(tex, plain, value) fprintf(fid, print_str, sprintf(print_str2, tex, sepa, plain), value);
            
            % intercept in utility
            for car=1:2
                this_car = cars{car};
                f_('$\hat{\delta}^0_j$', sprintf('Intercept, %s', this_car), mp.u_0{1,car});
            end
            assert(mp.u_0{1,1} == mp.u_0{2,1}, 'Expecting symmetric utility intercepts over consumers');
            
            % age profile
            for car=1:2
                this_car = cars{car};
                f_('$\hat{\delta}^a_j$', sprintf('Car age, %s', this_car), mp.u_a{1,car});
            end
            assert(mp.u_a{1,1} == mp.u_a{2,1}, 'Expecting symmetric utility car age profiles over consumers');
            
            % mu (marg. util of money)
            for tp=1:2
                phi = 1/( -db.pkm{tp}/mp.mum{tp}  );
                type = types{tp};
                %f_('$\phi_\tau$', sprintf('Driving curvature, %s', type), phi);
                f_('$\hat{\mu}_\tau$', sprintf('Utility of money, %s', type), mp.mum{tp});
            end
            
            f_('$\alpha$', 'Accident probability', mp.acc_0{1});
            
            print_header('Fixed parameters');
            
            f_('$\beta$', 'Discount factor', mp.bet);
            
            f_('', 'Transaction cost', mp.transcost);
            
            % end
            if DOLATEX
                fprintf(fid, '\\bottomrule\n\\end{tabular} \n');
            end
        end

        function print_outcomes(outcomes, mp, headerStr, fname, DOLATEX)
            % prints details from a simulation from the model
            
            if nargin < 6
                DOSHORT = 1; 
            end
            if nargin < 5 || isempty(DOLATEX) 
                DOLATEX = 0; 
            end
            if nargin<4 
                fname = ''; % print to screen
            end
            if nargin<3 || isempty(headerStr)
                headerStr = 'OUTCOMES';
            end
            
            dktax.print_outcomes_comparison({outcomes}, {mp}, headerStr, {headerStr}, fname, DOLATEX, DOSHORT);
        end
        
        function print_outcomes_comparison(outcomes_cell, mp_cell, headerStr, names_cell, outputFileName, DOLATEX, DOSHORT)
            
            if ~isempty(outputFileName)
                if DOLATEX
                    fExt = 'tex';
                else
                    fExt = 'txt';
                end
                
                if DOSHORT
                    outputFileName=sprintf('%s_short.%s', outputFileName,fExt);
                else
                    outputFileName=sprintf('%s.%s', outputFileName,fExt);
                end
                
                fid = fopen(outputFileName, 'w+');
                fprintf('Writing comparison table to <a href="%s">%s</a>.\n', outputFileName, outputFileName);
            else
                fid = 1;
            end
            
            if DOLATEX
                fExt = 'tex';
                sep = '&';
                endLine = '\\\\ \n';
                print_separator = @() fprintf(fid, '\\midrule \n');
                print_separator_top = @() fprintf(fid, '\\toprule \n');
                print_separator_bot = @() fprintf(fid, '\\bottomrule \n');
                
                print_sec_header = @(s) fprintf(fid, '\\multicolumn{5}{c}{\\emph{%s}} \\\\ \n\\midrule\n', s);
            else
                fExt = 'txt';
                sep = ' ';
                endLine = '\n';
                print_separator = @() fprintf(fid, '%s\n', repmat('-',1,90));
                print_separator_top = print_separator; 
                print_separator_bot = print_separator; 
                
                print_sec_header = @(s) fprintf(fid, '   --- %s ---\n', s);
            end
            
            % print imperfect passthrough information? 
            DOIMPERFECTPASSTHROUGH = false;
            for i=1:numel(mp_cell)
                mp = mp_cell{i}; 
                if isfield(mp, 'passthrough')
                    DOIMPERFECTPASSTHROUGH = true;
                end
            end
            
            N = numel(outcomes_cell);
            assert(numel(mp_cell) == N, 'Dimensions mismatch');
            
            % 1. header
            c=clock; d=date;
            
            if DOLATEX
                tabType = sprintf('l%s', repmat('r', 1, N));
                fprintf(fid, '\n\n\\begin{tabular}{%s} \n', tabType);
            else
                print_separator();
                fprintf(fid,'%s %1.0f:%02.0f:%02.0f\n%s\n',d,c(4:6));
                fprintf(fid,'%s\n', headerStr);
            end
            
            print_separator_top();
            
            % compute fuel tax in pct. of total price at the pump
            for j=1:numel(mp_cell)
                mp_cell{j}.p_fuel_tax_share = (mp_cell{j}.p_fuel - mp_cell{j}.p_fuel_notax)/mp_cell{j}.p_fuel;
                mp_cell{j}.p_fuel_dkk = mp_cell{j}.p_fuel * 1000; 
            end
            
            % 2. build print string
            % each row has structure "VARNAME = VAL VAL VAL ... \n"
            % VARNAME is a 25 long string, while each VAL has 12 spaces and
            % requires 3 significant digits.
            
            % 2.a content lines: "VARNAME = 30.1 4.2 412.3 ... \n"
            numStr = repmat(sprintf('%s %%12.3f', sep), 1, N);
            printStr = sprintf('%%40s %s %s', numStr, endLine);
            
            % 2.b header
            strStr = repmat(sprintf('%s %%12s', sep), 1, N);
            headerStr = sprintf('%40s %s %s', ' ', strStr, endLine); % empty string: nothing to write in the first
            
            % 2.c print header
            fprintf(fid, headerStr, names_cell{:});
            print_separator();
            
            % 3. print contents
            % toPrint is the name of a member of either mp or outcomes
            % labels has the corresponding human readable name of that
            % variable.
            
            toPrintPolicy = { 'cartax_lo',     'cartax_hi',     'p_fuel_tax_share'                }; 
            labelsPolicy =  { 'Registration tax (bottom rate)', 'Registration tax (top rate)', 'Fuel tax (share of pump price)' };  
            
            labs = { ... 
                ... 'variable name',     'readable label'; ... 
                'social_surplus_total',  'Social surplus (1000 DKK)'; ... 
                'total_revenue',         'Total tax revenue (1000 DKK)'; ... 
                'revenue_fuel',          'Fuel tax revenue (1000 DKK)'; ... 
                'revenue_car',           'Car tax revenue (1000 DKK)'; ... 
                'extern_vkt',            'Non-CO2 externalities (1000 DKK)'; ... 
                'total_extern',          'Externalities (1000 DKK)'; ... 
                'consumer_surplus',      'Consumer surplus (1000 DKK)'; ... 
                'total_co2',             'CO2 (ton)'; ... 
                'total_vkt',             'VKT (1000 km)'; ... 
                'Ecarage',               'E(car age)'; ... 
                'nocar_share',           'Pr(no car)' ... 
            }; 
        
            if DOIMPERFECTPASSTHROUGH
                labs(1, :) = {'social_surplus_total_incl_profits', 'Social surplus (1000 DKK)'}; 
                labs = [labs(1:end-4, :); ...
                       {'markup_car', 'Producer surplus (change, 1000 DKK)'}; ...
                        labs(end-3:end,:)];
            end
            
        
            toPrint = labs(:, 1); % varnames 
            labels  = labs(:, 2); % human readable labels
            
            % deal with "|"
            if DOLATEX
                for i=1:numel(labels)
                    labels{i} = strrep(labels{i}, '|', '$\vert$');
                end
            end
            
            % 3.a Print policy variables 
            
            print_sec_header('Policy choice variables');
            
            for i=1:numel(toPrintPolicy)
                % 3.a.i construct cell-vector of values
                x = toPrintPolicy{i};
                vec = cell(N,1);
                for j=1:N
                    % must find out whether x is in outcomes or mp
                    if isfield(outcomes_cell{j}, x)
                        vec{j} = outcomes_cell{j}.(x);
                    elseif isfield(mp_cell{j}, x)
                        vec{j} = mp_cell{j}.(x);
                    else
                        error('Cannot find field %s in either mp or outcomes.', x);
                    end
                end
                
                % 3.a.ii print
                fprintf(fid, printStr, labelsPolicy{i}, vec{:});
            end
            
            print_separator(); 
            
            % 3.b implications for exogenous prices 
            
            print_sec_header('Exogeneous prices');
            
            % 3.b.i copy into convenient cell 
            mp = mp_cell{1}; % for convenience: we will only be using fields that should not be varying 
            pnew_vec = cell(N,mp.ncartypes);
            p_fuel = cell(N, 1);
            for i=1:N
                p_fuel{i} = mp_cell{i}.p_fuel_dkk; 
                for j=1:mp.ncartypes
                    pnew_vec{i, j} = mp_cell{i}.pnew{j};
                end
            end

            % 3.b.ii print
            for j=1:mp.ncartypes
                s_ = sprintf('Price, %s (1000 DKK)', mp.lbl_cartypes{j}); 
                fprintf(fid, printStr, s_, pnew_vec{:, j});
            end
            
            % 3.b.iii fuel prices 
            fprintf(fid, printStr, 'Fuel price (DKK/l)', p_fuel{:});
            
            print_separator(); 
            
            print_sec_header('Outcomes');
            
            % 3.c outcome variables 
            for i=1:numel(toPrint)
                
                % 3.b construct cell-vector of values
                x = toPrint{i};
                vec = cell(N,1);
                for j=1:N
                    % must find out whether x is in outcomes or mp
                    if isfield(outcomes_cell{j}, x)
                        vec{j} = outcomes_cell{j}.(x);
                    elseif isfield(mp_cell{j}, x)
                        vec{j} = mp_cell{j}.(x);
                    else
                        error('Cannot find field %s in either mp or outcomes.', x);
                    end
                end
                
                % 3.c print
                fprintf(fid, printStr, labels{i}, vec{:});
                
            end
            
            print_separator_bot();
            
            % 4. end
            if DOLATEX
                fprintf(fid, '\\end{tabular}\n');
            end
            
            if ~isempty(outputFileName)
                fclose(fid);
            end
        end
        
    end % end of methods
end % end of class

