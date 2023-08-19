classdef toy
% Toy model code written by Anders Munk-Nielsen 

    methods (Static)
        function cat_within = split_on_x_var_prcts_within_category(x, category, prcts, sales)
            % INPUTS:
            %   x: J-vector of continuous values (does not make sense for integer
            %   variable
            %   category: vector of integer values.
            %   prcts: percentile lists. E.g. prcts=50 to split at above/below median
            %   sales: [optional] J-vector of counts by which to weight percnetiles 
            % OUTPUTS:
            %   cat_within: J-vector of intengers in {1,2,...,numel(prcts+1)}
            
            if nargin < 4 || isempty(sales)
                sales = ones(size(x)); % equally weighted
            end
            
            if prcts(1) ~= 0.0
                prcts = [0.0;prcts];
            end
            if prcts(end) ~= 100.0
                prcts = [prcts; 100.0];
            end
            
            cc = unique(category(not(isnan(category))));
            Nc = numel(cc);
            cat_within = nan(size(x));
            for i=1:Nc
                I = category == cc(i);
                x_ = repelem(x(I), sales(I)); % we want sales-weighted percentiles
                if all(isnan(x_))
                    warning 'Not a single unit with positive sales has observed value of characteristic x.'
                    keyboard; 
                end
                xx = prctile(x_, prcts);
                xx(end) = inf; % discretize will not include right-end-points, so the max should not be there
                cat_within(I) = discretize(x(I), xx);
            end
            
        end
        
        function x_imputed = impute_missings(x, carAll_t, regressors)
            if isempty(regressors)
                regressors = {'diesel', 'totalvaegt', 'vintage'};
            end
            
            assert(isfield(carAll_t, x), 'input x must be a field in carAll_t');
            for i=1:numel(regressors)
                X_ = regressors{i};
                assert(isfield(carAll_t, X_), 'input x must be a field in carAll_t');
            end
            
            X = nan(carAll_t.nCars, numel(regressors));
            for i=1:numel(regressors)
                X(:, i) = carAll_t.(regressors{i});
            end
            T = array2table([X, carAll_t.(x)], 'VariableNames', {regressors{:}, x});
            m = fitlm(T); % without a model specification, it assumes the *last* variable is Y
            x_imputed = m.predict(); % will inherit missings from missings in the regressors
        end
        
        function [var_num, types] = create_int_group(categorical_var)
            % INPUTS:
            %   categorical var: N*num_var, where num_var = the number of variables
            %   based on which we create unique combinations into a group
            % OUTPUTS: 
            %   var_num: N-vector of integers indicating the groups. Will
            %   be nan if any of the categorical var rows are missing 
            %   types: Ntypes*num_var matrix that maps from the types to
            %   the underlying categorical values 
            %
            var_num = nan(size(categorical_var,1), 1); 
            Iobs = not(any(isnan(categorical_var), 2)); 
            
            if size(categorical_var, 2) == 1 % vector
                [types,b_,var_num(Iobs)] = unique(categorical_var(Iobs));
            else % matrix of categories
                [types,b_,var_num(Iobs)] = unique(categorical_var(Iobs, :), 'rows');
            end
            
        end
        
        function id = create_decision_vector(dat, ntypes, nages) 
            % construct integer ID for discrete car choice (-1 means keep,
            % -2 means purge, 0 means new car of first type, etc.).
            % Construction relies on "dKeepCar", "dBuyCar", etc. 
            
            N = numel(dat.iecar); 
            
            % initialize 
            id = nan(N,1); 
            
            % buy a car 
            I = (dat.dBuyCar == 1) | (dat.dReplaceCar == 1); 
            assert(not(any(dat.icardrive(I) == -1)), 'Found obs. where a car was bought/replaced but icardrive was -1.');
            id(I) = toy.car_type_and_age_to_icar(dat.icardrive(I), dat.carageyrsout(I), ntypes, nages); 
            
            % check result 
            J = I & isnan(dat.carageyrsout); 
            if any(J)
                warning('Encountered %d obs (%5.2f%%) which replace car but where car age outgoing is not observed. Setting to purge', sum(J), mean(J)*100.0); 
                id(J) = -2; 
            end
            
            % keep choice
            I = dat.dKeepCar == 1; 
            id(I) = -1;%toy.car_type_and_age_to_icar(dat.icardrive(I), dat.carageyrsin(I), ntypes, nages); 
            
            % keep no car 
            I = dat.dKeepNoCar == 1; 
            id(I) = -1; 
            
            % purge 
            I = dat.dSellCar == 1; 
            id(I) = -2; 
            
        end
        
        function icar = car_type_and_age_to_icar(car_type, car_age, ntypes, nages)
            % INPUTS: 
            %   car_type: base 1 car type 
            %   car_age: age in years "base 0" so to speak (0 y/o is an
            %   option)
            % OUTPUT
            %   icar: base 1 car type variable
            
            N = numel(car_type); 
            assert(all(car_type ~= 0), 'assumes car_type base 1, found a zero.');
            nc = numel(unique(car_type(car_type ~= -1)));
            assert(nc <= ntypes, 'ntypes is %d but found %d types.', ntypes, nc);
            
            I = car_age >= nages;
            if any(I)
                warning('nages=%d, but found a car aged %d. setting all cars aged >= %d (%d obs) to %d.', nages, max(car_age), nages, sum(I), nages-1');
                car_age(I) = nages-1;
            end
            assert(not(any(car_age < 0)), 'negative car ages found')
            
            I = isnan(car_age); 
            assert(all(car_type(I) == -1))
            assert(all(car_type <= ntypes), 'Found car types above ntypes=%d', ntypes); 
            
            % initialize 
            icar = nan(N,1); 
            
            % no car
            I = car_type == -1; 
            icar(I) = -1; 
            
            % has car 
            I = car_type >= 1 & car_type <= ntypes; 
            icar(I) = (car_type(I)-1)*nages + car_age(I) + 1; % +1 to get to base 1 (recall that car age is base 0!
            
        end
        
        function group = create_int_group_invertible(iz,icar,id,nZ,nCar)
            %
            % INPUTS:
            %   iz: demographic group index, in [1, 2, ..., nZ]
            %   icar: in [1,...,nCar]
            
            % 0. checks
            assert(size(iz,2) == 1, 'iz should be a column vector'); 
            assert(min(iz) == 1, 'expecting iz to be base 1 but min(iz) is not 1')
            assert(all(iz == round(iz)), 'iz does not appear to be integer');
            assert(min(id) >= -2, 'id must have min >= -2');
            assert(min(icar) >= -1, 'icar must have min >= -1');
            assert(not(any(icar == 0)), 'icar = 0 detected, that choice should not occur');
            assert(not(any(id == 0)),   'id = 0 detected, that choice should not occur');
            
            if max(id) ~= max(icar)
                warning('id and icar have different max values... this could be a problem');
            end
            
            % 1. handle "-1", "-2" (keep, purge)
            
            % 1.a icar 
            I = icar == -1; 
            icar(I) = nCar+1; 
            
            % 1.b id
            I = id == -1;
            id(I) = nCar + 1;
            I = id == -2;
            id(I) = nCar + 2;
            
            % 2. create group variable
            group = sub2ind([nZ, nCar+1, nCar+2], iz, icar, id);
            
        end
        
        function [iz, icar, id] = create_group_vars_from_int(group, nZ, nCar)
            % 1. convert to index form
            [iz, icar, id] = ind2sub([nZ, nCar+1, nCar+2], group);
            
            % 2. handle keep/purge
            
            % 2.a icar
            I = icar == nCar+1; 
            icar(I) = -1; 
            
            % 2.b id 
            I = id == nCar+1;
            id(I) = -1;
            I = id == nCar+2;
            id(I) = -2;
        end
        
        function C = extract_car_attributes(dataAll, carAll, typeVar, USESALESWEIGHTEDAVGATTRIBUTES, numCarAges)
            % converts carAll (which has one row per unique car) to C
            % (which has one row per car type). Relies on the function
            % "get_car_attributes()", which should be in the Matlab path. 
            
            numT = numel(dataAll); 
            assert(numel(carAll) == numT); 
            assert(isfield(carAll{1}, typeVar), 'typeVar, "%s% not found in carAll{1} struct', typeVar);
            
            C = cell(numT, 1);
            for t=1:numT
                % 1. delete observations with illegal cars (type -99)
                
                % 1.a iecar
                this_count = 0;
                for i=1:dataAll{t}.NT
                    this_iecar = dataAll{t}.iecar(i);
                    switch this_iecar
                        case -1
                            % fine skip it
                        case 0 % should not happen
                            dataAll{t}.iecar(i) = -1; this_count = this_count + 1;
                            this_count = this_count + 1;
                        otherwise
                            type = carAll{t}.(typeVar)(dataAll{t}.iecar(i));
                            if type == -99
                                dataAll{t}.iecar(i) = -1;
                                this_count = this_count + 1;
                            end
                    end
                end
                fprintf('changing type for %d impossible cases for iecar, ', this_count);
                
                % 1.b icardrive 
                this_count = 0;
                for i=1:dataAll{t}.NT
                    this_icardrive= dataAll{t}.icardrive(i);
                    switch this_icardrive
                        case -1
                            % fine skip it
                        case 0 % should not happen
                            dataAll{t}.icardrive(i) = -1; this_count = this_count + 1;
                            this_count = this_count + 1;
                        otherwise
                            type = carAll{t}.(typeVar)(dataAll{t}.icardrive(i));
                            if type == -99
                                dataAll{t}.icardrive(i) = -1;
                                this_count = this_count + 1;
                            end
                    end
                end
                fprintf('and %d for icardrive\n', this_count);
                
                % 2. get car attributes
                C{t} = get_car_attributes(carAll{t}, dataAll{t}, typeVar, USESALESWEIGHTEDAVGATTRIBUTES, numCarAges);
                
                % 3. replace missings with imputed values
                ff = fieldnames(C{t}.imputed_pooled);
                for i=1:numel(C{t}.imputed_pooled)
                    v = ff{i};
                    if isfield(C{t}, v) % do not overwrite with stuff that does not belong
                        I = isnan(C{t}.(v));
                        C{t}.(v)(I) = C{t}.imputed_pooled.(v)(I);
                    end
                end
            end
            
        end % extract_car_attributes
        
        function tab = rename_groupsummary_vars(tab, rename_dict)
            % the output from "groupsummary()" attaces prefix "mean_" to
            % all variable names. this removes that from specific ones. 
            
%             tab = tab_in; 
            
            vars = tab.Properties.VariableNames; 
            
            for i=1:numel(vars)
                v = vars{i}; 
                if (numel(v) > 5) && (strcmp(v(1:5), 'mean_')) % variables like "mean_inc" -> "inc" 
                    % overwrite this
                    v_new = v(6:end); 
                    tab.Properties.VariableNames(v) = {v_new}; 
                end
            end
            
            % rename variables
            if nargin < 2 
            rename_dict = { ...
                'd_scrap',               'scrap'; ...
                'bigcity1',             'urban'; ...
                'wd_sum',               'dist'; ...
                'realinc',              'inc'; ... 
                'fuelpriceannual',      'fuelp_annual';...
                'fuel_price_t1',        'fuelp_realized'; ...
                }; 
            end
            
            
            for i=1:size(rename_dict, 1)
                from = rename_dict{i, 1};
                if ismember(from, tab.Properties.VariableNames)
                    to = rename_dict{i, 2};
                    tab.Properties.VariableNames(from) = {to};
                end
            end    
            
        end
        
        function d = convert_demo_struct2table(dataAll, carAll, typeVar)
            % returns demographic data in a table formatm
            % NOTE: removes zeros from the car type variables 
            
            numT = numel(dataAll);
            zdemoLabel = dataAll{1}.zdemoLabel;
            
            count_notfound = 0; % keep track of how many, we fix this year
            
            % preallocate 
            d = cell(numT, 1); 
            
            for t = 1:numT
                d_ = dataAll{t};
                
                % 2. remove fields that are scalars or not NT*1 vectors
                d_ = rmfield(d_, 'NT');
                d_ = rmfield(d_, 'zdemoLabel');
                d_ = rmfield(d_, 'zdemo');
                
                % 4. recode cartypes
                
                % 4.a iecar
                I = d_.iecar >= 1; % owns a car
                cars = d_.iecar(I);
                d_.iecar(I) = carAll{t}.(typeVar)(cars);
                
                % 4.b icardrive
                I = d_.icardrive >= 1; % owns a car
                cars = d_.icardrive(I);
                d_.icardrive(I) = carAll{t}.(typeVar)(cars);
                
                % 4.c id
                I = d_.id >= 1;
                cars = d_.id(I);
                d_.id(I) = carAll{t}.(typeVar)(cars);
                
                % 5. to table
                d{t} = struct2table(d_);
                
                % delete cars that were not found
                I = isnan(d{t}.iecar) | isnan(d{t}.icardrive) | isnan(d{t}.id) ...  
                    | (d{t}.iecar == 0) | (d{t}.icardrive == 0) | (d{t}.id == 0); 
                count_notfound = count_notfound + sum(I); 
                d{t}(I,:) = []; % delete rows 
            end
            
            % report status
            tot_rows = 0;
            for t=1:numT
                tot_rows = tot_rows + dataAll{t}.NT;
            end
            fprintf('Car type 0s found for %d rows (of %d across all years)\n', count_notfound, tot_rows);
        end
        
        function tab = convert_car_struct2table(C, VARTYPE, numCarAges, numCarTypes, years)
            % INPUTS: 
            % VARTYPE: should be
            %   toplevel: extract top-level variables, one obs. per car
            %   byage: age-varying
            % C: numT-cell array of structs with car attributes 
            % years: numT-vector of calendar years
            
            numT = numel(C);
            c = cell(numT, 1);
            
            for t = 1:numT
                
                switch VARTYPE
                    case 'toplevel'
                        c_ = C{t};
                        c_ = rmfield(c_, 'imputed_pooled');
                        c_ = rmfield(c_, 'imputed_crossect');
                        c_ = rmfield(c_, 'byage');
                        c_ = rmfield(c_, 'nobs');
                        c_ = rmfield(c_, 'ntypes');
                        
                        % this struct is now ready to be converted to a
                        % table
                        
                    case 'byage'
                        c__ = C{t}.byage;
                        c__ = rmfield(c__, 'nobs');
                        
                        deleteThese = {'typegodkendnrnum1', 'cartype'};
                        for i=1:numel(deleteThese)
                            this_v = deleteThese{i}; 
                            if isfield(c__, this_v)
                                c__ = rmfield(c__, this_v); 
                            end
                        end
                        
                        c_ = struct(); 
                        
                        % must go from wide to long first 
                        vv = fieldnames(c__); % loop through all variables, assuming they are all nages*ntypes
                        for i=1:numel(vv)
                            v = c__.(vv{i}); 
                            assert(all(size(v) == [numCarAges, numCarTypes+1]), 'unexpected dimensions for variable %s', vv{i}); 
                            
                            % missings, part I: take from imputed values if
                            % they are available 
                            if isfield(C{t}.imputed_crossect.byage, vv{i})
                                v2 = C{t}.imputed_crossect.byage.(vv{i}); 
                                I = isnan(v); 
                                J = not(isnan(v2)); 
                                v(I & J) = v2(I & J); 
                            end
                            
                            % delete the column for type "-99" (all
                            % missings)
                            assert(all(isnan(v(:,1))), 'assumed first car was -99 with all missings'); 
                            v = v(:,2:end); 
                            
                            % missings, part II: fill using mean of nearby
                            % and, if that fails, just take the closest
                            % value 
                            I = isnan(v); 
                            if any(I, 'all')
                                v = fillmissing(v, 'movmean', 5);
                                v = fillmissing(v, 'nearest', 1);
                            end
                            assert(not(any(isnan(v), 'all')), 'Still missings after imputing... '); 
                             
                            % reshape, long->wide 
                            c_.(vv{i}) = v(:); 
                        end
                        
                        % add variable for the age and type of the car 
                        types = repmat(1:numCarTypes, numCarAges, 1); 
                        ages  = repmat((1:numCarAges)', 1, numCarTypes); 
                        
                        c_.typelist = types(:); 
                        c_.carage   = ages(:); 
                        
                    case 'imputed_pooled'
                        c_ = C{t}.imputed_pooled;
                        c_.typelist = C{t}.typelist;
                        warning('Not tested');
                    case 'imputed_crossect'
                        c_ = C{t}.imputed_crossect;
                        c_.typelist = C{t}.typelist;
                        warning('Not tested');
                    otherwise
                        error('Unexpected VARTYPE, "%s"', VARTYPE);
                end
                
                % to table: assumes all fields are vectors with same length
                c{t} = struct2table(c_);
                
                % car variables
                c{t}.t = t * ones(size(c{t}, 1), 1, 'int64');
            end
            
            tab  = vertcat(c{:});
            tab.Properties.VariableNames{'typelist'} = 'icar'; % rename
            I = tab.icar == -99;
            if any(I) % delete
                tab = tab(not(I), :);
            end
            tab.year = years(tab.t)';
            
            switch VARTYPE
                case 'toplevel'
                    tab = tab(:, [1,end-1,end,2:end-2]); % move last column ("t") to the front for convenience
                case 'byage'
                    tab = tab(:, [end-3, end-2, end-1, end, 1:end-5]);
            end
            
            assert(not(any(sum(ismissing(tab)))));
            
        end
        
        function ISCOL = table_has_col(tab, colname)
            ISCOL = ismember(colname, tab.Properties.VariableNames);
        end
        
        function tab = censor_table(tab_in, do_not_censor_these, DOPRINT)
            % INPUTS:
            %   table_in: table, must have variable "count"
            %   do_not_censor_these: cell array of strings that should not be
            %   overwritten
            
            % defaults 
            if nargin < 3
                DOPRINT = true; 
            end
            
            % checks
            assert(toy.table_has_col(tab_in, 'count'), 'Table should have field "count"');
            assert(iscell(do_not_censor_these), 'Input, "do_not_censor_these", should be a cell array (possibly empy, "{}").'); 
            
            % copy input 
            tab = tab_in;
            
            % find offending rows 
            I = (tab.count > 0) & (tab.count <= 5);
            if any(I)
                if DOPRINT
                    fprintf('Censoring %d of %d cells (%5.2f%% of obs.)\n', sum(I), numel(I), sum(tab.count(I))/sum(tab.count));
                end
                
                % loop through variables to be censored
                vv = tab.Properties.VariableNames;
                for i=1:numel(vv)
                    v = vv{i};
                    if not(ismember(v, do_not_censor_these)) % skipping excluded variables 
                        tab(I, v).Variables = nan(sum(I), 1);
                    end
                end
            else
                if DOPRINT
                    fprintf('No censoring required!!\n');
                end
            end
        end
        
    end
    
end
