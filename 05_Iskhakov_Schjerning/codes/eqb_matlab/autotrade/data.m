%% Descriptive Statistics data sets
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef data
properties (Constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Default parameters to be occasionally changed by hand
  %directory for the results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %properties
methods (Static, Access=public) %methods callable from the outside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mp, s, dta, car, demo, fuelp, est]=setup(datadir, years, loadmatfiles)
  % data.setup: Process data from excel files or read in .mat files
  %
  % SYNTAX [mp, s, dta, car, demo, fuelp, est]=data.read(datadir, loadmatfiles);
  %
  % INPUTS: 
  %     datadir:        path to location of data (optional)
  %                     default is '../../../../data-new/Raw/low/'
  %                 
  %     years:          1 x T vector of years to select 
  %
  %     loadmatfiles:   true/false 
  %                     if false (default) data from server is processed and saved in datadir as dta.mat 
  %                     if true (FASTER) previously saved dta.mat file is loaded 
  % OUTPUT
  %    car:             table with car characteristics (by cartype, carage, year) and variables
  %   
  %    demo             table with information about avg income, age, ect. (demographic group and year)
  %
  %    dta              table tabulated data with states, choices and counts (by year, states, choices)
  %.            dta has variables:
  %                 'category','year','count','vkt','avg_fuel_price,'scrap','tau','iecar',
  %                 'id','s_car_age','s_car_type','d_car_age','d_car_type','s_vintage','fuelp_annual'
  %
  %    fuelp            table with annual fuel data (by year)
  %
  %    est              table with regression results from driving equation (log-log specification)
  %
  % REMARK: loadmatfiles=true requires that process_data(datadir) has been run previously
  %                    
  % Currently available options for datadir: 
  %     1. (default):  '../../../../data-new/Raw/low/'  (4 car types, 8 consumer types). Here low = low-dimensional (relative to datasets coming later)
  %     2. (add new data paths here when data is ready)
  %
  % Guide on adding new data sets: 
  %     The "Raw" sub-folder contains the raw data files in xlsx format (in our case exports from Matlab on the server)
  %     Different sub-folders in ../../../../data-new/Raw/ will eventually contain different tabulations of the data (HH or car types)
  %     When data is added - run process_data(datadir, false) and add data dir to the help file

  % -------------------------------------------------------------------
  % 1. Read data and set up indexing 
  % -------------------------------------------------------------------
  
  [car, demo, dta, fuelp, est]=data.read(datadir, loadmatfiles, years);

  % Renaming some variables 
  dta.Properties.VariableNames({'d_km'})={'vkt'};   
  car.Properties.VariableNames({'car_type', 'car_age'})={'type','age'};        

  % Extract sample period and pool data (car characteristics are averaged over time) 
  [car, demo, dta, fuelp] = data.pool(car, demo, dta, fuelp, years); 

  % expand data set by splitting up cell conditional scrap-page
  [dta] = data.expand_scrap(dta);

  %% select data
  % delete the oldest car age: pileup of vintage cars in the raw data by construction. 
  max_car_age=22;
  I = dta.s_car_age<max_car_age & dta.d_car_age<max_car_age & ~isnan(dta.count) & ~isnan(dta.count_scrap);
  fprintf('Delete %d obs (%5.2f%%) where car age is %i years or older\n', sum(~I), 100*mean(~I),max_car_age); 
  dta=dta(I,:);

  % setup indexing (see help for setup_index)
  [mp, s, dta]  = data.setup_index(dta, car, demo);


  % -------------------------------------------------------------------
  % 2. Car characteristics
  % -------------------------------------------------------------------
  % compute holdings-weighted avg. of fueleff (mp.fe) 
  % FIXME: When fuel efficiency is pooled over years - it is not weighted by holdings

  % fueleff = car.fueleff; 
  % fe = nan(1, mp.ncartypes); 
  % for j=1:mp.ncartypes
  %     is = dta(dta.s_car_type == j, :).is; 
  %     n = dta(dta.s_car_type == j, :).count; 
  %     fe(j) = sum(n.*fueleff(is))/sum(n); 
  %     % fe(j) = mean(fueleff(is)); 
  % end

  fe = nan(1, mp.ncartypes); 
  for j=1:mp.ncartypes
      fe(j) = mean(car.fueleff(car.type==j)); % note: car has just one row per car (type, age) combination
  end
  mp.fe = num2cell(fe); 

  % registration taxes values 
  mp.vat         =  .25;    % value added tax
  mp.cartax_lo   =  1.05;   % registration tax (below kink, K_cartax_hi)
  mp.cartax_hi   =  1.8;    % registration tax (above kink, K_cartax_hi)
  mp.K_cartax_hi =  81;     % mp.K_cartax_hi before mp.cartax_hi tax kicks in

  share_gas = 0.80; 
  tax_fuel_of_pump = share_gas*0.5979 +  (1.-share_gas)*0.4730; 
  mp.tax_fuel    =  tax_fuel_of_pump/(1.0-tax_fuel_of_pump);    % proportional fuel tax 
  
  % prices 
  pnew = car(car.age == 1, :).price_new/1000;
  mp.pnew = num2cell(pnew'); % crucial that it becomes a row 

  mp.p_fuel = fuelp/1000;   % all prices should be measured in 1000 DKK 
  [mp.p_fuel_notax, mp.pnew_notax] = trmodel.price_notax(mp);

  assert(isscalar(mp.p_fuel), 'Got multiple fuel prices.'); 

  % -------------------------------------------------------------------
  % 3. Reduced form driving parameter
  % -------------------------------------------------------------------
  [mp.db]= data.driving_param(est, mp.ncartypes, mp.ntypes);

end

function [db]= driving_param(est, ncartypes, ntypes)
  % Parse driving equation parameters from estimates store table "est" 
  
  % logarithmic specification 
  if ismember('logppk', est.Properties.RowNames)
    db.specification = 'log-log'; 
    price_var = 'logppk'; 
  elseif ismember('ppk', est.Properties.RowNames)
    db.specification = 'lin-lin';
    price_var = 'ppk'; 
  else
      error('Cannot find "ppk" or "logppk" in row names of table est.'); 
  end

  % driving parameters (ntypes x 1)
  db.a1 = repcell({est('pt_car_age',   :).Estimate},ntypes,1); 
  db.a2 = repcell({est('pt_car_age_sq', :).Estimate},ntypes,1); 
 
  % demographic coefficients
  db.tau =cell(ntypes, 1); 
  db.pkm =cell(ntypes, 1); 
  for iz=1:ntypes
      % intercept 
      b0 = est('(Intercept)', :).Estimate; 
      if iz > 1
          row = ['tau_cat_' num2str(iz)];
          b1 = est(row, :).Estimate; 
      else
          % tau=0 is reference category
          b1 = 0; 
      end
      db.tau{iz} = b0 + b1; 
      
      % price 
      pkm0 = est(price_var, :).Estimate; 
      if iz > 1 
          row = [price_var ':tau_cat_' num2str(iz)];
          pkm1 = est(row, :).Estimate;
      else
          % tau=0 is reference category
          pkm1 = 0;
      end
      db.pkm{iz} = pkm0 + pkm1; 
  end

  % car specific coefficients 
  db.car =cell(1,ncartypes); 
  db.car{1} = 0.0; % reference category
  for ic=2:ncartypes
      row = ['pt_car_type_cat_' num2str(ic)];
      db.car{ic} = est(row, :).Estimate;
  end
end

function [car, demo, dta, fuelp, est]=read(datadir, loadmatfiles, years)
  % data.read: Process data from excel files or read in .mat files
  %
  % SYNTAX [car, demo, dta, fuelp, est]=data.read(datadir, loadmatfiles);
  %
  % INPUTS: 
  %     datadir:        path to location of data (optional)
  %                     default is '../../../../data-new/Raw/low/'
  %                 
  %     loadmatfiles:   true/false 
  %                     if false (default) data from server is processed and saved in datadir as dta.mat 
  %                     if true (FASTER) previously saved dta.mat file is loaded 
  % OUTPUT
  %    car:             table with car characteristics (by cartype, carage, year) and variables
  %   
  %    demo             table with information about avg income, age, ect. (demographic group and year)
  %
  %    dta              table tabulated data with states, choices and counts (by year, states, choices)
  %.                    dta has variables:
  %                     'category','year','count','vkt','avg_fuel_price,'scrap','tau','iecar',
  %                     'id','s_car_age','s_car_type','d_car_age','d_car_type','s_vintage','fuelp_annual'
  %
  %    fuelp            table with annual fuel data (by year)
  %
  %    est              table with regression results from driving equation (log-log specification)
  %
  % REMARK: loadmatfiles=true requires that process_data(datadir) has been run previously
  %                    
  % Currently available options for datadir: 
  %     1. (default):  '../../../../data-new/Raw/low/'  (4 car types, 8 consumer types). Here low = low-dimensional (relative to datasets coming later)
  %     2. (add new data paths here when data is ready)
  %
  % Guide on adding new data sets: 
  %     The "Raw" sub-folder contains the raw data files in xlsx format (in our case exports from Matlab on the server)
  %     Different sub-folders in ../../../../data-new/Raw/ will eventually contain different tabulations of the data (HH or car types)
  %     When data is added - run process_data(datadir, false) and add data dir to the help file
 

  fname_dta = sprintf('%sdta.mat', datadir);
  if nargin < 1
    datadir = '../../../../data/8x4/';
  end

  if nargin < 2 
    loadmatfiles = false;
  end
  
  if nargin < 3
      years = (1996:2009)'; 
  end

  if loadmatfiles
    disp(sprintf('Data read from previously saved .mat file:  %s', fname_dta))
    load(fname_dta)
    return
  end

  % read demographics
  disp(sprintf('Reading data from xlsx files in directory:  %s', datadir))
  demo = readtable([datadir 'demo.xlsx']); 

  % read in choice data
  T = numel(years); 
  dta = cell(T,1); 
  for t=1:numel(years)
      y = years(t); 
      fname = sprintf('%s/counts_%d.xlsx', datadir, y); 
      dta{t} = readtable(fname); 
  end
  dta = vertcat(dta{:}); % dta is now a table containing the full tabulated data including all years 
  
  % scrappage 
  
  
  scrap = readtable([datadir 'counts_scrap.xlsx']); 
  scrap.Properties.VariableNames({'count'}) = {'count_scrap'};
  keys = {'tau', 's_car_type', 's_car_age', 'year'};
  
  % merge on to the counts data with a left join 
  dta = join(dta, scrap, 'Keys', keys); 
  
  % car attributes 
  car = readtable([datadir 'car_attributes.xlsx']); 

  % Fuel prices 
  fuelp = readtable([datadir 'fuelp.xlsx']);

  % Read driving estimates 
  % Specification: linear-linear, model 4 (i.e. tau-dummies interacted
  % with price-per-km). 
  est = readtable([datadir 'est_lin_4.xlsx'], 'readrownames', true);

  % Merge on the fuel price (left join)
  dta = join(dta, fuelp(:, {'year', 'fuelp_annual'}));

  % Convert NaN counts -> 2: this is the median number, and the avg. is
  % 2.13, and the only possible values for these cells are {1,2,3,4,5}.
  dta.count(isnan(dta.count)) = 2; 

  % Convert driving from km/day to 1000 km / year 
  % Driving is simple, we just multiply. Driving parameter estimates:
  % because the model is linear-linear, we just multiply all coefficients
  % and standard errors (as we are simply scaling the y-variable). 
  dta.d_km = dta.d_km * 365.25/1000; 
  pars = est.Properties.RowNames; 
  for i=1:numel(pars)
      v = pars{i}; 
      est(v,:).Estimate = 365.25/1000 * est(v,:).Estimate; 
      est(v,:).SE       = 365.25/1000 * est(v,:).SE; 
  end  
  
  % Save the dataset
  save(fname_dta, 'car', 'demo', 'dta', 'fuelp', 'est')
  disp(sprintf('Data saved saved as .mat file:  %s', fname_dta))
end

function [dta_new] = expand_scrap(dta)
  if any(strcmp('d_scrap',dta.Properties.VariableNames))
    warning('d_scrap exists - data set already expanded')
    dta_new=dta;
    return
  end

  I=~isnan(dta.pr_scrap); 
  dta.d_scrap(:)=nan;

  dta_noscrap=dta;
  dta_noscrap.count=round((1-dta.pr_scrap).*dta.count);
  dta_noscrap.d_scrap(:)=0;

  dta_scrap=dta;
  dta_scrap.count=round(dta.pr_scrap.*dta.count);
  dta_scrap.d_scrap(:)=1;

  dta_new=[dta(isnan(dta.pr_scrap),:); dta_noscrap(~isnan(dta.pr_scrap),:); dta_scrap(~isnan(dta.pr_scrap),:)];

  fprintf('number of observations in original data: %d (cells = %d)\n',sum(dta.count), size(dta.count,1));
  fprintf('number of observations in expanded data: %d (cells = %d)\n',sum(dta_new.count), size(dta_new.count,1)); 
end

function lbl_types = hh_lbl_types(ntypes)
    switch ntypes
        case 16 % FIXME: DATA NOT YET AVAIABLE
            splits = { ...
                {'Low WD', 'High WD'}, ...
                {'Couple', 'Single'}, ...
                {'Q1 income', 'Q2 income', 'Q3 income', 'Q4 income'};
                };
        case 8
            splits = { ...
                {'Low WD', 'High WD'}, ...
                {'Couple', 'Single'}, ...
                {'Poor', 'Rich'};
                };
        otherwise
            error('household type labels not yet implemented for %d types', mp.ntypes);
    end
    
    lbl_types = cell(ntypes,1);
    i = 1;
    fprintf('Hard-coded household type names with %d types:\n', ntypes);
    for i1 = splits{1}
        for i2 = splits{2}
            for i3=splits{3}
                lbl_types{i} = sprintf('%s, %s, %s', i1{1}, i2{1}, i3{1});
                fprintf('%2d: %s\n', i, lbl_types{i});
                i = i+1;
            end
        end
    end
end

function [mp, s, dta] = setup_index(dta, car, demo)
  % data.setup_index:  updates the table dta with variables for states and choses (values and indexes), 
  %                    initialize the struct mp with dimensions and labels for cars and consumers 
  %                    and create the struct (s) with model indexing
  %
  % SYNTAX [dta, s, mp=data.setup_index(dta, car, demo);
  %
  % INPUTS:
  %    dta              table tabulated data with states, choices and counts (by year, states, choices)
  %.                    dta has variables:
  %                     'category', 'year', 'count', 'vkt', 'avg_fuel_price', 'p_scrap', 'fuelp_annual', 
  %    car:             table with car characteristics (by cartype, carage, year) and variables
  %   
  %    demo             table with information about avg income, age, ect. (demographic group and year)
  %
  % OUTPUT
  %    mp               struct (dimensions and labels for cars and consumers)
  % 
  %    s                structure with model indexing (see trmodel.index)
  %
  %    dta:             updated table with states, choices and counts
  %                     setup_index add the following variables to dta:
  %                     's_car_age','s_car_type', 'd_car_age', 'd_car_type', 'is', 'id', 'tau' 
  %   

  % -------------------------------------------------------------------
  % 1. Consumer types, demographic labels and population shares
  % -------------------------------------------------------------------
  mp.ntypes = numel(unique(demo.tau)); 

  mp.lbl_types = data.hh_lbl_types(mp.ntypes);

  % -------------------------------------------------------------------
  % 2. Car types and characteristics
  % -------------------------------------------------------------------
  mp.ncartypes = numel(unique(car.type)); 
  mp.ncarages = numel(unique(car.age)); 
  mp.lbl_cartypes = {'light, brown', 'light, green', 'heavy, brown', 'heavy, green'}; 
  % mp.abar_j0 = repcell({max(car.age+1)}, 1, mp.ncartypes); 
  mp.abar_j0 = repcell({max(car.age+1)}, 1, mp.ncartypes); 
  mp.ncars=mp.ncartypes*mp.ncarages;

  % -------------------------------------------------------------------
  % 3. Consumer type shares 
  % -------------------------------------------------------------------
  cc_ = groupsummary(dta(:, {'tau', 'count'}), 'tau', 'sum'); 
  mp.tw = cc_.sum_count ./ sum(cc_.sum_count); 


  % -------------------------------------------------------------------
  % 4. indexing
  % -------------------------------------------------------------------
  s=trmodel.index(mp, mp.abar_j0);

  % -------------------------------------------------------------------
  % 5. add indexes to dta
  % -------------------------------------------------------------------
  
  % recode 0 year old cars to new cars (these are cars that are trade withing first year)
  if any(dta.s_car_age==0)
    fprintf('%f observations with s_car_age=0 was recoded to s_car_age=1\n', sum(dta.count(dta.s_car_age==0)))
    dta.s_car_age(dta.s_car_age==0) = 1; 
  end

  dta.is = data.state_index(s, dta.s_car_age, dta.s_car_type) ;
  dta.id = data.decision_index(s, dta.d_car_age, dta.d_car_type) ; 
 

  % N = numel(dta.tau);
  % % initialize 
  % dta.is=nan(N,1);
  % dta.id=nan(N,1);


  % recode keep to purge for consumers without cars 
  % (in model it is not possible to keep no car, but possible to purge car - both results in staying with no car)
  dta.id(dta.id==s.id.keep & dta.is==s.is.nocar,:)=s.id.purge;

  owner= (dta.s_car_age >0 );
  keep= (dta.d_car_age == -1);  % index in id indicating keep 
  purge=(dta.d_car_age == -2);  % index in id indicating purge 
  
  % post-trade car 
  dta.pt_car_type = dta.d_car_type; 
  dta.pt_car_type(keep) = dta.s_car_type(keep); 
  dta.pt_car_type(purge) = -1; 
  
  dta.pt_car_age = dta.d_car_age; 
  dta.pt_car_age(keep & ~owner) = -1; 
  dta.pt_car_age(keep & owner) = dta.s_car_age(keep & owner) ;


  % dta.pt_car_age(keep & owner) = min(dta.s_car_age(keep & owner)+1, [s.abar_j{dta.pt_car_type(keep & owner)}]'); 


  dta.pt_car_age(purge) = -1; 
  
  % dta.ipt = data.state_index(s, dta.pt_car_age, dta.pt_car_type);
  dta.ipt = data.post_trade_index(s, dta.pt_car_age, dta.pt_car_type); 

  % -------------------------------------------------------------------
  % 5. checks
  % -------------------------------------------------------------------

  assert(not(any(dta.s_car_type == 0)), 'Car types should be -1 for "no car" and base 1 otherwise: found a zero'); 
  assert(not(any(isnan(dta.s_car_type))), 'Car types should be -1 for "no car" and base 1 otherwise: found a nan'); 
  
end

function id = decision_index(s, d_car_age, d_car_type) 
  % syntax id = data.decision_index(s, d_car_age, d_car_type) 
  % state indexes of car state
  id=nan(size(d_car_age,1),1);

  trader= (d_car_age >= 0);
  keep= (d_car_age == -1);  % index in id indicating keep 
  purge=(d_car_age == -2);  % index in id indicating purge 
  
  % decisions indexes (same index as in model)
  id(trader)= 2+d_car_age(trader) + [s.abar_j{d_car_type(trader)}]'.*(d_car_type(trader)-1);
  id(keep)=s.id.keep;
  id(purge)=s.id.purge;

  if max(id(trader))>=s.id.purge
    error('decision index for car traders is out of bounds, must be smaller s.is.purge - s.abar_j may be too small')
  end
end

function ipt = post_trade_index(s, pt_car_age, pt_car_type) 
  % syntax ipt = data.post_trade_index(s, pt_car_age, pt_car_type) 
  % post trade state indexes of car state
  ipt=nan(size(pt_car_age,1),1);

  owner= (pt_car_type >0);
  nocar= (pt_car_type == -1);   % index in is indicating no car
  
  % state index, dta.is (same index as in model)
  ipt(owner) = 1+pt_car_age(owner) + [s.abar_j{pt_car_type(owner)}]'.*(pt_car_type(owner)-1);
  ipt(nocar)=s.ipt.nocar;

  if max(ipt(owner))>=s.ipt.nocar
    error('state index for car owners is is out of bounds, must be smaller than s.is.nocar - s.abar_j may be too small')
  end

end

function is = state_index(s, s_car_age, s_car_type) 
  % syntax is = data.state_index(s, s_car_age, s_car_type) 
  % state indexes of car state
  is=nan(size(s_car_age,1),1);

  owner= (s_car_age >0);
  nocar=(s_car_age == -1);   % index in is indicating no car
  
  % state index, dta.is (same index as in model)
  is(owner) = s_car_age(owner) + [s.abar_j{s_car_type(owner)}]'.*(s_car_type(owner)-1);
  is(nocar)=s.is.nocar;

  if max(is(owner))>=s.is.nocar
    error('state index for car owners is is out of bounds, must be smaller than s.is.nocar - s.abar_j may be too small')
  end

end

function [car, demo, dta, fuelp] = pool(car, demo, dta, fuelp, years)
  % data.pool: pools data over selected years and aggregate car characteristics
  %      demographics and choice are not aggregated

  % pick out years 
  dta = dta(any(dta.year == years,2), :); 
  car = car(any(car.year == years,2), :); 
  demo = demo(any(demo.year == years,2), :);
  fuelp = fuelp(any(fuelp.year == years,2), :).fuelp_annual; 

  if numel(years)>1
      % average cars characteristics over years *unweighted* 
      car = groupsummary(car, {'type', 'age'}, 'mean'); 
      car = toy.rename_groupsummary_vars(car); 
      car(:, 'count') = []; 
      fuelp = mean(fuelp); % unweighted time average 
  end
end

function print_indexes(mp, s, dta) 
  format short g
  display('Randomly selected rows from from final data')
  dta_non_cesored=dta(dta.count>2.5,:);
  disp(sortrows(dta_non_cesored(randi(numel(dta_non_cesored.count),10,1),:)))

  display('Model indexes')
  disp(s)
  display('State space')
  disp(s.is)
  display('Decision space')
  disp(s.id)

  disp('Tabulation of key variables in the data')
  format long g
  disp('Variable defining consumer type')
  data.tabulate(dta, {'tau'}, mp.lbl_types);

  disp('Variables defining scrap decision - tabulation only for car-owners')
  if any(strcmp('d_scrap',dta.Properties.VariableNames))
    data.tabulate(dta, {'d_scrap'}, {'sell', 'scrap'}', dta.is ~= s.is.nocar);
    out=data.mean(dta, @(data) data.d_scrap, {'s_car_age'});
    disp('Scrap by car-age - for car owners only')
    disp(out)
  else
    disp('FIXME: d_scrap is missing from dta!!!')
  end

  disp('Distribution of car types (state)')
  data.tabulate(dta, {'s_car_type'}, {'nocar', mp.lbl_cartypes{:}}');
  disp('Distribution of car age (state)')
  data.tabulate(dta, {'s_car_age'});

  disp('Distribution of car types (decision)')
  data.tabulate(dta, {'d_car_type'}, {'purge', 'keep', mp.lbl_cartypes{:}}');
  disp('Distribution of car age (decision)')
  data.tabulate(dta, {'d_car_age'});

  disp('Distribution of car types (post trade)')
  data.tabulate(dta, {'pt_car_type'}, {'nocar', mp.lbl_cartypes{:}}');
  disp('Distribution of car age (post trade)')
  data.tabulate(dta, {'pt_car_age'});

  disp('Distribution of ipt (post trade index)')
  data.tabulate(dta, {'ipt'});

  disp('Index for states and decisions defined by car type and car age')
  if any(strcmp('is',dta.Properties.VariableNames))
    data.tabulate(dta, {'is'});
  else
    disp('WARNING: is is missing from dta!!!')
  end

  data.tabulate(dta, {'id'});
end

function [out] = mean(data, datavar, byvar, ifvar)

  if nargin>=4 & numel(ifvar)~=1
    data=data(ifvar==1,:);
  end

  data.datavar=datavar(data).*data{:,'count'};
  sumd= grpstats(data,byvar,'sum','DataVars',{'count' 'datavar'});
  sumd.mean_datavar= sumd.sum_datavar./sumd.sum_count; 
  
  out=table; 

  for i=1:numel(byvar)
    out.(byvar{i}) = sumd{:,i};  
  end
  out.mean = sumd.mean_datavar; 
end % end of table

function [out] = tabulate(data, byvar, lbl, ifvar)
  if nargin>=4
    data=data(ifvar==1,:);
  end
  sumd= grpstats(data,byvar,'sum','DataVars',{'count'});
  
  out=table; 
  for i=1:numel(byvar)
    out.(byvar{i}) = sumd{:,i};  
  end
  if nargin>=3
    if ~isempty(lbl)
      out.label=lbl;
    end
  end
  out.count =  round(sumd.sum_count, 0); 
  out.share = round(sumd.sum_count/sum(sumd.sum_count), 10); 

  if nargout==0
    fprintf('tabulation of %s\n', byvar{1})
    disp(out)
  end
end % end of table

function sumstat_table_households(dta, demo, mp, fname)
    % sumstat_table_demo
    if nargin<4 || isempty(fname) 
        fid = 1; 
    else 
        fid = fopen(fname, 'w'); 
    end
    tt = groupsummary(dta(:, {'tau', 'year', 'count'}), {'tau', 'year'}, 'sum');
    nobs = groupsummary(dta(:, {'tau', 'count'}), {'tau'}, 'sum');
    tt = tt(:, {'tau','year','sum_count'});
    j = join(tt, demo);
    j = removevars(j, {'count'});
    j.Properties.VariableNames('sum_count') = {'count'};
    j.real_inc = j.real_inc / 1000;

    vv = {'real_inc', 'single', 'wd', 'age', 'bigcity1', 'nkids'};
    xx = nan(mp.ntypes, numel(vv));
    for tau=1:mp.ntypes
        for i=1:numel(vv)
            v=vv{i};
            x = data.mean(j, @(data)data.(v), {'tau'});
            xx(:, i) = x.mean;
        end
    end
    
    % print 
    fprintf(fid, '%7s & %23s & %10s & %10s   & %10s   & %10s   & %10s   & %10s   & %10s   \\\\ \n', '$\\tau$', 'Name', 'N', vv{:});
    for tau=1:mp.ntypes
        fprintf(fid, '%7d & %23s & %10d & %10.2f & %10.2f & %10.2f & %10.2f & %10.2f & %10.2f \\\\ \n', tau, mp.lbl_types{tau}, nobs.sum_count(tau), xx(tau,:));
    end
    
    if fid ~= 1
        fclose(fid); 
    end
end % sumstat_table_households()

function sumstat_table_cars(dta, car, mp, fname) 
  % INPUTS: 
  %     fname: (optional) Filename to write to, e.g.
  %     'results/tables/sumstats_cars.tex' 
  
  agg_car_shares = data.tabulate(dta, {'s_car_type'}, {'nocar', mp.lbl_cartypes{:}}');

  % create a quick summary stats table 
  if nargin < 4 || isempty(fname) 
      fid = 1; 
  else 
      fid = fopen(fname, 'w');
  end
  
  car_tab = groupsummary(car, 'type', 'mean');
  fprintf(fid, '\\begin{tabular}{lrrrrr} \n\\toprule \n');
  fprintf(fid, ' %25s & %12s &  %12s   & %12s   & %12s & %12s \\\\ \n\\midrule \n', '', 'No car',                      mp.lbl_cartypes{:});
  fprintf(fid, ' %25s & %12d &  %12d   & %12d   & %12d & %12d \\\\ \n',                 'Obs.',                        agg_car_shares.count);
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Diesel share',           '',  car_tab.mean_diesel);
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Depreciation Factor',    '',  car_tab.mean_beta);  
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Weight (tons)',          '',  car_tab.mean_total_weight_tons);
  
  fprintf(fid, '\\midrule \n\\multicolumns{6}{c}{ \\emph{Variables used in the model} } \\\\ \n\\midrule \n');
  %fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n ',             'Price, new (1000 DKK)',  '',  car_tab.mean_price_new/1000); 
  %fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n ',             'Fuel efficiency (km/l)', '',  car_tab.mean_fueleff);
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Price, new (1000 DKK)',  '',  [mp.pnew{:}]); 
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Price, new excl. tax (1000 DKK)',  '',  [mp.pnew_notax{:}]); 
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Price, scrap (1000 DKK)','',  [mp.pscrap{:}]); 
  fprintf(fid, ' %25s & %12s &  %12.2f & %12.2f & %12.2f & %12.2f \\\\ \n',             'Fuel efficiency (km/l)', '',  [mp.fe{:}]);
  
  
  fprintf(fid, '\\bottomrule \n\\end{tabular} \n');
  
  if fid~=1
      fprintf('Wrote file: %s\n', fname);
      fclose(fid);
  end
end

end % end of methods
end % end of data
