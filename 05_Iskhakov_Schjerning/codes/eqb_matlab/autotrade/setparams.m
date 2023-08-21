% The main static class for setting paremeters 
% April 2021

classdef setparams
methods (Static)
  function mp = default()
    % Default parameters used for examples
    % SYNTAX :
    %   [mp] =setparams.default
    %
    % See also trmodel.setparams to get list of acctive parameters used in trmodel (all set to zero)
    
    mp.modeltype = 'reducedform'; % uses u_0 and u_a, not the "deep" structural parameters

    mp.lbl_cartypes={'Luxury', 'Normal'}; % car type labels (for plotting)
    mp.abar_j0 = {25,25};                 % car specific value of abar 
    mp.ncartypes=numel(mp.abar_j0);       % number of car types 

    mp.lbl_types={'Rich'; 'Poor'};  % consumer type labels (for plotting)
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    mp.sigma=1.;             % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=0.5; 
    mp.acc_0={-5.}; % = log(0.01), i.e. 1% accidents per year 
    mp.acc_a={0.0};  % 

    mp.mum = {.1; .3};     % marginal utility of money
    mp.psych_transcost = {0}; 
    mp.transcost = 1.5;      % fixed transaction cost for byuing car cars
    mp.ptranscost=0.0;     % proportional component of transaction costs, as a fraction of the price of the car the consumer buys
    mp.u_0={6  6.5};        % intercept in utility function (ntypes x ncartypes cell) 
    mp.u_a={-.5 -.475};      % slope coefficient on car-age in utility function (ntypes x ncartypes cell) 

    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat       =   25/100;   % value added tax
    mp.cartax_lo =  105/100;   % registration tax (below kink, K_cartax_hi)
    mp.cartax_hi =  180/100;   % registration tax (above kink, K_cartax_hi)
    mp.tax_fuel  =  100/100;   % proportional fuel tax 
    mp.K_cartax_hi = 81;       % mp.K_cartax_hi before mp.cartax_hi tax kicks in
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew={200,260};      % new car prices (cartype specific)
    mp.pscrap={1,1};        % scrapcar prices (cartype specific)
    mp.p_fuel=10.504;       % pricer of fuel (DKK/liter) from average fuel price in 2008 from our dataset
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
    
    % coefficent with price per kilometer by consumer type
    mp.fe     = {25,20};              % car specific fuel efficency (km/liter)
    mp.db.specification = 'linear';
    mp.db.pkm = {-26.2549; -26.3131}; % coefficient on pkm; 
    mp.db.car = {35.74};              % car type specific coefficient on car; 
    mp.db.tau = {-0.0968; -0.0723};   % coefficent with consumer type fixed effect by consumer type            
    mp.db.a1  = {-0.3444};            % coefficient on a*1; 
    mp.db.a2  = {0.00246};            % coefficient on a^2; 

    mp=trmodel.setparams(mp);
  end % end of setparams

  function mp = template()
    % Default parameters used for examples
    % SYNTAX :
    %   [mp] =setparams.template
    %
    % See also trmodel.setparams to get list of acctive parameters used in trmodel (all set to zero)

    % Copy this function above and parameters adjustments 

    mp=trmodel.setparams(mp)
  end % end of setparams

  function mp = first_submission()
    % Default parameters used for examples in the first submission of the paper 
    % ("Figure 3: Equilibrium ownership in a two consumer type economy")
    
    mp.lbl_cartypes={'Luxury', 'Normal'}; % car type labels (for plotting)
    mp.abar_j0 = {25,25};                 % car specific value of abar 
    mp.ncartypes=numel(mp.abar_j0);       % number of car types 

    mp.lbl_types={'Rich', 'Poor'};  % consumer type labels (for plotting)
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    mp.sigma=5;             % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.acc_0={-4.6};   % = log(0.01), i.e. 1% accidents per year 
    mp.acc_a={0.299};  % original had acc_a

    mp.mum = {1; 1.75};     % marginal utility of money
    mp.transcost =1.5;      % fixed transaction cost for byuing car cars
    mp.ptranscost=0.03;     % proportional component of transaction costs, as a fraction of the price of the car the consumer buys
    mp.u_0={60  65};        % intercept in utility function (ntypes x ncartypes cell) 
    mp.u_a={-5 -4.75};      % slope coefficient on car-age in utility function (ntypes x ncartypes cell) 

    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat       =   25/100;   % value added tax
    mp.cartax_lo =  105/100;   % registration tax (below kink, K_cartax_hi)
    mp.cartax_hi =  180/100;   % registration tax (above kink, K_cartax_hi)
    mp.tax_fuel  =  100/100;   % proportional fuel tax 
    mp.K_cartax_hi = 81;       % mp.K_cartax_hi before mp.cartax_hi tax kicks in
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew={200,260};      % new car prices (cartype specific)
    mp.pscrap={1,5};        % scrapcar prices (cartype specific)
    mp.p_fuel=10.504;       % pricer of fuel (DKK/liter) from average fuel price in 2008 from our dataset
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
    
    % coefficent with price per kilometer by consumer type
    mp.fe     = {25,20};              % car specific fuel efficency (km/liter)
    mp.db.specification = 'linear';
    mp.db.pkm = {-26.2549 -26.3131};  % coefficient on pkm; 
    mp.db.car = {35.74};              % car type specific coefficient on car; 
    mp.db.tau = {-0.0968; -0.0723};    %coefficent with consumer type fixed effect by consumer type            
    mp.db.a1  = {-0.3444};            % coefficient on a*1; 
    mp.db.a2  = {0.00246};            % coefficient on a^2; 

    mp=trmodel.setparams(mp);
  end % end of setparams

end %methods

end %class
