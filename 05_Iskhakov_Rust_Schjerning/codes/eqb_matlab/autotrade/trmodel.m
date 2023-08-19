% The main static class for the AUTOMOBILE TRADE model.
% Contains all needed model parts.
% 
% April 2021
%

classdef trmodel

methods (Static)
  function [mp] =setparams(mp0)
    % Standard parameters of the model
    % SYNTAX :
    %   [mp] = trmodel.setparams       set all parameters to zeros
    %   [mp] = trmodel.setparams(mp0)  update patameters with parameters in mp0
    %
    % NOTE: When adding new parameters to trmodel, they must be added here 
    % DO NOT ENTER VALUES OF PARAMETERS - JUST ZEROS (or values that turn off the parameter)
    % See also setparams.default to get a default parameters

    % ******************************************
    % switches
    % ******************************************
    mp.es=1;          % use model with endogenous scrapage if 1 
    mp.fixprices=0;   % set to zero, to compute model solution at a given set of prices without solving for equilibrium 
    mp.ll_scrap=true; % add scrap decisions to likelihood if true
    % ******************************************
    % misc. parameters
    % ******************************************
    mp.dk_population =2.5;     % number of danish households in millions

    % ******************************************
    % consumer types and car types
    % ******************************************
    % Some parameters below are consumer type specific (in rows) and car type specific (in columns) 
    % If number of entries is smaller than ntypes, last entry is repeated)

    mp.lbl_cartypes={'Luxury', 'Normal'}; % car type labels (for plotting)
    mp.abar_j0 = {25,25};                 % car specific value of abar 
    mp.ncartypes=numel(mp.abar_j0);       % number of car types 

    mp.lbl_types={'Rich', 'Poor'};  % consumer type labels (for plotting)
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
    mp.ntypes=numel(mp.tw);         % number of types of consumers

    % ******************************************
    % discount factor (beta)
    % ******************************************
    mp.bet=.95;  

    % ******************************************
    % accident parameters (alpha)
    % ap=1./(1+exp(-apv)), where apv = mp.acc_0{ct} + mp.acc_a{ct}.*a + mp.acc_even{ct}*(1-mod(a,2));     
    % ******************************************
    mp.acc_0={-10};
    mp.acc_a={0};
    mp.acc_even={0};

    % ******************************************
    % parameters of GEV distribution (theta)
    % ******************************************
    mp.sigma=1;        % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=1;       % the extreme value scale parameter for the idiosyncratic extreme value shocks affecting the sale vs scrappage decision
 
    % ******************************************
    % utility patameters (theta)
    % ******************************************
    mp.mum = {1};   % marginal utility of money

    mp.phi = {0.01}; % scale of driving utilty (set low to exclude driving)

    % transactions costs parameters 
    mp.transcost =0;      % fixed transaction cost for byuing car cars
    mp.ptranscost=0;      % proportional component of transaction costs, as a fraction of the price of the car the consumer buys
    mp.psych_transcost      = {0};      % utility cost of bying a car
    mp.psych_transcost_nocar= {0};      % additional utility cost of bying a car, when not owning a car
    mp.nocarsc     =0;    % additional search cost incurred by a consumer who has no car
    mp.tc_sale     =0;    % set the fixed component of the transaction cost to a seller (inspection/repair cost to make a car qualified to be resold)
    mp.tc_sale_age =0;    % coefficient of age on the total sales costs
    mp.ptc_sale    =0;    % set the proportional component o the transaction cost to a seller (inspection/repair cost to make a car qualified to be resold)
    mp.tc_sale_even=0;    % incremental sales transactions costs during inspection years (even years) necessary to pass inspection before car can be sold

    % car utility parameters (see also trmodel.u_car) 
    % Reduced form car utility mp.u_0{tp, car}+mp.u_a{tp, car}*car_age + mp.u_a_sq{tp, car}*car_age.^2
    mp.u_0={0};         % intercept in utility function (ntypes x ncartypes cell) 
    mp.u_a={0};         % slope coefficient on car-age in utility function (ntypes x ncartypes cell) 
    mp.u_a_sq = {0};    % car age squared 

    mp.convexutility=0; % this is a switch, if 1 then mp.u.a_sq is forced to be positive via the transformation exp(mp.u_a_sq) in u_car

    mp.u_og={0};        % utilty of outside good  (ntypes x 1 cell) 
    mp.u_even={0};      % utility during even inspection years (expect to be a disutility or negative value due to the disutility of having car inspected)        
    % tax policy parameters (must be before specification of prices of fuel and cars)
    mp.vat         =  0;   % value added tax
    mp.cartax_lo   =  0;   % registration tax (below kink, K_cartax_hi)
    mp.cartax_hi   =  0;   % registration tax (above kink, K_cartax_hi)
    mp.tax_fuel    =  0;   % proportional fuel tax 
    mp.K_cartax_hi =  0;   % mp.K_cartax_hi before mp.cartax_hi tax kicks in

    % ******************************************
    % Prices before taxes (prices after taxes are computed when solving the model)
    % if car or consumer type parameters have size smaller than ncartypes and ntypes, the last element is repeated.
    % ******************************************
    mp.pnew_notax={100};    % new car prices (cartype specific)
    mp.pscrap_notax={0};    % scrapcar prices (cartype specific)
    mp.p_fuel_notax=5;      % pricer of fuel (DKK/liter) from average fuel price before tax
    
    % ******************************************
    % Parameters of reduced form driving equation 
    % x = mp.db.pkm{tp}*pkm{car} + mp.db.car{car}+ mp.db.tau{tp} +mp.db.a1{tp}*a+mp.db.a2{tp}*a.^2;
    % ******************************************
    % Reduced form driving equation coefficients are stored in db structure 
    mp.fe={20};       % car specific fuel efficency (km/liter)
    mp.db.specification = 'linear'; % 'linear' or 'log-log' 
    mp.db.pkm = {0};  % coefficient on pkm; 
    mp.db.pkm = {0};  % coefficient on pkm; 
    mp.db.car = {0};  % car type specific coefficient on car; 
    mp.db.tau = {0};  % coefficent with car type fixed effect by consumer type            
    mp.db.a1  = {0};  % coefficient on a*1; 
    mp.db.a2  = {0};  % coefficient on a^2; 
    
    % ******************************************
    % Reduced form or structural form
    % ******************************************
    % Structural parameters can be identified from reduced from parameters for car demand and driving
    % by using trmodel.update_structural_par to obtain "structural" utility parameters mp.sp. 
    % Because of implicit dependence on type specific marginal utility of money, all these coefficients have to be consumer type specific
    % coefficent with price per kilometer by consumer type
    
    mp.modeltype = 'reducedform'; 
    % if mp.modeltype == 'reducedform': reudced form parameters mp.u_0, mp.u_a, and mp.u_a_sq are used in trmodel.u_car. 
    % if mp.modeltype == 'structuralform': structural parameters mp.sp, mp.mum are used in trmodel.u_car
    % 
    % If the model allows for driving, the reduced form parameters mp.u_0, mp.u_a, and mp.u_a_sq 
    % are not policy invariant since they on fuel-prices (and marginal utility of money and other structural parameters)
    % So to run counterfactuals where fuelprices are changed you need to set mp.modeltype = 'structuralform';
    %
    % Estimation strategy: 
    %   First set mp.modeltype = 'reducedform' to estimate reduced form parameters during estimation. 
    %   Then set mp.modeltype = 'structuralform' for running counter-factuals that changes fuel prices. 

    % ********************************************************************************
    % update parameters
    % ********************************************************************************
    % update mp with mp0 (default values are over written by input values mp0)
    if nargin>0
      mp.db=combinestructs(mp.db, mp0.db);
      mp=combinestructs(mp, rmfield(mp0,'db'));
    end

    % update endogenous parameters
    mp=trmodel.update_mp(mp);
    [mp.p_fuel_notax, mp.pnew_notax, mp.pscrap_notax] = trmodel.price_notax(mp);
  end % end of setparams

  function [p_fuel_notax, pnew_notax, pscrap_notax] = price_notax(mp)

    % This function computes prices before taxes given after tax prices and tax rates in mp
    % When implementing new taxes/subsidies remember to update the corresponding after tax function below. 

    % fuel prices before taxes
    p_fuel_notax=mp.p_fuel/(1+mp.tax_fuel);  
          
    pnew_notax = cell(1, mp.ncartypes); 
    for car=1:mp.ncartypes 

      pnew_notax{car} = trmodel.pcar_notax(mp.pnew{car}, mp.K_cartax_hi, mp.cartax_lo, mp.cartax_hi, mp.vat);
    
      % scrap price before tax: not currently implemented 
      if nargout >=3 
        pscrap_notax{car} = mp.pscrap{car};
      end
    end 
    
  end

  function pnew_notax = pcar_notax(pcar, kink_threshold, cartax_low, cartax_high, vat)
    % pcar_notax(): new car prices before registration taxes and VAT.
    %
    % INPUTS: (all scalar floats)
    %   pcar: consumer price of the car 
    %   kink_threshold: the threshold where the registration tax rate shifts 
    %   cartax_high: marginal tax rate above kink
    %   cartax_low: marginal tax rate below kink
    %   vat: Value Added Tax (should be in [0; 1).)
    % 
    % OUTPUT: 
    %   pnew_notax: car price without any taxes

    assert(cartax_low <= cartax_high, 'Low-bracket rate is above high-bracket rate: this sounds like an error (unless you are doing a fancy counterfactual, then comment this assertion out)');
    assert(vat < 1.0, 'Expecting VAT = 0.25 (or at least < 1.0). Delete this assertion if you are trying a crazy counterfactual and know what you are doing.')
    
    % cutoff = tax amount paid if a car price is precisely equal to the kink threshold value
    
    % how much is paid if the pre-tax car price (with VAT) puts it
    % precisely at the cutoff 
    price_paid_at_cutoff = (1+cartax_low)*kink_threshold;

    % compute the final inversion separately depending on whether we are
    % above or below the cutoff 
    if pcar <= price_paid_at_cutoff 
        pnew_notax = pcar / ((1+cartax_low)*(1+vat)); 
    else 
        numer = pcar + (cartax_high - cartax_low)*kink_threshold;
        denom = (1+vat)*(1+cartax_high); 
        pnew_notax = numer/denom; 
    end
        
  end

  function [p_fuel, pnew, pscrap] = price_aftertax(mp)
    % price_aftertax(): This function computes prices after taxes given before tax prices and tax rates in mp
    % When implemneting new taxes/subisdies remember to update the correspoding after tax function below. 

    p_fuel  = mp.p_fuel_notax*(1+mp.tax_fuel); % scalar 
          
    pnew    = cell(1, mp.ncartypes); 
    pscrap  = cell(1, mp.ncartypes); 
    
    for icar=1:mp.ncartypes

      pnew_notax = trmodel.pcar_after_passthrough(mp.pnew_notax{icar}, mp, icar); 

      pnew{icar} = trmodel.pcar_aftertax(pnew_notax, mp.K_cartax_hi, mp.cartax_lo, mp.cartax_hi, mp.vat);

      if nargout >= 3
          pscrap{icar} = mp.pscrap_notax{icar}; % currently, no special treatment 
      end
    end
  end

  function passthrough = set_up_passthrough(mp_baseline, rate)

    passthrough = struct(); 
    passthrough.pnew_notax_baseline = mp_baseline.pnew_notax;
    passthrough.pnew_baseline = mp_baseline.pnew; 
    passthrough.rate = rate; 
    
    passthrough.cartaxes_baseline = struct(); 
    passthrough.cartaxes_baseline.K_cartax_hi   = mp_baseline.K_cartax_hi; 
    passthrough.cartaxes_baseline.cartax_lo     = mp_baseline.cartax_lo; 
    passthrough.cartaxes_baseline.cartax_hi     = mp_baseline.cartax_hi; 
    passthrough.cartaxes_baseline.vat           = mp_baseline.vat; 
    
  end

  function pcar = pcar_after_passthrough(pcar_raw, mp, icar, DOPRINT)
      % pcar_after_passthrough(): adds a mechanical firm-response to the
      % raw pre-tax price. The amount to add is set so that a target
      % passthrough rate is achieved 
      %
      % INPUTS: 
      %     pcar_raw (double): car price, pre-tax and pre-passthrough 
      %     mp (struct): model parameters 
      %     icar (int): index for the car [XXX drop input??] 
      %     DOPRINT (bool, optional): print details 
      % 
      % OUTPUT: 
      %     pcar (double): car price after firm "markup" has adjusted 
      % 
      
      if nargin < 4
          DOPRINT = false;
      end
      
      if isfield(mp, 'passthrough')
          assert(icar>=1 && icar<=mp.ncartypes);
          assert(isfield(mp.passthrough, 'pnew_notax_baseline'));
          assert(isfield(mp.passthrough, 'rate'));
          
          % 1. compute price that *would* have prevailed absent any changes
          % in firm behavior 
          pcar_at_full_passthrough = trmodel.pcar_aftertax(pcar_raw, mp.K_cartax_hi, mp.cartax_lo, mp.cartax_hi, mp.vat);
          mp0 = mp.passthrough.cartaxes_baseline; % tax rates in the baseline
          pcar_baseline_full_pasthrough = trmodel.pcar_aftertax(pcar_raw, mp0.K_cartax_hi, mp0.cartax_lo, mp0.cartax_hi, mp0.vat);
          delta_tax = pcar_at_full_passthrough - pcar_baseline_full_pasthrough;
          
          % 2. change in manufacturer price
          % NOTE: "-1" because manufacturers move *opposite* of the policy makers intended direction
          delta_firm_price = (-1) * (1-mp.passthrough.rate) * delta_tax;
          
          % 3. final price before taxes get applied
          pcar = pcar_raw + delta_firm_price;
          
          if DOPRINT
              fprintf('At full passthrough: p0 = %5.2f to p1 = %5.2f: implied delta tax = %5.2f\n', pcar_baseline_full_pasthrough, pcar_at_full_passthrough, delta_tax); 
              fprintf('Requested passthrough = %5.2f%% => delta p raw = %5.2f\n', 100.0*mp.passthrough.rate, delta_firm_price); 
              fprintf('Final result: p0 = %5.2f -> p with passthrough = %5.2f\n', pcar_raw, pcar); 
          end
          
      else
          % nothing to do
          pcar = pcar_raw;
      end
  end

  function pcar_incl_taxes = pcar_aftertax(pcar_notax, K_cartax_hi, cartax_lo, cartax_hi, vat)
      % pcar_notax(): new car prices before registration taxes
      %
      % INPUTS: 
      %   pcar_notax: consumer price of the car 
      %   K_cartax_hi: the threshold where the registration tax rate shifts 
      %   cartax_hi: marginal tax rate above kink
      %   cartax_lo: marginal tax rate below kink
      %   vat: Value Added Tax 
      % 
      % OUTPUT: 
      %   pcar_incl_taxes: car price including all taxes

      assert(cartax_lo <= cartax_hi, 'hi/low bracket rates reversed: this could be an error!'); 
      assert(vat <= 1.0, 'VAT should be in [0;1] (unless you are doing a crazy large counterfactual)');
      
      if pcar_notax*1.25 <= K_cartax_hi
          % no top-tax will be paid 
          pcar_incl_taxes = (1+cartax_lo)*1.25*pcar_notax; 
      else
          % price is in the top-bracket 
          pcar_incl_taxes = (1+cartax_hi)*1.25*pcar_notax - (cartax_hi - cartax_lo)*K_cartax_hi; 
      end
      
  end

  function [s] = index(mp, abar_j)
    % This function set up the state space for the trmodel. 
    % 
    % SYNTAX: [s] = trmodel.index(mp, abar_j)
    % 
    % INPUT: 
    %   mp:     structure with model parameters (see trmodel.setparams)
    %   abar_j: mp.ncartypes cell array with max age (forced scrap age) for each car j=1,..,mp,ncartypes
    %
    % OUTPUT:  
    %   s:  a struct with the following elements:
    %
    %     id: struct that holds decision indexes groups of decisions. 
    %     is: struct that holds state indexes groups of states
    %     ip: mp.ncartypes x 1 cell of indexes used car prices. 
    %     ns: number of states
    %     nd: number of decisions
    %     np: number of price parameters
    %     tr: struct with sub-indexes for transitions
    %     ipt: struct that holds indexes for the post trade distribution
    %
    %  Example: 
    %  s=trmodel.index(mp, {16,22}) gives
    %   s = 
    %     struct with fields:
    %       abar_j: {[16]  [22]}
    %           id: [1x1 struct]
    %           is: [1x1 struct]
    %           ip: {[1 2 ... 15]  [16 17 ... 36]}
    %           ns: 39
    %           nd: 40
    %           np: 36
    %           tr: [1x1 struct]
    %          ipt: [1x1 struct]
    %  where
    % 
    %   s.id has the following fields
    %                keep: 1
    %               trade: {[2 3 ... 17]  [18 19 ... 39]}
    %          trade_used: {[3 4 ...  17]  [19 20 ...  39]}
    %           trade_new: {[2]  [18]}
    %                 age: [NaN 0 1 ...  15 0 1 2 ...  21 NaN]
    %               purge: 40
    %    
    %   s.is has the following fields
    %                  car: {[1 2  ... 16]  [17 18 ... 38]}
    %       car_ex_clunker: {[1 2  ... 15]  [17 18 ... 37]}
    %              clunker: {[16]  [38]}
    %                  age: [1 2 ... 16 1 2 ... 22 NaN]
    %                nocar: 39
    %       
    %               choice: [3 4 ... 17 2 19 20 ... 39 18 40]
    %                state: [1 2 ...  39]
    %  
    %  s.ip mp.ncartypes x 1 cell with the following elements
    %       s.ip{1} = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
    %       s.ip{2} = [16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36]
    %  
    %  s.ipt has the following fields
    %       car: {[1 2 ... 16]  [17 18 ... 38]}
    %       nocar: 39
    %       new: {[1]  [17]}
    %       used: {[2 3 ... 16]  [18 19 ... 38]}
    %       age: [0 1 ... 15 0 1 2 ... 21 NaN]    %
    % 
    %  usage: 
    %   s.id.keep, s.id.trade{j},s.trade_used{j}, s.trade_new{j}, and s.id.purge are choice indexes that 
    %   for example used for indexing columns in matrices such as ccp_tau, and decision specific utility
    %
    %   s.is.car{j}, s.is.car_ex_clunker{j},s.is.clunker{j}, and s.is.nocar are state indexes that 
    %   for example used for indexing rows in matrices such as ccp_tau, ev_tau, q, q_tau{tau}.
    %
    %   For example, ccp_tau{tau}(s.is.car_ex_clunker{j}, s.id.keep) gives the vector of probabilities 
    %   of keeping a car of type j with ages 1,2,..,abar_j{j}-1. The age of such a car can also be accessed by 
    %   s.is.age(s.is.car_ex_clunker{j})  
    %
    % Some example descriptions of subindices: 
    % 
    % s.ip  for each car type, the indices for each traded car age (1,....abar_j{i}-1), i=1...ncartypes sequentially ordered
    %       for that the union of these indices over all car types is the sequence 1,2,....,abar_j{1}+...abar_j{ncartype}-ncartypes, 
    %
    % s.is.car (mp.ncartypes cell array): 
    %   s.is.car{j} is a vector holding the indexes in the overall ordering of car states for car type j - except the no car state. 
    % s.is.car_ex_clunker (mp.ncartypes cell array): 
    %   s.is.car_ex_clunker{j} is a vector holding the indices in the overall ordering of car states for car type j
    %    - except the oldest car (clunker) and the no car state. 
    % s.is.clunker (mp.ncartypes cell array): 
    %   s.is.clunker{j} is a scalar holding the index for the oldest car of type j
    % s.is.nocar (scalar): 
    %   index for the state of having no car
    % s.is.age (s.ns vector): 
    %   s.is.age vector of car ages corresponding to state indexes. Age of no car is coded as NaN


    if nargin==1
        s.abar_j=mp.abar_j0; 
    else
        s.abar_j=abar_j;
    end
 
    % state and decision indexes
    s.id.keep=1;    % keep decision
    ei=0;
    eip=0;
    for j=1:mp.ncartypes
      % indexes states and trade decisions (including new car)
      si=ei+1;
      ei=ei+s.abar_j{j};

      % car owner state
      s.is.car{j}=[si:ei];
      s.is.car_ex_clunker{j}=s.is.car{j}(1:end-1);
      s.is.clunker{j}=s.is.car{j}(end);

      % trade decision 
      s.id.trade{j}=1+(si:ei); % add 1 to account for keep beeing the first index
      s.id.trade_used{j}=s.id.trade{j}(2:end);
      s.id.trade_new{j}=s.id.trade{j}(1);

      % corresponding car age values
      s.is.age(s.is.car{j})=(1:s.abar_j{j})';   % not possible to enter period with new car
      s.id.age(s.id.trade{j})=(0:s.abar_j{j}-1)'; % not possible to buy clunker
      
      %indexes for used car prices and excess demand
      sip=eip+1;
      eip=eip+s.abar_j{j}-1; 
      s.ip{j}=(sip:eip); 
    end

    s.ns=ei+1;         % number of states      (add one for no car state)
    s.nd=ei+2;             % number of decisions   (add one for keep and purge)
    s.np=eip;              % number prices          

    s.is.nocar=s.ns;    % nocar state
    s.id.purge=s.ns+1;    % purge decision

    s.is.age(s.is.nocar)=nan;
    s.id.age(s.id.keep)= nan;
    s.id.age(s.id.purge)=nan;

    % Indices for trade transition matrix
    % build index for chocies correnspoding to columns in deltaT matrix (see definition of deltaK in paper)
    s.tr.choice=[];
    for j=1:mp.ncartypes; % loop over cars 
        s.tr.choice=[s.tr.choice s.id.trade{j}(2:end) s.id.trade{j}(1)]; % first used, then new
    end
    s.tr.choice= [s.tr.choice s.id.purge];  % purge in last column 

    % index for all rows and columns in delta
    s.tr.state=[s.is.car{:} s.is.nocar]; 

    % index for next period car state 
    s.tr.next_acc=nan(1, s.ns);    % state after accident
    s.tr.next_keep=s.tr.state+1;   % state after keep (current car ages one period) - no accident
    s.tr.next_trade=s.tr.state;    % state after trading - no accident
    for j=1:mp.ncartypes;
      s.tr.next_acc(s.is.car{j})=s.is.clunker{j};
      s.tr.next_keep([s.is.clunker{j}])=[s.is.car{j}(1)];
    end
    s.tr.next_keep(s.is.nocar)=s.is.nocar;
    s.tr.next_acc(s.is.nocar)=s.is.nocar;

    % Age transition matrices conditional on accidents
    s.Q.no_acc= sparse(s.tr.state,s.tr.next_keep,1,s.ns,s.ns);
    s.Q.acc= sparse(s.tr.state,s.tr.next_acc,1,s.ns,s.ns);
    s.dQ=s.Q.acc-s.Q.no_acc;
    s.F.no_acc= sparse(s.tr.state,s.tr.state,1,s.ns,s.ns);
    s.F.acc=s.Q.acc;
    s.dF=s.F.acc-s.F.no_acc;

    % indexes for post trade distribution (before aging)
    s.ipt.car=s.is.car; 
    s.ipt.nocar=s.is.nocar; 
    for j=1:mp.ncartypes;
      s.ipt.new{j}=s.ipt.car{j}(1);
      s.ipt.used{j}=s.ipt.car{j}(2:end);
    end 
    s.ipt.age=s.is.age -1; 
  end % end of index

  function [u, ev_scrap, ccp_scrap]=utility(mp, s, tau, price_vec)
    %   utility:  Function to compute utility for each consumer in each state and for each decision 
    %             and its derivative wrt to used car prices. 
    %
    %   SYNTAX:     [u, ev_scrap, ccp_scrap]=trmodel.utility(mp, s, price_j)
    %
    %   INPUTS:
    %     mp          structure with model parameters (see trmodel.setparams)
    %     s:          structure with indexes for dimensions. see trmodel.index
    %     tau:        household type 
    %     price_j0    (optional) starting values for price_j. price_j0 is a mp.ncartypes dimensional cell array 
    %
    %
    %   OUTPUTS          
    %   u:      a cell array of dimension mp.ntypes whose elements are matrices that store the utility
    %           received from trading for a car of type j and age a for a car of type j' and age d, or "purge"
    %           trade the current car for the outside good, or to keep the car, for each consumer type.
    %
    %           u_tau{t} has s.ns rows and s.nd columns and rows is indexed by s.is and columns by s.id. 
    %
    %           Rows correspond to the different states 
    %           i.e the type/age of the currently held car or the state of no car. 
    %           The last row corresponds to the state of having no outside good. ()
    %
    %           The columns correspond to the different possible decisions a consumer has.
    %           The different possible decisions are: keep, trade and purge. 
    %
    %           The first column is the utility of the decision to keep (s.id.keep=1) 
    %           for each possible age of car (or outside good in last column)
    %
    %           NOTE: The utility of keep contain a missing value in in rows where the age equals the scrappage age (is.clunker{j})
    %           for each car type, since a consumer is not allowed to keep a car once it reaches the scrappage age, and for a consumer
    %           who has no car, the law row of the u_tau provides the utility of continuing to stay with the outside good (keeping outside good)
    % 
    %           The last column, s.id.purge=s.nd, is the utility of the decision to purge 
    %            (i.e. choose the outside good, i.e. move to having no car)
    %
    %           The columns in between, from 2 to s.nd-1 or [s.id.trade{:}], correspond to the decision to trade 
    %           for a car of age a, a=0,...,abar{j}-1 for each of the ncartypes of cars
    % 
    %  ev_scrap  the extra option value from being able to scrap instead of sell a used car (only non-zero if mp.es=1)
    %            This is the s.ns x 1 vector of integrated value of scrap choice relative to selling choice
    %            i.e. logsum(mp.mum{tau}*(pscrap-psell), 0)  
    %  ccp_scrap: 
    %           s.ns x 1 vector of scrapping probabilities (conditional on not keeping)
    %
    %           Gillingham, Iskhakov, Munk-Nielsen, Rust and Schjerning, September 2020

    % update the prices in the vector price_vec to the cell array price_j
    for i=1:mp.ncartypes;
      price_j{i}=price_vec(s.ip{i});
    end;

    % net sales price after transactions cost 
    [psell]=trmodel.psell(mp, s, price_j); 

    
    % net purchase price after transactions cost 
    [pbuy]=trmodel.pbuy(mp, s, price_j); 
  
    u=nan(s.ns,s.nd); % utility: number of states by number of decisions (see state_decision_index for information about of indexing)

    % utility of scrap option (when selling)
    [ccp_scrap, ev_scrap]=trmodel.p_scrap(mp, s, price_j, tau);

    for j=1:mp.ncartypes
      % utility of keeping (when feasible)
      u(s.is.car_ex_clunker{j}, s.id.keep) =  trmodel.u_car(mp, s.is.age(s.is.car_ex_clunker{j}), tau, j);
      
      % utility of trading 

      u(: , s.id.trade{j})           =  trmodel.u_car(mp, s.id.age(s.id.trade{j}), tau, j) ...
                                      + ev_scrap-mp.mum{tau}*(pbuy(:,s.id.trade{j})-psell)-mp.psych_transcost{tau};
    end % end loop over car-types

    % utility of purging (i.e. selling/scrapping the current car and not replacing it, i.e. switching to the outside good)
      
    u(:, s.id.purge)                  =  mp.u_og{tau} + mp.mum{tau}*psell + ev_scrap;   

    % additional psych transactions cost and monetary search cost incurred by a consumer who has no car
    u(s.is.nocar , [s.id.trade{:}])   =  u(s.is.nocar , [s.id.trade{:}]) ... 
                                             -  mp.psych_transcost_nocar{tau} - mp.mum{tau}*mp.nocarsc; 
    
    u([s.is.clunker{:}] , s.id.keep)  =  nan; % not possible to keep clunker
    u(s.is.nocar , s.id.keep)         =  nan; % not possible to keep no car
  end % end of utility

  function [ccp_scrap, ev_scrap]=p_scrap(mp, s, price_j, tau)
    %   p_scrap:  Function to compute expected utility of the option to scrap rather than sell  
    %             for each consumer in each car state, the probability of scrapping,  and the
    %             derivatives with respect to mum{t} and the sell-side transaction cost parameters
    %
    %   SYNTAX:   [ccp_scrap, ev_scrap]=trmodel.p_scrap(mp, s, psell, dpsell, tau)

    psell=trmodel.psell(mp, s, price_j); 

    pscrap=nan(s.ns,1);
    for j=1:mp.ncartypes
      % tc=[price_j{j}; mp.pscrap{j}] *mp.ptc_sale+mp.tc_sale+mp.tc_sale_age*s.is.age(s.is.car{j})'+mp.tc_sale_even*(1-mod(s.is.age(s.is.car{j})',2));
      % pscrap(s.is.car{j})=mp.pscrap{j} - tc;
      pscrap(s.is.car_ex_clunker{j})=mp.pscrap{j};
    end

    ev_scrap = (mp.es==1)*logsum([mp.mum{tau}*(pscrap-psell), zeros(s.ns,1)], mp.sigma_s); 
    ccp_scrap=1-exp(-(ev_scrap)/mp.sigma_s);

    for j=1:mp.ncartypes
      ccp_scrap(s.is.clunker{j})=1;
      ev_scrap(s.is.clunker{j})=0;
    end

  end % end of p_scrap
  
  function [uv]=u_car(mp, car_age, tau, car)
    % INPUTS: 
    %     car_age: vector 
    %     tau: household type, scalar index
    %     car: car type, scalar index 

    pkm = mp.p_fuel./mp.fe{car}*1000; % p_fuel is measured in 1000 DKK/l, mp.fe{j} is measured as km/l, but pkm was in DKK/km in regression)
    
    switch mp.modeltype
      case 'structuralform' % structural form - requires that mp.sp are set by trmodel.update_structural_par
        % evaluate utility 
        uv = mp.sp.alpha_0(tau,car) + mp.sp.alpha_a(tau,car) .* car_age + mp.sp.alpha_a_sq(tau,car) * car_age.^2 ... 
        - 1./(2*mp.sp.phi(tau)) .* (max(0,mp.sp.gamma_0(tau,car) + mp.sp.gamma_a(tau) .* car_age - pkm .* mp.mum{tau})).^2;
      case 'reducedform'
        if (mp.convexutility)
        uv=mp.u_0{tau, car}+mp.u_a{tau, car}*car_age + exp(mp.u_a_sq{tau, car})*car_age.^2;
        else
        uv=mp.u_0{tau, car}+mp.u_a{tau, car}*car_age + mp.u_a_sq{tau, car}*car_age.^2;
        end

      otherwise 
      error('Unexpected reduced form type, "%s".', mp.modeltype); 
    end

    % add the (dis)utility from car inspections in even years 
    % (after the first inspection at age 4)
    inspection=(1-mod(car_age,2)).*car_age>=4; % dummy for inspection year
    uv= uv+mp.u_even{tau,car}.*inspection; 
  end % end of u

  function [pbuy] = pbuy(mp, s, price_j)     
    pbuy=nan(1, s.nd); 
    for j=1:mp.ncartypes
      pbuy(s.id.trade_new{j})=  mp.transcost + [mp.pnew{j}*(1+mp.ptranscost)];  
      pbuy(s.id.trade_used{j})= mp.transcost + [price_j{j}*(1+mp.ptranscost)];  
      pbuy(s.id.purge)=0;
    end
  end

  function [psell] = psell(mp, s, price_j)     
    % psell is the selling price net of seller-side transactions costs
    psell=nan(s.ns,1);
    for j=1:mp.ncartypes
      % car_age=s.is.age(s.is.car{j})'; 
      % inspection=(1-mod(car_age,2)).*car_age>=4;
      % tc=[price_j{j}; mp.pscrap{j}] *mp.ptc_sale+mp.tc_sale+mp.tc_sale_age*car_age+mp.tc_sale_even*inspection;
      % psell(s.is.car{j})=[price_j{j}; mp.pscrap{j}]-tc;

      car_age=s.is.age(s.is.car_ex_clunker{j})';
      inspection=(1-mod(car_age,2)).*car_age>=4; % dummy for inspection year
      psell(s.is.car_ex_clunker{j})=[price_j{j}*(1-mp.ptc_sale)-mp.tc_sale-mp.tc_sale_age*car_age-mp.tc_sale_even*inspection];
      psell(s.is.clunker{j})=[mp.pscrap{j}];
      psell(s.is.nocar)=0;

    end
  end % end of trade_cost

  function [x]=driving(mp, a, tau, j)
    % Reduced form or structural form driving equation 
    % x = mp.db.pkm{tp}*pkm{car} + mp.db.car{car}+ mp.db.tau{tp} +mp.db.a1{tp}*a+mp.db.a2{tp}*a.^2;
    % trmodel.driving: Driving in km/day 
    % Syntax: [x]=trmodel.driving(mp, a, tau, j)
    % OUTPUTS: 
    %   x: optimal driving 1000 km / year 
    %
    % NOTE: if mp.sp parameters are updated appropriately, the two
    % modeltypes should give identical driving (by construction). 
    pkm = mp.p_fuel./mp.fe{j}*1000; % p_fuel is measured in 1000 DKK/l, mp.fe{j} is measured as km/l, but pkm was in DKK/km in regression)
    switch mp.modeltype          
      case 'structuralform' % structural driving equation   
        x = -1./mp.sp.phi(tau).*max(0,mp.sp.gamma_0(tau,j) + mp.sp.gamma_a(tau) .* a - pkm .* mp.mum{tau}); 
      case 'reducedform' % use the estimated OLS equation directly 
        x =  mp.db.car{j} + mp.db.tau{tau} + mp.db.a1{tau}*a + mp.db.pkm{tau}*pkm; 
      otherwise 
        error('Unexpected value in mp.modeltype, "%s".', mp.modeltype); 
    end
  end

  function [ap,dap]=acc_prob_j(mp,a,ct,abar)
    % acc_prob.m: probability of an accident as a function of age and type of vehicle
    %             
    %  apv is either scalar (if a is scalar integer) or vector (if a is a vector of integer car ages)
    %      containing the accident probabilities for a car of car type ct (integer index) and 
    %      a consumer of type tau (integer index). The parameters are in cell arrays mp.acc_0,mp.acc_a,mp.acc_even
    %      which have mp.ntypes of rows (ntypes number of consumer types, if this is 1, then all consumer types
    %      have the same accident probability) and mp.ncartypes columns (if there is only 1 column all cartypes have
    %      the same accident probability parameters).  mp.acc_even is coefficient on a dummy for whether the car age
    %      is even, to capture higher likelihood a car involved in an accident is scrapped rather than repaired due to 
    %      the even-year inspections of cars in Denmark, which may impose a higher cost of repairs to make vehicle pass inspection
    %      NOTE: dependence of accident probabilities on consumer type tau not yet implemented due to ramifications elsewhere in the code:
    %      the physcial_transition function would also have to be consumer-type specific, as well as its three outputs.
    %
    % dap is the corresponding vector of gradients with respect to the parameters of the accident probability model described 
    %       above
    %         
    % Note: all cars are scrapped at age abar and an accident that totals a car
    %       is treated as equivalent to scrappage of the car. So we set the accident
    %       probability for a car of age abar-1 to 1 since it will become a car of
    %       age abar next period and be scrapped.

    f_apv= @(mp) mp.acc_0{ct} + mp.acc_a{ct}.*a + mp.acc_even{ct}*(1-mod(a,2));
    ap=1./(1+exp(-f_apv(mp)));
    ap(a >= abar-1)=1.0;

    if (nargout > 1)
      [pvec0, mp]=estim.mp2pvec(mp, mp.pnames);
      dap=ap.*(1-ap).*estim.matgrad_mp(f_apv, mp, pvec0, 'alpha');
    end

  end % end of acc_prob

  function [F] = age_transition(mp, s)
    % age_transition.m: age transition probability matrices 
    % 
    %  SYNTAX:   [F]= trmodel.age_transition(mp, s)
    %
    % INPUTS: 
    %   mp:     structure with model parameters
    %   s:      structure with state and decision indexes (generated by trmodel.states) 
    %
    % OUTPUS: 
    %   F:      structure with age transition matrices - taking account for accidents
    %
    %   The fields of F are block-diagonal matrices with a block for each car and 1 for no car  
    %     F.notrade:   s.ns x s.ns extended physical transition matrix (Similar to Q matrix in paper)
    %     F.trade:     s.ns x s.ns state transition matrix conditional on trading a car 

    accprob=zeros(s.ns,1);
    for j=1:mp.ncartypes;
        accprob(s.is.car{j})=trmodel.acc_prob_j(mp,(0:s.abar_j{j}-1)',j, s.abar_j{j});
    end

    F.notrade =s.Q.no_acc+s.dQ.*accprob(s.tr.next_keep);
    F.trade   =s.F.no_acc+s.dF.*accprob(s.tr.state);
    F.accprob=accprob;
  end

  function [delta, deltaK, deltaT, delta_scrap] = trade_transition(mp, s, ccp, ccp_scrap)
    % Keeping transition probabilities
    deltaK=sparse(1:s.ns,1:s.ns,ccp(:,s.id.keep),s.ns,s.ns);

    % Trading transition probabilities
    deltaT     = ccp(s.tr.state,s.tr.choice);

    %  trade transition probability matrix
    delta=deltaT +  deltaK; 

    if nargout>3
      delta_scrap =(1-ccp(:, s.id.keep)).*ccp_scrap;
    end

  end % end of trade_transition

  function [ev1, ccp, dev] = bellman(mp, s, util, F, ev)
    % bellman.m:    Implements the Bellman operator ev1=Gamma(ev) for the consumer's problem of trading cars,
    %               to maximize discounted utility in the presence of a secondary market for cars
    %               This function returns the *expected value function* ev and is a fast vectorized
    %               program that also calculates choice probabilities implied by current expected value ev
    %               and the derivative of the Bellman operator
    %
    % INPUTS: 
    %   mp:     structure with model parameters
    %   s:      structure with state and decision indexes (generated by trmodel.states) 
    %   util:   utility matrix (n.s x n.s+1) (states in rows, decision in columns)
    %   Q:      s.ns x s.ns extended physical transition matrix (block-diagonal matrix with a block for each car and 1 for no car)
    %   F:      s.ns x s.ns state transition matrix conditional on trading a car (block-diagonal matrix with a block for each car and 1 for no car)
    %   ev:     s.ns x 1 vector of expected values 
    %
    % OUTPUS: 
    %   ev1:    s.ns  x 1 vector of updated expected values after evaluating the Bellman operator 
    %   ccp:    s.ns  x s.ns + 1 matrix of conditional choice probabilities  
    %   dev:    s.ns  x s.ns matrix of derivatives of Bellman operator
    % 
    %  Fixed point of Bellman can be found by calling dpsolver
    % 
    %  SYNTAX: 
    %     [ev,ccp,dev]= dpsolver.poly(@(ev) trmodel.bellman(mp, s, util, F, ev), ev0, mp, mp.bet);
    % 
    % See also:
    %   dpsolver.poly
    %
    %   John Rust, Georgetown University, January 2019

    v=nan(s.ns,s.nd);        % choice specific value functions (states by decision)
    
    % Calculate the expected discounted utility of keeping car
    % Keeping is only valid for car ages  a=1,...abar-1 (i.e not in clunker state and nocar state)
    % Infeasible choices are coded as NaNs in values and as 0 for ccps
    % ev_keep  =  ev(s.tr.next_keep)   + (ev(s.tr.next_acc)-ev(s.tr.next_keep)). *accprob(s.tr.next_keep);
    % ev_trade =  ev(s.tr.next_trade)  + (ev(s.tr.next_acc)-ev(s.tr.next_trade)).*accprob(s.tr.next_trade)

    v([s.is.car{:}],s.id.keep)=    util([s.is.car{:}],s.id.keep) + mp.bet*F.notrade([s.is.car{:}],:)*ev; 

    % calculate the values for all options involving buying a car
    v(:, [s.id.trade{:}])         =util(:,[s.id.trade{:}]) + (mp.bet*F.trade([s.is.car{:}],:)*ev)';   

    % calculate the values for the purge decision (ev of no car)
    v(:,s.id.purge)               =util(:,s.id.purge)   + mp.bet*ev(s.is.nocar);                     % calculate the values for the purge decision (or keep having
    v(isnan(util))=nan;
    ev1=logsum(v, mp.sigma); 
    if nargout>1
      ccp=exp((v-ev1)/mp.sigma);
      ccp(isnan(ccp))=0;  % restore the nans in the elements of ccp_cell corresponding to infeasible choices
    end

    if nargout>2
      %  compute trade transition probability matrix
      delta=trmodel.trade_transition(mp, s, ccp);

      % derivative of Bellman operator wrt ev
      dev=mp.bet*delta*F.notrade;
    end
  end % end of bellman

  function mp = update_mp(mp)

    if (~iscell(mp.u_even))
         mp.u_even={mp.u_even};
    end


    % update model dependent parameters when changing mp.ntypes or mp.ncartypes
    % if cells are smaller than mp.ntypes or mp.ncartypes the last element is repeated

    mp.db.car               =repcell(mp.db.car,1,mp.ncartypes); 
    mp.db.pkm               =repcell(mp.db.pkm,mp.ntypes,1); 
    mp.db.tau               =repcell(mp.db.tau,mp.ntypes,1); 
    mp.db.a1                =repcell(mp.db.a1,mp.ntypes,1); 
    mp.db.a2                =repcell(mp.db.a2,mp.ntypes,1); 

    % car-specific parameters
    param_j       = {'lbl_cartypes', 'abar_j0','acc_0','acc_a','acc_even','fe', 'pnew_notax', 'pscrap_notax'};
    for k=1:numel(param_j)
      mp.(param_j{k}) =repcell(mp.(param_j{k}),1, mp.ncartypes);
    end  

    % household type-specific parameters
    param_tau     = {'lbl_types','mum','phi','u_og','psych_transcost','psych_transcost_nocar'};
    for k=1:numel(param_tau)
      mp.(param_tau{k}) =repcell(mp.(param_tau{k}),mp.ntypes,1);
    end  

    % car-and-household specific parameters
    param_tau_j   = {'u_0','u_a','u_a_sq','u_even'};
    for k=1:numel(param_tau_j)
      mp.(param_tau_j{k}) =repcell(mp.(param_tau_j{k}),mp.ntypes, mp.ncartypes);
    end  

    % the fundamental prices are p_fuel_notax, pnew_notax, plus the tax rates 
    [mp.p_fuel, mp.pnew, mp.pscrap] = trmodel.price_aftertax(mp);
  
    if (mp.ntypes == 1)
       mp.tw=1;
    end

    if (abs(sum(mp.tw)-1.0) > 1e-12)
       fprintf('Error in trmodel.setup: mp.tw vector does not sum to 1\n');
    end
    
  end % end of setup

  function sp = update_structural_par(mp)
      % SYNTAX: sp = trmodel.update_structural_par(mp)
      % NOTE: this should be called after estimation, but not whenever
      % changing e.g. fuel taxes (then the structural parameters are
      % fixed). 
      
      sp = struct(); 
      
      % 1. phi coefficients (on squared driving)
      mum = cell2mat(mp.mum);
      sp.phi = mum ./ cell2mat(mp.db.pkm); % ntypes*1 vector
      
      % 2. compute gamma coefficients
      d0 = cell2mat(mp.db.car) + cell2mat(mp.db.tau); % adding a row and column vector -> a ntypes*ncartypes matrix
      d1 = cell2mat(mp.db.a1);
      sp.gamma_0 = - d0 .* sp.phi; % ntypes*ncartypes
      sp.gamma_a = - d1 .* sp.phi; % ntypes*ncar
      
      % 3. compute alpha coefficients
      pkm=mp.p_fuel./cell2mat(mp.fe)*1000; % p_fuel is measured in 1000 DKK/l, mp.fe{j} is measured as km/l, but pkm was in DKK/km in regression)
      if (mp.convexutility)
      sp.alpha_a_sq = exp(cell2mat(mp.u_a_sq)) + 1./(2*sp.phi)      .* (sp.gamma_a.^2);
      else
      sp.alpha_a_sq = cell2mat(mp.u_a_sq) + 1./(2*sp.phi)      .* (sp.gamma_a.^2);
      end
      sp.alpha_a    = cell2mat(mp.u_a)    + sp.gamma_a./sp.phi .* (sp.gamma_0 - mum.*pkm);
      sp.alpha_0    = cell2mat(mp.u_0)    + 1./(2*sp.phi)      .* (sp.gamma_0 - mum.*pkm).^2;
      
      % 4. evaluate utility
      %uv = sp.alpha_0 + sp.alpha_a * car_age + sp.alpha_a_sq * car_age.^2 ...
      %    - 1/(2*sp.phi) * (sp.gamma_0 + sp.gamma_a * car_age - pkm * mum);

      % 5. check to make sure no negative driving:  issue warnings if this is found

      for t=1:mp.ntypes
       for c=1:mp.ncartypes
          if (sp.gamma_0+sp.gamma_a*(mp.abar_j0{c}-1) < mum(t)*pkm)
            fprintf('Warning: update_structural_par  negative driving predicted for household type %i and cartype %i\n',t,c);
          end
       end
      end
  end
  

end %methods

end %class
