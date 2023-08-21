% Class that collects functions related to computing trade equilibrium at the secondary market for cars.
%
% This class relies on: 
%   trmodel     : defines model   
%   g_trmodel   : defines gradients of key functions in trmodel   
%   spp         : Solves the social planners problem to get initial gues of abar and prices for each car type   
%  
% January 2021

classdef equilibrium

properties (Constant)
  % Constant properties of the equilibrium class (can only be modified here)
  verbosity=0;        % 0 (no output), 1 (output on)
  options=optimoptions('fsolve','Jacobian','on','TolFun',1e-12,'TolX',1e-10, 'Display','off');
  tolerance=1e-6;     % solution tolerance: equilibrium not considered solved if max(abs(ed)) > equilibrium.tolerance
end

methods (Static)
  function [sol]=solve(mp, s, p0);
    % equilibrium.solve: solves for the auto market equilibrium given parameter structure 
    %
    % This function checks if norm of excess demand is zero and resolve using social planner prices if no equilibrium was found
    % (equilibrium.solve repeatedly calls function equilibrium.solve_once)
    %
    % SYNTAX: 
    %
    %   [sol]=equilibrium.[sol]=solve(mp, s, p0);  
    %   [sol]=equilibrium.[sol]=solve(mp, s);
    %   [sol]=equilibrium.[sol]=solve(mp);
    %
    % INPUTS:
    % mp          structure with model parameters (see trmodel.setparams)
    %
    % s:          (optional) structure with indexes for model. If s is not supplied it is computed using s=trmodel.index(mp)  
    %           
    % p0          (optional) starting values for p:  s.np x 1 price vector holding prices for each car-type and car age combination 
    %             where s.np  = sum(mp.abar_j0(:)-1) is the number of used cars (exluding clunkers). 
    %             p is indexed using s.ip: For example p0(s.ip{j}) gives the price vectire for car j             
    %
    %             If p0 is not supplied it can be initialized using p0 = equilibrium.price_init(mp,s) that computes 
    %             the prices from the homogeneous consumer economy. 
    %
    % OUTPUT:    
    %  sol: solution structure 
    %   struct with fields:
    %             p: [s.np x1 double]
    %            ed: [s.np x1 double]
    %             F: [1x1 struct]
    %         q_tau: {mp.ntypes x 1 cell}
    %             q: [s.ns x 1 double]
    %        ev_tau: {mp.ntypes x 1 cell}
    %       ccp_tau: {mp.ntypes x 1 cell}
    % ccp_scrap_tau: {mp.ntypes x 1 cell}
    %       ctp_tau: {mp.ntypes x 1 cell}
    %     delta_tau: {mp.ntypes x 1 cell}
    %         h_tau: {mp.ntypes x 1 cell}
    %  marketshares: [mp.ntypes x mp.ncartypes+1 double]
    %      exitflag: 1
    %
    %  where:
    %
    %   p         s.np x 1 dimensional equilibrium price vector (sol.p have the same structure as p0)
    %
    %   ed        s.np dimensional vector of excess demand for each car type/age combinations at prices, sol.p 
    %
    %   F:        struct with fields F.notrade and F.notrade holding sparse [s.ns x s.ns] extended age transition matrices  
    %             F.notrade is a sparse block-diagonal matrix with a block for each car and 1 for no car 
    %              - F.notrade governs aging for current car (identical to the Q matrix in the paper). 
    %              - F.trade governs aging for traded car 
    %              - see also F=trmodel.age_transition(mp, s)
    %
    %   q_tau     mp.ntypes dimensional cell array holding s.ns dimensional vectors of equilibrium 
    %             car quantity distribution for each consumer type in each state (not weighted with consumer types)
    %
    %   q         s.ns dimensional vector with aggreagte equilibrium car quantity distribution aggregated over consumers 
    %
    %   ev_tau    mp.ntypes dimensional cell array holding s.ns x 1 vectors with expected value function 
    %             for each consumer tau=1,..,mp.ntypes
    %
    %   ccp_scrap_tau: 
    %             mp.ntypes dimensional cell array holding s.ns x 1 vectors with scrapping probabilities 
    %             conditional on not keeping car
    %   ccp_tau   mp.ntypes dimensional cell array holding s.ns x s.nd matrix conditional choice probabilities 
    %             for each consumer tau=1,..,mp.ntypes, where
    %             the ccp matrix ccp_tau{t} for consumer t gives ccps for all combinations of states (rows) and choice (columns) 
    %              
    %   delta_tau mp.ntypes dimensional cell array holding s.ns x s.ns "trade transition probability matrices" for consumer tau
    %
    %   ctp_tau   mp.ntypes dimensional cell array holding s.ns x s.ns matrices ctp_tau{tau}=delta{tau}*F.notrade  
    %             
    %   h_tau     mp.ntypes dimensional cell array holding s.ns dimensional vectors of post trade equilibrium 
    %             holdings distribution for each consumer type in each state. h_tau{tau} is indexed by s.ipt
    %
    %  marketshares: 
    %             is mp.ntypes x mp.ncartypes+1 matrix of market shares held by consumers tau=1,..,mp.ntypes
    %             The first mp.ncartypes columns holds the market shares in car type j= 1: mp.ncartypes, 
    %             the last column holds the fraction of consumers in the outside good (i.e. with no car)
    %
    %   exitflag: 1 if the equilibrium is correctly solved to within equilibrium.tolerance, 0 otherwise
    %             
    %  see also equilibrium.solve_once, equilibrium.edf, equilibrium.print

    mp=trmodel.update_mp(mp);    % update parameters dependencies 


    if nargin==1
      s=trmodel.index(mp);
    end
    if nargin<=2
      p0 = equilibrium.price_init(mp,s);
    end

    for attempt=1:2 
      [sol]=equilibrium.solve_once(mp, s, p0);
      
      if mp.fixprices
          return
      end
      
      if (max(abs(sol.ed)) > equilibrium.tolerance || sol.exitflag<0) 
         p0 = equilibrium.price_init(mp,s); % reset prices and use spp solution

        if (equilibrium.verbosity>0)
          fprintf('Problem in the equilibrium solution: maximum absolute excess excess is %g\n',max(abs(sol.ed)));
          fprintf('exitflag from fsolve: %i\n',exitflag);
          if (attempt < 2)
            fprintf('Re-trying using the social planning solution as initial guess for prices\n');
          else
            fprintf('No equilibrium solution found in attempt %i: neither initial price guess or social planning starting guess worked.\n', attempt);
          end
        end
      else
        if (equilibrium.verbosity>0)
            fprintf('Equilibrium solution found in attempt %i: maximum absolute excess excess is %g\n',attempt, max(abs(sol.ed))); 
        end
        break; 
      end
    end
  end % end of while (~solved) loop

  function [p0] = price_init(mp, s)
    % Initialize price vector, p0, using social planner solution of used car prices
    % SYNTAX:    [p0] = equilibrium.price_init(mp, s)

    p0=nan(s.np,1);
    [abar_spp, p_spp, w_spp]=spp.solve_spp(mp);
    % use social planner solution as starting values for p0
    for j=1:mp.ncartypes
        [abar_j_spp, tau_j]=max([abar_spp{:,j}]);
        if abar_j_spp > mp.abar_j0{j};
          fprintf('WARNING: Social planner solution not consistent with mp.abar_j0: abar_j_spp=%d > mp.abar_j0=%d \nTruncating price vector - consider increasing mp.abar_j0\n', abar_j_spp, mp.abar_j0{j})
        end
        pj_spp=ones(1000,1)*mp.pscrap{j};
        pj_spp(1:abar_j_spp-1)=p_spp{tau_j, j}(2:abar_j_spp);
        p0(s.ip{j})=pj_spp(1:mp.abar_j0{j}-1); 
    end
  end

  function [sol]=solve_once(mp, s, p0);
    % equilibrium.solve_once: solves for the equilibrium (with no check for convergence)
    % equilibrium.solve_once us called equilibrium.solve and has same syntax, inputs and outputs
    % 
    % SYNTAX, INPUTS and OUTPUTS: same equilibrium.solve          
    % see also equilibrium.solve, equilibrium.print, equilibrium.edf

    mp=trmodel.update_mp(mp);  % update parameters dependencies 


    % starting values for p0
    % this creates the J*(abar-1) x 1 vector p containing the secondary market prices for each car type traded
    if nargin==1
      s=trmodel.index(mp);
    end
    if nargin<=2 || isempty(p0)
      p0 = equilibrium.price_init(mp,s);
    end
    
    if (equilibrium.verbosity>0);
      fprintf('Solving equilibrium with s.ns=%i car states, %i car types and %i consumer types\n',s.ns,mp.ncartypes,mp.ntypes);
    end

    % Initialize ev0 (if empty or if dimensions does not match)
    % Need to be global. 
    global ev0;
    if isempty(ev0) || numel(ev0)~=mp.ntypes || numel(ev0{1})~=s.ns;
      ev0=cell(mp.ntypes,1);
      for t=1:mp.ntypes;      
        ev0{t}=zeros(s.ns,1);
      end
    end

    % find the heterogeneous agent equilibrium for given abar_cell values and p0 starting values via Newton algorithm
    if mp.fixprices
      p=p0;
      exitflag = 1; 
      output=nan;
    else
      [p,ed,exitflag,output]=fsolve(@(p) equilibrium.edf(mp, s, p),p0, equilibrium.options);
    end

    [ed, ded, sol]=equilibrium.edf(mp, s, p);
    sol.p=p;
    sol.exitflag=exitflag;
    sol.output=output;

  end % end of solve

  function [ed, ded, sol]=edf(mp, s, p)
    % equilibrium.edf: evaluation of excess demand function (ed), its derivatives (ded) and the solution structure (sol)
    %        at a given price vector, p
    %        ed is evaluated by first solving the consumer DP problem and evaluating the choice probabilities
    %        and then determine the excess demand function. 
    %      
    % SYNTAX: 
    %   [ed, ded, sol]=equilibrium.edf(mp, s, p)
    %   [ed, ded]=equilibrium.edf(mp, s, p)
    %   [ed]=equilibrium.edf(mp, s, p)
    %
    % INPUTS:
    %   mp          structure with model parameters (see trmodel.setparams)
    %   s:          structure with indexes for model (see trmodel.index)
    %   p:          s.np x 1 dimensional equilibrium price vector (sol.p have the same structure as p0)
    % 
    % OUTPUTS:
    %   ed          structure with model parameters (see trmodel.setparams)
    %   ded:        structure with indexes for model (see trmodel.index)
    %   sol:        solution evaluated at prices p 
    %   
    %   sol is structure with fields
    %   
    %             p: [s.np x1 double]
    %            ed: [s.np x1 double]
    %             F: [1x1 struct]
    %         q_tau: {mp.ntypes x 1 cell}
    %             q: [s.ns x 1 double]
    %        ev_tau: {mp.ntypes x 1 cell}
    %       ccp_tau: {mp.ntypes x 1 cell}
    % ccp_scrap_tau: {mp.ntypes x 1 cell}
    %       ctp_tau: {mp.ntypes x 1 cell}
    %     delta_tau: {mp.ntypes x 1 cell}
    %         h_tau: {mp.ntypes x 1 cell}
    %  marketshares: [mp.ntypes x mp.ncartypes+1 double]
    %
    %  For detailed description fields in sol, see equilibrium.solve
    %
    %  see also equilibrium.solve, equilibrium.solve_once, equilibrium.print

    % note ev0 is a global cell array that stores the last values of the expected value
    % function to provide a starting point for the solution of the consumer problem
    global ev0

    % step 1: excess demand and it's subcomponents

    % compute the extended extended physical transition matrices F 
    [F] = trmodel.age_transition(mp, s);

    % initialize space
    ed_tau=nan(size(p,1),mp.ntypes);
    ev_tau=cell(mp.ntypes,1);
    ccp_tau=cell(mp.ntypes,1);
    ctp_tau=cell(mp.ntypes,1);
    q_tau=cell(mp.ntypes,1);
    delta_tau=cell(mp.ntypes,1); 
    deltaK_tau=cell(mp.ntypes,1); 
    deltaT_tau=cell(mp.ntypes,1); 
    delta_scrap=cell(mp.ntypes,1); 
    util_tau=cell(mp.ntypes,1); 
    ev_scrap_tau=cell(mp.ntypes,1); 
    ccp_scrap_tau=cell(mp.ntypes,1); 

    for t=1:mp.ntypes;
      if (mp.tw(t) > 0);
        % precompute utility 
        [util_tau{t}, ev_scrap_tau{t}, ccp_scrap_tau{t}] = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau

        % step 1: solve the DP problem for the price vector p for each consumer type
        bellman=@(ev) trmodel.bellman(mp, s, util_tau{t}, F, ev);

        [ev_tau{t}, ccp_tau{t}]= dpsolver.poly(bellman, ev0{t}, [], mp.bet);

        % trade, keep and scrap transition matrices
        [delta_tau{t}, deltaK_tau{t}, deltaT_tau{t}, delta_scrap{t}]=trmodel.trade_transition(mp, s, ccp_tau{t}, ccp_scrap_tau{t});    

        % state transition matrix
        ctp_tau{t}=delta_tau{t}*F.notrade; 

        % calculate the type-specific car distributions
        q_tau{t}=ergodic(ctp_tau{t}); 

        % step 2: compute contributions to excess demand for each consumer type       
        idx=[s.is.car_ex_clunker{:}];
        ed_tau(:,t)=deltaT_tau{t}(:,idx)'*q_tau{t} - (1-ccp_tau{t}(idx, s.id.keep)).*(1-ccp_scrap_tau{t}(idx,:)).*q_tau{t}(idx);

      end;

    end;

    % cumulate excess demand for all relevant ages of used cars that can be purchased, over all types of consumers
    ed=ed_tau*mp.tw; 

    % update starting values for ev. (used for next evaluation of edf)
    ev0=ev_tau;

    if (nargout >= 2);
      % equilibrium holdings distribution aggregated over consumers
      q=[q_tau{:}]*mp.tw(:); 

      % compute market shares and equilibrium holdings distribution post trade
      [h_tau] = equilibrium.h_tau(mp, s, q_tau, delta_tau);
      [marketshares] = equilibrium.marketshares(mp, s, h_tau);

      fields={'p', 'ed', 'F', 'q_tau', 'q', 'ev_tau', 'ccp_tau', 'ccp_scrap_tau', 'ctp_tau', 'delta_tau', 'h_tau', 'marketshares'};
      sol=cell2struct({p, ed, F, q_tau, q, ev_tau, ccp_tau, ccp_scrap_tau, ctp_tau, delta_tau, h_tau, marketshares}, fields, 2);

      % return excess demand and derivative of excess demand for the price vector p for each consumer type
      grad_dp=g_trmodel.partial_deriv(mp, s, sol,'prices');
      ded=grad_dp.ed;
    end 
  end % end of edf

  function [F, ev_tau, ccp_tau, ccp_scrap_tau]=solve_dp(mp, s, p)

    % solve_dp.m: solves for each consumer type's optimal trading strategy via DP at a given price vector
    %        It first solves the consumer DP problem and then evaluates the choice probabilities
    %      
    % SYNTAX: [F,  ev_tau, ccp_tau, ccp_scrap_tau]=equilibrium.solve_dp(mp, s, p)

    % note ev0 is a golbal cell array that stores the last values of the expected value
    % function to provide a starting point for the solution of the consumer problem

    global ev0

    % compute the physical transition probability matrices for all car types (conditional on trading / or not trading)
    [F] = trmodel.age_transition(mp, s);

    % initialize space
    ev_tau=cell(mp.ntypes,1);
    ccp_tau=cell(mp.ntypes,1);

    for t=1:mp.ntypes;
      if (mp.tw(t) > 0);
        % precompute utility 
        [util_tau{t}, ev_scrap_tau{t}, ccp_scrap_tau{t}] = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau

        % step 1: solve the DP problem for the price vector p for each consumer type
        bellman=@(ev) trmodel.bellman(mp, s, util_tau{t}, F, ev);

        [ev_tau{t}, ccp_tau{t}]= dpsolver.poly(bellman, ev0{t}, [], mp.bet);

      end;

    end;

    % update starting values for ev. (used for next evaluation of edf)
    ev0=ev_tau;
  end % end of solve_dp

  function [h_tau] = h_tau(mp, s, q_tau, delta_tau)
    % equilibrium.h_tau: distribution of ownership after trading 
    %
    % SYNTAX: [h_tau] = equilibrium.h_tau(mp, s, q_tau, delta_tau)
    % INPUTS: see definitions under equilibrium.solve
    % OUTPUT:
    %  h_tau  mp.ntypes dimensional cell array holding s.ns dimensional vectors of post trade equilibrium 
    %         holdings distribution for each consumer type in each state. h_tau{tau} is indexed by s.ipt

    h_tau=cell(mp.ntypes,1);    
    for tau=1:mp.ntypes;
      qpt=delta_tau{tau}'*q_tau{tau};
      for j=1:mp.ncartypes;
        h_tau{tau}(s.ipt.new{j},1)=mp.tw(tau)*qpt(s.is.clunker{j});
        h_tau{tau}(s.ipt.used{j},1)=mp.tw(tau)*qpt(s.is.car_ex_clunker{j});
        h_tau{tau}(s.ipt.nocar,1)=mp.tw(tau)*qpt(s.is.nocar);
      end
    end
  end

  function [marketshares] = marketshares(mp, s, h_tau)
    % equilibrium.marketshares: market shares after trading 
    %
    % SYNTAX: [marketshares] = equilibrium.marketshares(mp, s, h_tau)
    % INPUTS: see definitions under equilibrium.solve
    % OUTPUT:
    % mp.ntypes x mp.ncartypes+1 matrix of market shares held by consumers tau=1,..,mp.ntypes
    %           The first mp.ncartypes columns holds the market shares in car type j= 1: mp.ncartypes, 
    %           the last column holds the fraction of consumers in the outside good (i.e. with no car)

    marketshares=nan(mp.ntypes,mp.ncartypes+1);
    for tau=1:mp.ntypes;
      for j=1:mp.ncartypes;
        marketshares(tau,j)=nansum(h_tau{tau}(s.ipt.car{j}));
      end
      marketshares(tau,mp.ncartypes+1)=h_tau{tau}(s.ipt.nocar,1);
    end
  end

  function []=print(mp, s, p);
    % SYNTAX equilibrium.print(mp, s, p)
    [~, ~, sol]=equilibrium.edf(mp, s, p);  % compute excess demand one final time

    for i=1:mp.ncartypes;
        fprintf('market share of car type %i %g\n',i,sum(sol.q(s.is.car{i})));
    end

    fprintf('market share of outside good: %g\n',sol.q(s.ns));
    fprintf('%15s %15s\n', 'prices','excess demands');
    fprintf('%15g %15g\n', [p sol.ed]');
    
    if (max(abs(sol.ed)) > equilibrium.tolerance);
      warning('Maximum value of excess demand is %g possibility that this is not an equilibrium\n',max(sol.ed));
    end
  end

end % end of methods
end % end of classdef
