% Class that collects gradients necessary to compute equilibrium trade at the 
% secondary market: 
% functions: derivatives of excess demand and its subcomponents
%
% This class relies on: 
%   trmodel       : defines model   
% 
% April 2021

classdef g_trmodel

methods (Static)

  function [grad_dp_dmp, grad_dp, grad_dmp] = total_deriv(mp, s, sol)
    % total_deriv: compute total and partial derivatives of a number of elements in solutions struct 
    % Derivatives of utility and accident probabilities are computed numerically 
    % and thus this function is fully general wrt. to changes in specification of utility

    % SYNTAX 1: [grad_dp_dmp, grad_dp, grad_dmp] = g_trmodel.total_deriv(mp)
    % SYNTAX 2: [grad_dp_dmp, grad_dp, grad_dmp] = g_trmodel.total_deriv(mp, s)
    % SYNTAX 3: [grad_dp_dmp, grad_dp, grad_dmp] = g_trmodel.total_deriv(mp, s, sol)
    %
    % Syntax 1 and 2 require solving for equilibrium prices (time consuming)
    % 
    % INPUTS: 
    %   mp          structure with model parameters (see trmodel.setparams for definitions)
    %   s:          structure with indexes for model (see s=trmodel.index for definitions) 
    %   sol:        solution evaluated at prices p (struct)
    %             (for detailed description fields in sol, see equilibrium.solve)
    %
    % OUTPUTS:
    %   grad_dp_dmp: (struct) total derivatives of change in model parameters (mp) taking into account the effect parameter changes has on equilirbrium prices
    %   grad_dp:     (struct) partial derivatives of change in princes for a given set of parameters    
    %   grad_dmp     (struct) partial derivatives of change in model parameters (mp) holding fixed prices
    %
    %   fields in grad_dp_dmp includes:
    %           'q','ccp_tau','ccp_scrap_tau','q_tau','h_tau' where as usual fields ending with "_tau" are mp.ntypes dimensional cell arrays 
    %           Dimensions of the fields of grad_dp_dmp is similar to sol, but with an additional dimension (the last one)
    %           that equals the number of parameters (np). Some examples below       
    %
    %             q:    [s.ns x np double]
    %       ccp_tau{t}: [s.ns x s.nd x np double]  
    % ccp_scrap_tau{t}: [s.ns x np double]  
    %         q_tau{t}: [s.ns x np double]  
    %         h_tau{t}: [s.ns x np double] 
    %
    % where np is the total number of parameters stored in mp.pvec_info.np
    %
    %   
    %   fields in grad_dp and grad_dmp includes:
    %           'u', 'ccp_scrap','ev', 'v', 'logccp', 'ccp', 'ctp', 'q', 'ed' (see description in g_trmodel.partial_deriv)
    %            (see g_trmodel.partial_deriv for description of the individual fields)
    %
    %  See also: g_trmodel.partial_deriv 

    if nargin==1
      s=trmodel.index(mp);
    end
    if nargin<3
      [sol]=equilibrium.solve(mp);
    end

    % analytical partial derivatives wrt prices
    grad_dp=g_trmodel.partial_deriv(mp, s, sol, 'prices');
    
    % analytical derivatives wrt model parameters
    grad_dmp=g_trmodel.partial_deriv(mp, s, sol, 'model');

    % derivatives of equilibrium prices wrt model parameters
    grad_dmp.p          =   -inv(grad_dp.ed)*grad_dmp.ed;

    grad_dp_dmp.q=0;

    for tau=1:mp.ntypes
        nmp=size(grad_dmp.u_tau{tau},3); 
        % derivative of conditional choice probabilities, ccp_tau
        dccp_dp_reshaped=reshape(grad_dp.ccp_tau{tau}, s.ns*s.nd, s.np);
        grad_dp_dmp.ccp_tau{tau}=reshape(dccp_dp_reshaped*grad_dmp.p,[s.ns,s.nd,nmp]) +grad_dmp.ccp_tau{tau};

        % derivative of scrap probabilities, ccp_scrap_tau
        grad_dp_dmp.ccp_scrap_tau{tau}=grad_dp.ccp_scrap_tau{tau}*grad_dmp.p +grad_dmp.ccp_scrap_tau{tau};

        % derivative of stationary_distribution, q_tau
        grad_dp_dmp.q_tau{tau}=grad_dp.q_tau{tau}*grad_dmp.p + grad_dmp.q_tau{tau};
        grad_dp_dmp.q=grad_dp_dmp.q+mp.tw(tau)*grad_dp_dmp.q_tau{tau}; 

        % derivative of stationary derivative of distribution of ownership after trading, h_tau
        grad_dp_dmp.h_tau{tau} = g_trmodel.dh_dp(mp, s, sol.q_tau{tau}, sol.delta_tau{tau}, grad_dp_dmp.ccp_tau{tau}, grad_dp_dmp.q_tau{tau}, tau);

    end
  end

  function [grad]=partial_deriv(mp, s, sol, parameter)
    % g_trmodel.partial_deriv: compute partial derivatives of a number of elements in solutions struct, sol
    %                          
    % SYNTAX 1: [grad] = g_trmodel.partial_deriv(mp, s, sol, parameter)
    %
    % 
    % INPUTS: 
    %   mp          structure with model parameters (see trmodel.setparams for definitions)
    %   s:          structure with indexes for model (see s=trmodel.index for definitions) 
    %   sol:        solution evaluated at prices p (struct)
    %               (for detailed description fields in sol, see equilibrium.solve)
    %   parameter:  (string) can be either 'prices' or 'model' 
    %               if parameter='prices': derivatives are taken wrt. prices of used cars 
    %               if parameter='model':  derivatives are taken wrt. model parameters defined in mp.pnames
    %               NOTE to use partial_deriv with parameter='model' you must first call [pvec0,mp] = estim.mp2pvec(mp,pnames)
    %               to populate mp with the field mp.pvec_info.
    %               
    % OUTPUT: 
    %   grad:       (struct) with partial derivatives of change in prices or model parameters (holding fixed prices in the latter case)
    %
    %   fields in grad includes: 'q','ccp_tau','ccp_scrap_tau','q_tau','h_tau' where as usual fields ending with "_tau" are mp.ntypes dimensional cell arrays 
    %   Some examples below       
    %
    %          u_tau{t}: [s.ns x s.nd x np double]
    %  ccp_scrap_tau{t}: [s.ns x np double] 
    %         ev_tau{t}: [s.ns x np double] 
    %          v_tau{t}: [s.ns x s.nd x np double] 
    %     logccp_tau{t}: [s.ns x s.nd x np double] 
    %        ccp_tau{t}: [s.ns x s.nd x np double] 
    %        ctp_tau{t}: [s.ns x s.ns x np double] 
    %          q_tau{t}: [s.ns x np double] 
    %             ed: [s.np x np double]
    %              p: [s.np x np double]
    %
    %   Dimensions of the fields of grad are similar to sol, but with an additional dimension (the last one)
    %   that equals the number of parameters (np). 
    %   If parameter = 'prices' np equals the number of used car prices (s.np)
    %   If parameter = 'model' np=mp.pvec_info.np
    %
    %  See also: g_trmodel.total_deriv and g_trmodel.gradients_tau

    % partial derivatives wrt parameter 
    ded=0;
    for tau=1:mp.ntypes;
      [du{tau}, dccp_scrap{tau}, dev{tau}, dv{tau}, dlogccp{tau}, dccp{tau}, dctp{tau}, dq{tau}, ded_tau, daccprob] ...
        = g_trmodel.gradients_tau(mp, s, tau, sol, parameter);
      ded=ded+mp.tw(tau)*ded_tau;
    end
    fields={'u_tau', 'accprob', 'ccp_scrap_tau','ev_tau', 'v_tau', 'logccp_tau', 'ccp_tau', 'ctp_tau', 'q_tau', 'ed'};
    grad=cell2struct({du, daccprob, dccp_scrap, dev, dv, dlogccp, dccp, dctp, dq, ded},fields, 2);
  end

  function [du, dccp_scrap, dev, dv, dlogccp, dccp, dctp, dq, ded, daccprob]=gradients_tau(mp, s, tau, sol, parameter)
    % gradients_tau: calculates the Jacobian matrix of the excess demand and contributions for consumer tau
    %                and all it's sub-derivatives wrt to changes in either prices (if parameter=='prices') or model parameters (if parameter=='model')
    %
    % SYNTAX: 
    %   [du, dccp_scrap, dev, dv, dlogccp, dccp, dctp, dq, ded]  = g_trmodel.gradients_tau(mp, s, tau, ccp, ccp_scrap, ctp, F, 'prices')
    %   [du, dccp_scrap, dev, dv, dlogccp, dccp, dctp, dq, ded]  = g_trmodel.gradients_tau(mp, s, tau, ccp, ccp_scrap, ctp, F, 'model')
    % 

    ev=sol.ev_tau{tau};
    ccp=sol.ccp_tau{tau};
    ccp_scrap=sol.ccp_scrap_tau{tau};
    ctp=sol.ctp_tau{tau};
    delta=sol.delta_tau{tau};
    F=sol.F;
    p=sol.p;
    
    if strcmp(parameter, 'prices')
      daccprob=[];
      parameter_types={'prices'};
      [du, dccp_scrap]=g_trmodel.du(mp, s, tau, ccp_scrap, p, 'prices');
      dbellman=g_trmodel.dbellman(mp, s, ev, ccp, ccp_scrap, ctp, delta, du, daccprob, 'prices'); 
    else
      parameter_types=unique(mp.pvec_info.pnames_type);

      daccprob=g_trmodel.daccprob(mp, s, parameter_types);
      np=mp.pvec_info.np;
      du=zeros(s.ns, s.nd, np); 
      dccp_scrap=sparse(s.ns, np); 
      dbellman=sparse(s.ns, np);
      for k=1: numel(parameter_types)
        idx=getfield(mp.pvec_info, sprintf('%s%s', parameter_types{k},'_idx'));
        [du(:,:,idx), dccp_scrap(:,idx)]=g_trmodel.du(mp, s, tau, ccp_scrap, p, parameter_types{k});
        dbellman(:,idx)=g_trmodel.dbellman(mp, s, ev, ccp, ccp_scrap, ctp, delta, du(:,:,idx), daccprob(:,idx), parameter_types{k}); 
      end
    end
    
    dev=(eye(s.ns)-mp.bet*ctp)\dbellman;
        
    % calculate the gradient of the choice specific value function
    dv = g_trmodel.dv(mp, s, F, dev, du, daccprob, ev, parameter_types);
    
    % calculate the gradient of the conditional choice probabilities
    dlogccp=g_trmodel.dlogccp_dp(mp, s, ccp, dv);

    % calculate the gradient of the conditional choice probabilities
    dccp=g_trmodel.dccp_dp(mp, s, ccp, dlogccp); 

    % % gradient of the  trade transition probability matrix
    % ddelta=g_trmodel.ddelta_dp(mp, s, dccp);

    % calculate the gradient of the controlled state transition probability matrix
    dctp=g_trmodel.dctp_dp(mp, s, F.notrade, dccp, delta, daccprob);

    % compute holdings distribution and gradient of holdings distribution
    [q, dq]=ergodic(ctp,dctp); 

    % compute the derivative of excess demand  for all relevant ages of used cars that can be purchased
    ded=g_trmodel.ded(mp, s, ccp, dccp, q, dq, ccp_scrap, dccp_scrap);
  end % end of ded_tau

  function [du, dccp_scrap]=du(mp, s, tau, ccp_scrap, p, parameter_type)
    if strcmp(parameter_type, 'prices') 
      [du, dccp_scrap]=g_trmodel.du_dp(mp, s, tau, ccp_scrap);
    elseif strcmp(parameter_type, 'theta')
      [du, dccp_scrap]=g_trmodel.du_dtheta(mp, s, tau, ccp_scrap, p);
    else % beta, alpha
      idx=getfield(mp.pvec_info, sprintf('%s%s', parameter_type,'_idx'));
      np=numel(idx);
      du=zeros(s.ns,s.nd,np);
      dccp_scrap=zeros(s.ns,np);
    end
  end

  function [du, dccp_scrap]=du_dtheta(mp, s, tau, ccp_scrap, p)
    %   du_dtheta:  Function to compute derivative of utility wrt to structural parameters
    %               for consumer tau in each state and for each decision 
    %
    %   SYNTAX:     [du, dccp_scrap]=g_trmodel.du_dtheta(mp, s, tau, ccp_scrap)
    %
    %   INPUTS:
    %     mp              structure with model parameters (see trmodel.setparams)
    %     s:              structure with indexes for dimensions. see trmodel.index
    %     ccp_scrap:      (s.ns x 1) vector of scrap probabilities. 
    %     tau:            (scalar interger) consumer type 
    %
    %   OUTPUTS          
    %   du:     (s.ns x s.nd x mp.nparam) with derivatives of utility wrt np x 1 parameter vector
    %           the last dimension is the price we take derivative wrt to
    %   dccp_scrap:     (s.ns x mp.nparam) with derivatives of ccp_scrap wrt np x 1 parameter vector
    %                   the last dimension is the price we take derivative wrt to
    %

    [pvec0, mp]=estim.mp2pvec(mp, mp.pnames);
    du=estim.matgrad_mp(@(mp) trmodel.utility(mp, s, tau, p), mp, pvec0,'theta');

    du(isnan(du))=0;
    du([s.is.clunker{:}] , s.id.keep,:)  =  0; % not possible to keep clunker
    du(s.is.nocar , s.id.keep,:)         =  0; % not possible to keep no car

    for i=1:mp.ncartypes;
      price_j{i}=p(s.ip{i});
    end;

    if mp.es
      dccp_scrap=estim.matgrad_mp(@(mp) trmodel.p_scrap(mp,s,price_j,tau), mp, pvec0,'theta');
    else
      dccp_scrap=zeros(s.ns, numel(pvec0)); 
    end
  end % end of du_dtheta

  function [du, dccp_scrap]=du_dp(mp, s, tau, ccp_scrap)
    %   du_dp:      Function to compute derivative of utility wrt to used car prices
    %               for consumer tau in each state and for each decision 
    %
    %   SYNTAX:     [du, dccp_scrap]=g_trmodel.du_dp(mp, s, tau, ccp_scrap)
    %
    %   INPUTS:
    %     mp              structure with model parameters (see trmodel.setparams)
    %     s:              structure with indexes for dimensions. see trmodel.index
    %     ccp_scrap:      (s.ns x 1) vector of scrap probabilities. 
    %     tau:            (scalar interger) consumer type 
    %
    %   OUTPUTS          
    %   du:     (s.ns x s.nd x s.np) with derivatives of utility wrt to prices. 
    %           the last dimension is the price we take derivative wrt to
    %           and is is indexed by s.ip
    %   ccp_scrap:  (s.ns x s.np) with derivatives of scrap probability wrt to prices. 
    %           the last dimension is the price we take derivative wrt to
    %           and is is indexed by s.ip

    dpbuy=zeros(s.nd, s.np); 
    dpbuy([s.id.trade_used{:}],:)= eye(s.np)*(1+mp.ptranscost);

    % derivative of net sales price after transactions cost 
    dpsell=zeros(s.ns,s.np); 
    dpsell([s.is.car_ex_clunker{:}],:)=eye(s.np)*(1-mp.ptc_sale);

    % derivative of value of scrap option
    dev_scrap =-ccp_scrap.*dpsell*mp.mum{tau}; 

    % derivative of scrap probability
    dccp_scrap=(1-ccp_scrap).*(dev_scrap)/mp.sigma_s;

    % derivative of utility
    du=mp.mum{tau}*(reshape(dpsell.*(1-ccp_scrap) ,[s.ns, 1, s.np])-reshape(dpbuy, [1, s.nd, s.np]));
    du(:,s.id.keep,:)=0; 
  end % end of utility

  function [daccprob]=daccprob(mp, s, parameter_types) 
    % daccprob: function to return the gradient of the accident probability with respect to accident parameters
    daccprob=sparse(s.ns,mp.pvec_info.np);
    if any(strcmp(parameter_types, 'alpha'))
      for j=1:mp.ncartypes;
          [~, daccprob(s.is.car{j},mp.pvec_info.alpha_idx)]=trmodel.acc_prob_j(mp,(0:s.abar_j{j}-1)',j, s.abar_j{j});
      end
    end
  end % end of daccprob

  function [dev]=dbellman(mp, s, ev, ccp, ccp_scrap, ctp, delta, du, daccprob, parameter_type) 
    % dev: function to return the gradient of the smoothed Bellman operator 
    % with respect to parameters types: 
    %   alpha (accident parameters), beta (discount factor), theta (utility parameters)

    if strcmp(parameter_type, 'beta') 
      dev=ctp*ev;
    elseif strcmp(parameter_type, 'alpha')
      dev=mp.bet*delta*(s.dQ*ev.*daccprob(s.tr.next_keep,:)); 
    else % theta, prices
      dev=permute(sum(ccp.*du, 2), [1,3,2]); % choice probability weighted sum of derivatives
    end
  end % end of dbellman

  function [dv]=dv(mp, s, F, dev, du, daccprob, ev, parameter_types)
    % dv   computes the gradients of the choice specific value functions 
    % 
    % OUTPUT
    %    dv:    s.ns x s.nd x s.np array of derivatives wrt to changes in value function (dv) 
    %           The first two dimensions the same as the util matrix (states and decisions)
 
    np=size(du, 3);
    nusedcars=numel([s.is.car_ex_clunker{:}]);
    dv=zeros(s.ns,s.nd,np); 

    % Step 1: take derivatives of F.notrade*ev and F.trade*ev (holding fixed F, varying ev)
    dev_keep=F.notrade*dev; 
    dev_trade_purge=F.trade*dev; 

    % Step 2: add derivative of F.notrade*ev and F.trade*ev (holding fixed ev, varying F)
    if numel(daccprob)>1
      dev_keep=dev_keep + (s.dQ*ev).*daccprob(s.tr.next_keep,:);
      dev_trade_purge=dev_trade_purge + (s.dF*ev).*daccprob(s.tr.next_trade,:);
    end 

    % gradient of v when keeping
    dv([s.is.car_ex_clunker{:}], s.id.keep,:) =du([s.is.car_ex_clunker{:}],s.id.keep,:) ... 
                           + reshape(mp.bet*dev_keep([s.is.car_ex_clunker{:}],:), [nusedcars, 1, np]);

    % gradient of v when trading and purging
    % we reshape to add a singleton choice-dimension to be able to broadcast gradient to all states
    % for both trade and purge, gradient is identical across states
    dv(:,[[s.id.trade{:}] s.id.purge],:) = du(:,[[s.id.trade{:}] s.id.purge],:) ...
                                 + reshape(mp.bet*dev_trade_purge(:,:), [1, s.ns, np]);

    if any(strcmp(parameter_types, 'beta')) 
      idx=getfield(mp.pvec_info, 'beta_idx');
      % dv(:,:, idx)=zeros(s.ns,s.nd,1); 
      dv([s.is.car_ex_clunker{:}], s.id.keep,idx)   =  dv([s.is.car_ex_clunker{:}], s.id.keep,idx)+  reshape(F.notrade([s.is.car_ex_clunker{:}],:)*ev, [nusedcars, 1, 1]); 
      dv(:, [[s.id.trade{:}] s.id.purge], idx)      =  dv(:, [[s.id.trade{:}] s.id.purge], idx) +  reshape(F.trade*ev, [1, s.ns, 1]);   
    end

  end % end of dvdp_updater

  function  dccp_dp=dccp_dp(mp, s, ccp, dlogccp);
    % SYNTAX: dccp_dp=g_trmodel.dccp_dp(mp, s, ccp, dvdp)
    % dccp_dp.m: updates the cell array dccp_dp_tau for the current type and price 
    dccp_dp=ccp.*dlogccp; 
  end % end of dccp_dp

  function  dlogccp=dlogccp_dp(mp, s, ccp, dvdp);
    % SYNTAX: dlogccp_dp=g_trmodel.dccp_dp(mp, s, ccp, dvdp)
    % dlogccp_dp.m: updates the cell array dlogccp_dp_tau for the current type and price 
    %
    % dlogccp_dp is computed as a choice probability weighted sum of derivatives (sum over choice alternatives)
    % assumes that infeasible choices have ccps and dvps = 0;
    dlogccp=(dvdp-sum(ccp.*dvdp, 2))/mp.sigma;
  end % end of dlogccp_dp


  function [ddelta_dp] = ddelta_dp(mp, s, dccp_dp)
    % ddelta_dp_dp: derivative of trade transition probability matrix with respect to used car prices. 
    % 
    % SYNTAX:     [ddelta_dp] = ddelta_dp(mp, s, dccp_dp)
    %             
    % OUTPUT: 
    %  ddelta_dp:   3 dimensional array (s.ns, s.ns, np) with derivatives of delta 
    %             wrt to change in np parameters such as pirces, p 
    % 
    % INPUT: 
    %   mp:       structure with model parameters (see trmodel.setparams)
    %   s:        structure with indexes for model (see trmodel.index)
    %   Q:        s.ns x s.ns extended physical transition matrix 
    %             (block-diagonal matrix with a block for each car and 1 for no car)
    % dccp_dp:    3 dimensional array (s.ns, s.nd, s.np) with derivatives of ccp wrt to pirces, p 
    % see also g_trmodel.dccp_dp

    % compute derivative of delta
    np=size(dccp_dp, 3);
    ddelta_dp=zeros(s.ns, s.ns, np); 
    for ip=1:np;
      ddelta_dp(:,:,ip)= dccp_dp(:,s.id.keep,ip)   + dccp_dp(s.tr.state,s.tr.choice,ip); 
    end 
  end

  function [dh] = dh_dp(mp, s, q, delta, dccp, dq, tau)
    % derivative of distribution of ownership after trading 
    % syntax: [dh] = g_trmodel.dh_dp(mp, s, q, delta, ddelta, dq, tau)
    % indexed by s.ipt

      np=size(dq,2);
      dqpt=zeros(s.ns, np);
      ddelta=zeros(s.ns, s.ns, np); 

      for ip=1:np;        
        ddelta(:,:,ip)=diag(dccp(:,s.id.keep,ip))   + dccp(s.tr.state,s.tr.choice,ip); 
        dqpt(:,ip)= ddelta(:,:,ip)'*q + delta'*dq(:,ip);
      end 

      for j=1:mp.ncartypes;
        dh(s.ipt.new{j},:)   =mp.tw(tau)*dqpt(s.is.clunker{j},:);
        dh(s.ipt.used{j},:)  =mp.tw(tau)*dqpt(s.is.car_ex_clunker{j},:);
        dh(s.ipt.nocar,:)    =mp.tw(tau)*dqpt(s.is.nocar,:);
      end
  end

  function [dmarketshares] = dmarketshares_dp(mp, s, dh_tau)
    % derivatives of market shares after trading 
    % syntax: [dmarketshares] = g_trmodel.dmarketshares_dp(mp, s, dh)
    % indexed by s.ipt
    np=size(dh_tau{1},2);
    dmarketshares=nan(mp.ntypes,mp.ncartypes+1, np);
    for tau=1:mp.ntypes;
      for j=1:mp.ncartypes;
        dmarketshares(tau,j,:)=sum(dh_tau{tau}(s.ipt.car{j},:),1);
      end
      dmarketshares(tau,mp.ncartypes+1,:)=dh_tau{tau}(s.ipt.nocar,:);
    end
  end

  function [dctp] = dctp(mp, s, Q, ddelta, delta, daccprob)
    % dctp_dp_dp: derivative of choice transition probability matrix with respect to parameters or used car prices. 
    % 
    % SYNTAX:     [dctp_dp] = dctp(mp, s, Q, dccp_dp, delta, daccprob)
    %             
    % OUTPUT: 
    %  dctp_dp:   3 dimensional array (s.ns, s.ns, s.np) with derivatives of cpt wrt to pirces, p 
    % 
    % INPUT: 
    %   mp:       structure with model parameters (see trmodel.setparams)
    %   s:        structure with indexes for model (see trmodel.index)
    %   Q:        s.ns x s.ns extended physical transition matrix 
    %             (block-diagonal matrix with a block for each car and 1 for no car)
    % dccp_dp:    3 dimensional array (s.ns, s.nd, s.np) with derivatives of ccp wrt to pirces, p 
    % see also g_trmodel.dccp_dp

    % compute derivative of deltaK*Q + derivative of deltaT*Q

    accprob=zeros(s.ns,1);
    for j=1:mp.ncartypes;
        accprob(s.is.car{j})=trmodel.acc_prob_j(mp,(0:s.abar_j{j}-1)',j, s.abar_j{j});
    end

    % F.notrade =s.Q.no_acc+s.dQ.*accprob(s.tr.next_keep);
    % s.Q.no_acc= sparse(s.tr.state,s.tr.next_keep,1,s.ns,s.ns);
    % s.Q.acc= sparse(s.tr.state,s.tr.next_acc,1,s.ns,s.ns);
    % dQ=s.Q.acc-s.Q.no_acc;

    dctp= ddelta(s.tr.state,s.tr.next_keep,:) + ddelta(s.tr.state,s.tr.next_acc,:).*reshape(accprob(s.tr.next_keep), [s.ns,1,1]);
    if numel(daccprob)>1
      for ip=1:np;
        dctp(:,:,ip)=dctp(:,:,ip) + delta*(s.dQ.*daccprob(s.tr.next_keep,ip)); 
      end 
    end
  end

  function [dctp_dp] = dctp_dp(mp, s, Q, dccp_dp, delta, daccprob)
    % dctp_dp_dp: derivative of choice transition probability matrix with respect to parameters or used car prices. 
    % 
    % SYNTAX:     [dctp_dp] = dctp_dp(mp, s, Q, dccp_dp)
    %             
    % OUTPUT: 
    %  dctp_dp:   3 dimensional array (s.ns, s.ns, s.np) with derivatives of cpt wrt to pirces, p 
    % 
    % INPUT: 
    %   mp:       structure with model parameters (see trmodel.setparams)
    %   s:        structure with indexes for model (see trmodel.index)
    %   Q:        s.ns x s.ns extended physical transition matrix 
    %             (block-diagonal matrix with a block for each car and 1 for no car)
    % dccp_dp:    3 dimensional array (s.ns, s.nd, s.np) with derivatives of ccp wrt to pirces, p 
    % see also g_trmodel.dccp_dp

    % compute derivative of deltaK*Q + derivative of deltaT*Q

    % accprob=zeros(s.ns,1);
    % for j=1:mp.ncartypes;
    %     accprob(s.is.car{j})=trmodel.acc_prob_j(mp,(0:s.abar_j{j}-1)',j, s.abar_j{j});
    % end

    % F.notrade =s.Q.no_acc+s.dQ.*accprob(s.tr.next_keep);

    np=size(dccp_dp, 3);
    dctp_dp=zeros(s.ns, s.ns, np); 
    for ip=1:np;
      dctp_dp(:,:,ip)= dccp_dp(:,s.id.keep,ip).*Q   + dccp_dp(s.tr.state,s.tr.choice,ip)*Q; 
    end 
    if numel(daccprob)>1
      for ip=1:np;
        dctp_dp(:,:,ip)=dctp_dp(:,:,ip) + delta*(s.dQ.*daccprob(s.tr.next_keep,ip)); 
      end 
    end
  end

  function dedvt=ded(mp, s, ccp, dccp, q, dq,ccp_scrap, dccp_scrap)
    % derivative of excess demand for all relevant ages of used cars that can be purchased for consumer t
    np=size(dccp, 3);
    dedvt=zeros(s.np,np);

    dDdp=sum(q.*dccp(:,[s.id.trade_used{:}],:),1);
    dSdp=(dccp([s.is.car_ex_clunker{:}],s.id.keep,:)).*(1-ccp_scrap([s.is.car_ex_clunker{:}],:)).*q([s.is.car_ex_clunker{:}]) ; % Add es here

    dedvt([s.ip{:}],:)= permute(dDdp, [2, 3, 1]) +permute(dSdp,[1, 3, 2]) ...    % 1st term: partial derivative of excess demand wrt to change in ccps
    +ccp(:,[s.id.trade_used{:}])'*dq ... % 2nd term can be done as a single matrix x matrix multiplication of ccp and dq 
    -(1-ccp([s.is.car_ex_clunker{:}],s.id.keep)).*(1-ccp_scrap([s.is.car_ex_clunker{:}],:)).*dq([s.is.car_ex_clunker{:}],:) ... 
    +(1-ccp([s.is.car_ex_clunker{:}], s.id.keep)).*dccp_scrap([s.is.car_ex_clunker{:}],:).*q([s.is.car_ex_clunker{:}]);

    if 0 % new version (identical)
      % ed=sum(ccp(s.tr.state,s.tr.choice(idx)).*q, 1)' - (1-ccp(idx, s.id.keep)).*(1-ccp_scrap(idx,:)).*q(idx);
      ded=zeros(s.np,np);
      idx=[s.is.car_ex_clunker{:}];
      ded([s.ip{:}],:)= permute(sum(...
                 q.*dccp(:,s.tr.choice(idx),:) + ccp(:,s.tr.choice(idx)).*permute(dq, [1, 3, 2]), ...
                  1), [2, 3, 1]) ...
              - (0-permute(dccp(idx, s.id.keep,:), [1, 3, 2])).*(1-ccp_scrap(idx,:))   .*q(idx) ...
              - (1-ccp(idx, s.id.keep))                       .*(0-dccp_scrap(idx,:))  .*q(idx) ...
              - (1-ccp(idx, s.id.keep))                       .*(1-ccp_scrap(idx,:))   .*dq(idx,:);

      % keyboard
      dedvt=ded; 
    end

  end % end of ded

end % end of methods
end % end of classdef
