%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for carmodel
% HOWTO:
% tests.run       Does standard battery of tests at default parameters
% tests.run(mp)   Does standard battery of tests at parameters mp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef tests
properties (Constant)
    output_derivative_check_tables = 'detailed'; % 'on', 'off', or 'detailed'   
   check_aggderivs = true
  end%properties
methods (Static)

  function [num_failed, fails] =run(mp, onlypartial)
    % SYNTAX: 
    %   Run tests at parameters mp: 
    %       [num_failed, fails] =tests.run(mp)
    %   Run tests at default parameters : 
    %         [num_failed, fails] =tests.run
    %
    %   OUTPUTS: 
    %     num_failed: number of tests fails given tolerances specified below (0=passed, 1=failed)
    %     fails:      struct with information about individual tests (0=passed, 1=failed)
    %
    %   INPUT: 
    %     mp (optional): parameter structure that index where to evaluate model
    %     
    %     mp.p0 and mp.pnames.
    %     To compare numerical to analytical derivatives mp.pnames and mp.p0 must be non-empty
    %     and must be mutually consistent. 
    %
    %     Example: Set parameters to default values 
    %              and compute derivatives of parameters in mp.u_0 and mp.mum
    %       mp=trmodel.setparam(); 
    %       mp.pnames={'u_0', 'mum'};
    %       mp.p0.u_0={60; 60}; 
    %       mp.p0.mum={1; 1.5};
    %       [num_failed, fails] =tests.run(mp); 
    %
    %     Here the total number of parameters is 4 since we allow mp.u_0 and mp.mum to be 
    %     heterogeneous across consumers, but homogeneous across cars. 
    %     Yet, we evaluate derivatives in the same value 60 for mp.u_0, but at distinct values for mum 
    %  
    %     Here mp.pnames={'u_0', 'mum'} instructs tests.run(mp) to perform comparison 
    %     of derivatives for parameters mp.u_0 and mp.mum and the user must provide 
    %     the structure corresponding elements in mp.p0. 
    %
    %     The elements in mp.p0 both specifies parameter values where to compute derivatives and 
    %     the structure of how they enter in the specification. 
    %
    %     For example 
    %               mp.p0.u_0={60}; % mp.u_0 is homogeneous across consumer and car types
    %
    %               mp.p0.u_0={60; 60}; % mp.u_0 is heterogeneous across the 2 (first) consumer types, 
    %                                     but homogeneous across car types
    %
    %               mp.p0.u_0={60,60; 60,60}; % mp.u_0 is homogeneous across 2 (first) consumer and car types  
    

    if nargin<1
      % standard model parameters
      mp=trmodel.setparam(); 
    end

    if nargin<2
      mp.onlypartial=false; 
    else
      mp.onlypartial=onlypartial;
    end


    fails =tests.deriv(mp);

    % add more tests here 

    % print summary output
    fprintf('\n%s\n',repmat('-', 100,1));
    disp('Test summary (0=passed, 1=failed)')
    fprintf('%s\n',repmat('-', 100,1));

    fld=fields(fails);
    for i=1:numel(fld)
      any_fails.(fld{i})=max(fails.(fld{i}));
    end
    disp(struct2table(any_fails));
    num_failed=sum(struct2array(any_fails));
    if num_failed==0
      fprintf("All tests PASSED\n");
    else
      fprintf("%Summary: d test FAILED\n", num_failed);
    end
    fprintf('%s\n',repmat('-', 100,1));
  end

  function [fails] = deriv(mp)

    % elements in solution to consider (for comparison of partial derivatives)
    fcts={'u_tau','ev_tau','ccp_tau','ctp_tau', 'ccp_scrap_tau', 'q_tau'};
    lbl={'Utility function','Expected value','Conditional choice probability','State transition probability', 'Scrap probability', 'Car quantity'};
    
    % elements in solution to consider (for comparison of total derivatives)
    fcts_eq={'ccp_tau','ccp_scrap_tau', 'q_tau', 'h_tau'};
    lbl_eq={'Conditional choice probability', 'Scrap probability', 'Car quantity', 'Post trade ownership distribution'};

    % create parameter vector and update mp with values in mp.p0
    [pvec0, mp]=estim.mp2pvec(mp, mp.pnames);
    mp = estim.update_mp(pvec0, mp);

    s=trmodel.index(mp);

    % solve model at mp
   [sol]=equilibrium.solve(mp);

   % analytical total and partial derivatives 
   [grad_dp_dmp, grad_dp, grad_dmp]=g_trmodel.total_deriv(mp, s, sol);

   % Compare with numerical derivatives
   if ~strcmp(tests.output_derivative_check_tables, 'off')
      fprintf('\n%s\n',repmat('-', 100,1));
      fprintf('Partial derivatives with respect to model parameters, mp\n')
      fprintf('%s\n',repmat('-', 100,1));
   end

   for tau=1:mp.ntypes;
      if ~strcmp(tests.output_derivative_check_tables, 'off')
        fprintf('\nConsumer type %d\n', tau);
        fprintf('%s\n',repmat('-', 50,1));
      end

      fx=      @ (mp, fct)  tests.eval_fct_fixp(mp, s, sol.p, sol.ev_tau{tau}, fct, tau);
      
      for i=1:numel(fcts)
        % partial derivatives of fcts wrt model parameters
        numgrad_dmp.(fcts{i}){tau}   =   estim.matgrad_mp(@(mp) fx(mp, fcts{i}),mp,pvec0,'all');
        % print output
        fails.((sprintf('d%s_dmp', fcts{i})))(tau)=tests.deriv_compare(mp, grad_dmp.(fcts{i}){tau}  , numgrad_dmp.(fcts{i}){tau},  sprintf('%s, %s{%d}', lbl{i}, fcts{i}, tau));
      end

   end

   if (tests.check_aggderivs)
    if ~strcmp(tests.output_derivative_check_tables, 'off')
      fprintf('\nAggregate quantities\n');
      fprintf('%s\n',repmat('-', 50,1));
    end
    
    numgrad_dmp.ed       =   estim.matgrad_mp(@(mp) fx(mp, 'ed'),mp,pvec0,'all');
    fails.ed_dmp=tests.deriv_compare(mp, grad_dmp.ed, numgrad_dmp.ed,   sprintf('Excess demand, ed'));

    if mp.onlypartial
      return
    end

    if ~strcmp(tests.output_derivative_check_tables, 'off')
      fprintf('\n%s\n',repmat('-', 100,1));
      fprintf('Derivatives of prices, p, with respect to model parameters, mp\n')
      fprintf('%s\n',repmat('-', 100,1));
    end

    fx_equil=@ (mp, fct)  tests.eval_fct_equil(mp, fct, nan);
    numgrad_dmp.p             =   estim.matgrad_mp(@(pvec) fx_equil(pvec, 'p'),mp,pvec0,'all');
    fails.dp_dmp=tests.deriv_compare(mp, grad_dmp.p, numgrad_dmp.p, 'Prices, p');

    if ~strcmp(tests.output_derivative_check_tables, 'off')
      fprintf('\n%s\n',repmat('-', 100,1));
      fprintf('Total derivatives with respect to model parameters, mp, and their effect though prices, p\n')
      fprintf('%s\n',repmat('-', 100,1));
    end

    for tau=1:mp.ntypes;
      if ~strcmp(tests.output_derivative_check_tables, 'off')
        fprintf('\nConsumer type %d\n', tau);
        fprintf('%s\n',repmat('-', 50,1));
      end

      fx_equil=@ (mp, fct)  tests.eval_fct_equil(mp, fct, tau);
      for i=1:numel(fcts_eq)
        % partial derivatives of fcts_eq wrt model parameters
        numgrad_dp_dmp.(fcts_eq{i}){tau}  =    estim.matgrad_mp(@(mp) fx_equil(mp, fcts_eq{i}),mp,pvec0,'all');

        % print output
        fails.(sprintf('d%s_dp_dmp', fcts_eq{i}))(tau)=tests.deriv_compare(mp, grad_dp_dmp.(fcts_eq{i}){tau}  , numgrad_dp_dmp.(fcts_eq{i}){tau},  sprintf('%s, %s{%d}', lbl_eq{i}, fcts_eq{i}, tau));
      end

    end
   end
  end


  function [out] = eval_fct_fixp(mp, s, p, ev0, fct, tau)
      % evaluates outcomes at fixed price vector, p for a given consumer type tau
      % used for computing numerical derivatives
      
      % compute utility and scrap probs
      [u, ev_scrap, ccp_scrap] = trmodel.utility(mp, s, tau, p);   % update the current period, post-trade payoffs stored in the global util_tau

      if strcmp(fct, 'u_tau')
        out=u;
        return
      elseif strcmp(fct, 'ccp_scrap_tau')
        out=ccp_scrap;
        return
      end

      % update transition matrices
      [F] = trmodel.age_transition(mp, s);

      % bellman equation
      bellman=@(ev) trmodel.bellman(mp, s, u, F, ev);

      if strcmp(fct, 'bellman_tau')
        out=bellman(ev0);
        return
      end

      % solve bellman equation      
      [ev, ccp]= dpsolver.poly(bellman, ev0, [], mp.bet);  
      [ev, ccp]=bellman(ev); % also return value and policy function   

      if strcmp(fct, 'logccp_tau')
        out=log(ccp);
        return
      elseif strcmp(fct, 'ev_tau')
        out=ev;
        return
      elseif strcmp(fct, 'ccp_tau')
        out=ccp;
        return
      elseif strcmp(fct, 'delta_tau')
        [delta]=trmodel.trade_transition(mp, s, ccp);    
        out=delta;
        return
      elseif strcmp(fct, 'ctp_tau')
        [delta]=trmodel.trade_transition(mp, s, ccp);    
        ctp=delta*F.notrade; 
        out=ctp;
        return
      elseif strcmp(fct, 'q_tau')
        [delta]=trmodel.trade_transition(mp, s, ccp);    
        ctp=delta*F.notrade; 
        q=ergodic(ctp); 
        out=q;
        return
      elseif strcmp(fct, 'ed')
         [ed]=equilibrium.edf(mp, s, p);
        out=ed;
        return
      end

      error('unknown fct');

    end

    function [out] = eval_fct_equil(mp, fct, tau)

      % solve for equilibrium 
      [sol]=equilibrium.solve(mp);

      if any(strcmp(fct, {'q_tau','ev_tau','ccp_tau','ccp_scrap_tau','ctp_tau','delta_tau','h_tau'}))
        out=sol.(fct){tau};
      elseif strcmp(fct, 'logccp')
        out=log(sol.ccp_tau{tau});
      elseif any(strcmp(fct, {'p', 'q','ed', 'F'}))
        out=sol.(fct);
      else 
        error('unknown fct')
      end
    end


  function [exitflag] = deriv_compare(p, g, g_num, header)
    % Compute and compare analytical vs numerical derivatives

    function lbl=f_param_lbl(p,pnames)
      if nargin<2; 
        pnames=p.pnames;
      end
      ij=1;
      for i=1:numel(pnames);
        for j=1:numel(p.p0.(char(pnames(i))));
          lbl{ij} = sprintf('%s(%i)',char(pnames(i)),j);
          ij=ij+1;
        end;
      end; 
    end
 
    if (numel(g) ~= numel(g_num))
        exitflag=0;
        fprintf('Error: tests.deriv_compare:  sizes of g and g_num arrays are different size(g)=%i size(g_num)=%i\n',numel(g),numel(g_num));
        return;
    end

    theta=struct2vec(p.p0, p.pnames);

    sz=size(g);
    g=reshape(g, prod(sz)/numel(theta), numel(theta));
    g_num=reshape(g_num, prod(sz)/numel(theta), numel(theta));

    % generate output
    out_mat=[theta max(abs(g))' max(abs(g_num))'  max(abs(g-g_num))' max(abs((g+1e-12)./(g_num+1e-12)-1))'];
    VariableNames={'param', '||analytical||','||numerical||','||an-num||','||an/num-1||'}; 
    param_lbl=f_param_lbl(p)';

    MAE=max(max(abs(g-g_num))');
    MAE_rel= max(max(abs((g+1e-12)./(g_num+1e-12)-1))'); 

    tol=1e-3;
    rtol=1e-4;
    exitflag=(MAE>tol).*(MAE_rel>rtol); 

    out_tab=array2table(out_mat, 'RowNames',param_lbl, 'VariableNames', VariableNames);

    % OUTPUT
    if ~strcmp(tests.output_derivative_check_tables, 'off')
        fprintf('\n### %s  ###\n', header); 

        if strcmp(tests.output_derivative_check_tables, 'detailed') || exitflag==1
            disp(out_tab)
        end

        fprintf('Max supnorm(err) over all parameters.    : %16g\n', MAE);
        fprintf('Max supnorm(rel_err) over all parameters : %16g\n', MAE_rel);
        if exitflag==0
            disp('Derivative check PASSED');  
        else            
            disp('Derivative check FAILED');

            fprintf('\nParameters where ||g-g_num||>%g\n', tol);
            disp(out_tab((out_tab{:,'||an-num||'}>tol),:));  

            fprintf('\nParameters where ||(g+esp)/(g_num+eps)-1||>%g\n', rtol);
            disp(out_tab((out_tab{:,'||an/num-1||'}>rtol),:));
        end 
        fprintf('### END %s  ###\n', header); 
    end
  end

end%methods    
end%class
