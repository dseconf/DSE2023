%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specification of cost, profits etc. for car producers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef firms
properties (Constant)
    debg = 'on'; % 'on', 'off', or 'detailed'   
  end%properties
methods (Static)

        function [pvec1,resids,exitflag]=primary_market_equilibrium(mp,pvec0,solution_method);

           nfirms=max(mp.firm);
           cars=(1:mp.ncartypes)';

           nfirm_check=unique(mp.firm);
           if (size(nfirm_check,2) ~= nfirms)
              fprintf('Error in indexing of firms: check mp.firm to make sure there are no gaps in the firm ownership indexing\n');
              mp.firm
              pvec1=pvec0;
              resids=NaN;
              exitflag=0;
              return;
           end;

           if (nfirms == 1)
                fprintf('Only 1 firm in the market: computing the monopoly solution\n');
		opt = optimoptions('fminunc','Display','off','Algorithm','trust-region','SpecifyObjectiveGradient',true);
		[pvec1,mprofit,exitflag]=fminunc(@(pvec) firms.profits(mp,pvec),pvec0,opt);
                [mprofit,resids]=firms.profits(mp,pvec1);
                resids=resids';
                return;
           end;
 

           switch solution_method

             case 'brz'

                disp('find Bertrand-Nash oligopoly prices of new cars in primary market by using fsolve to set best response residuals to 0');
                tic
                opt = optimoptions('fsolve','Display','iter','Algorithm','trust-region','SpecifyObjectiveGradient',false);
                [pvec1,resids,exitflag]=fsolve(@(pvec) firms.br_diff(mp, pvec),pvec0,opt);
                toc

                if (exitflag < 0)
                fprintf('Equilibrium not found: fsolve exitflag: %i\n',exitflag);
                end

             case 'sa'

                disp('find Bertrand-Nash oligopoly prices of new cars in primary market by successive approximations on the best response function');
                tic
                cars=(1:mp.ncartypes)';
                for i=1:100
	          pvec1=firms.br(mp, pvec0);
	          fprintf('iteration %d\n', i);
	          err=pvec1-pvec0;
	          disp(array2table([mp.firm(cars)' cars pvec0 pvec1 err], 'VariableNames', {'firm','car type', 'pvec0', 'pvec1','err'}))
	          pvec0=pvec1;
	          if max(abs(err))<1e-6
		     fprintf('convergence achieved after %d iterations\n', i);
		     break
	          end
                end
                resids=firms.foc(mp,pvec0);
                if (i == 100)
                   exitflag=0;
                   fprintf('Equilibrium not found: successive approximations terminated after %i iterations but did not converge: max absolute change in prices is %g\n',i,max(abs(err)));
                else
                   exitflag=1;
                end
                toc

             case 'foc'

                disp('find Bertrand-Nash oligopoly prices of new cars in primary market using fsolve to find a zero of the firm first order conditions');
                tic
                opt = optimoptions('fsolve','Display','iter','Algorithm','trust-region','SpecifyObjectiveGradient',false,'FunctionTolerance',1e-7);
                [pvec1,resids,exitflag]=fsolve(@(pvec) firms.foc(mp, pvec),pvec0,opt);
                toc

                if (exitflag < 0)
                fprintf('Equilibrium not found: fsolve exitflag: %i\n',exitflag);
                end

             case 'focsa'   % start using the first order conditions approach then switch to successive approximations

                disp('find Bertrand-Nash oligopoly prices of new cars in primary market using fsolve to find a zero of the firm first order conditions, then switch to successive approximations');
                tic
                opt = optimoptions('fsolve','Display','iter','Algorithm','trust-region','SpecifyObjectiveGradient',false,'FunctionTolerance',1e-7);
                [pvec1,resids,exitflag]=fsolve(@(pvec) firms.foc(mp, pvec),pvec0,opt);
                toc
                pvec0=pvec1;
                tic
                cars=(1:mp.ncartypes)';
                for i=1:100
	          pvec1=firms.br(mp, pvec0);
	          fprintf('iteration %d\n', i);
	          err=pvec1-pvec0;
	          disp(array2table([mp.firm(cars)' cars pvec0 pvec1 err], 'VariableNames', {'firm','car type', 'pvec0', 'pvec1','err'}))
	          pvec0=pvec1;
	          if max(abs(err))<1e-6
		     fprintf('convergence achieved after %d iterations\n', i);
		     break
	          end
                  resids=firms.foc(mp,pvec0);
                end
                if (i == 100)
                   exitflag=0;
                   fprintf('Equilibrium not found: successive approximations terminated after %i iterations but did not converge: max absolute change in prices is %g\n',i,max(abs(err)));
                end
                toc
     
             otherwise

                fprintf('unrecognized solution method  %s was supplied: please use either sa foc brz or focsa\n',solution_method);

           end

        end

	function [pvec1]=br(mp, pvec0);
		opt = optimoptions('fminunc','Display','off','Algorithm','trust-region','SpecifyObjectiveGradient',true);
		%opt.MaxIterations=1;
		mp=estim.pvec2param(pvec0, mp);
		nfirms=max(mp.firm);
                distinct_firms=size(unique(mp.firm),2);
                if (distinct_firms ~= nfirms)
                  fprintf('Error br: gaps in the numbering of firms in mp.firm:  largest firm is %i but there are only %i distinct firms in mp.firm\n',nfirms,distinct_firms);
                  return;
                end
		pvec1=pvec0; 
		for i=1:nfirms
			[price_i, profits_i, exitflag]=fminunc(@(pvec) firms.profits_i(mp, pvec, i), pvec0(mp.firm==i), opt);
			pvec1(mp.firm==i)=price_i; 
		end 

	end

        function [diff]=br_diff(mp,p)

                  diff=p-firms.br(mp,p);

        end

	function [profits_i, dprofits_i]=profits_all(mp, p_all);
		% profits of all firms
		[profits, dprofits]=firms.profits(mp, p_all);
		profits_i=sum(profits); 		
		dprofits_i=sum(dprofits,1); 
	end

	function [profits_i, dprofits_i]=profits_i(mp, p_i, firm_idx);
		% profits of firm if and drivatives wrt to the prices controlled by firm i
		p_all=struct2vec(mp,{'pnew_notax'}); 
		p_all(mp.firm==firm_idx)=p_i;
		[profits, dprofits]=firms.profits(mp, p_all);
		profits_i=profits(firm_idx); 		
		dprofits_i=dprofits(firm_idx,find(mp.firm==firm_idx)); 

	end

        function [foc]=foc(mp,p);

          % computes a vector of first order conditions for profit maximization for all products sold in the market

          foc=zeros(mp.ncartypes,1);

          [profits,dprofits]=firms.profits(mp, p);

	  nfirms=size(profits,1);

          for i=1:mp.ncartypes

              foc(i)=dprofits(mp.firm(i),i);

          end

        end

	function [profits, dprofits]=profits(mp, p);

        % computes profits for all firms allowing for any firm to sell multiple car types
        % returns profits, a vector of dimension mp.firm (number of firms) x 1 where the ith element is total profits
        %                  earned by firm i for all of the car types it sells/produces
        %         dprofits, a matrix of dimension mp.firm x mp.ncartypes  with the gradients of profits for each firm
        %                   for all mp.ncartypes cars produced in the market
     

		nfirms=max(mp.firm);

		% update prices
		mp.pnames = {'pnew_notax'}; 
                mp=estim.pvec2param(p, mp);
                [mp.p_fuel, mp.pnew, mp.pscrap] = trmodel.price_aftertax(mp);

		% solve for equilibrium at secondary market
		global price_j0;   % this global remembers a recent previous equilibrium solution in the hopes it is near the equilibrium for this solution
                [s,F,price_j,q_tau,q,ev_tau,ccp_tau,ccp_scrap_tau,ctp_tau,ed,h_tau,marketshares,exitflag]=equilibrium.solve(mp, price_j0);
                if (exitflag < 0)
                  fprintf('equilibrium failed to solve after 2 tries from starting prices. Terminating search for optimal price for firms in primary market.\n');
                   return;
                end
		price_j0=price_j;

		profits=zeros(nfirms,1);
                demand=zeros(mp.ncartypes,1);

		for j=1:mp.ncartypes
			% loop over consumer types to calculate demand for new cars weighted by q_tau, the pre-trade holdings of cars by consumers of type tau
			for tau=1:mp.ntypes
			   demand_nc_j=ccp_tau{tau}(:,s.id.trade_new{j})'*q_tau{tau}*mp.tw(tau);          
			   demand(j)=demand(j)+demand_nc_j;
			end
			i=mp.firm(j);
			profits(i)=profits(i)+([mp.pnew_notax{j}]-[mp.mc{j}]).*demand(j);
			%fprintf('profits for car %i %g  demand %g\n',j,([mp.pnew_notax{j}]-[mp.mc{j}]).*demand(j),demand(j));
		end

                profits= - profits;

                % compute derivatives of profits and demand with respect to prices
		if nargout>1
	   	  [grad_dp_dmp,grad_dp,grad_dmp]=g_trmodel.total_deriv(mp,s,price_j,ccp_tau,ccp_scrap_tau,ctp_tau,Q,F);

	    	  ddemand=zeros(mp.ncartypes,numel(p));
		  dprofits=zeros(nfirms,numel(p));
		  for j=1:mp.ncartypes
                    for tau=1:mp.ntypes
        	       ddemand_nc_j=ccp_tau{tau}(:,[s.id.trade_new{j}])'*grad_dp_dmp.q_tau{tau};
                       for ip=1:numel(p)
                          ddemand_nc_j(ip)=ddemand_nc_j(ip)+grad_dp_dmp.ccp_tau{tau}(:,[s.id.trade_new{j}],ip)'*q_tau{tau};              							 
                       end
	               ddemand(j,:)=ddemand(j,:)+ddemand_nc_j*mp.tw(tau);
                    end

                    % derivative of firm profits
	      	    i=mp.firm(j);
	  	    dprofits(i,:)=dprofits(i,:)+ (mp.pnew_notax{j}-mp.mc{j})*ddemand(j,:);
	   	    dprofits(i,j)=dprofits(i,j) + demand(j);
		  end
		  dprofits=- dprofits;
                end

   end  % end of profits function


end%methods    
end%class
