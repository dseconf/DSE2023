% The main static class for the homogeneous consumer economy and social planners problem.
% 
% July 2020
%

classdef spp
properties (Constant)
    verbosity=0;
    abar_max=100;
    solmeth={'poly'} %  solution method for spp DP problem
      % 'poly' (fastere) for poly algorithm that starts with successive approximations
      % 'policy' pure policy iteration

  % Constant properties social planners problmes
end

methods (Static)
  function [abar_spp, p_spp, w_spp]=solve_spp(mp);
    % solve_spp    solves the social planning problem to determine the optimal scrapping
    %              threshold by either successive approximations of policy iteration, depending
    %              on choice of solution method solmeth switch below.
    % INPUT 
    %   mp:   (optional) Structure with model parameters (see trmodel.setparam)
    %         if called without input, defaul parameters from trmodel.setparam is used. 
    %
    % OUTPUT
    %  W__spp_tau_j, abar_spp, and p_spp are all mp.ntypes x mp.ncartypes cell arrays with elements: 
    %  (where tau and j refer to consumer type and car type respectively) 
    %
    %   abar_spp{tau, j}  : optimal scrap dates for consumer types tau and car types j 
    %   p_spp{tau, j}     : (Shadow) Price function of social planner (abar+1 x 1). 
    %   W_spp{tau, j}     : Value function of social planner (abar+1 x 1) 
    %                       p_spp(i, j) includes new price and scrap price as first and last element   
    %
    %  John Rust, Georgetown University January 2019

    mp=trmodel.update_mp(mp); % update dependent parameters  

    abar=spp.abar_max;  % pick a large upper bound for the possible scrap age of a car
    abar_spp = cell(mp.ntypes, mp.ncartypes);
    
    for car=1:mp.ncartypes
      for t=1:mp.ntypes
        % abar_spp{t,car}=abar;
        %fprintf('Solving the social planning problem and computing homogeneous agent, zero transactions cost equilibrium for type %i consumers\n',t);
        
        v0=zeros(abar,1);

        switch spp.solmeth{1}
          case 'poly'   % poly algorithm - start with successive approximations        
            [w, dr]= dpsolver.poly(@(V) spp.bellman_spp(mp, V, t, car), v0, mp, mp.bet);
          case 'policy' % policy iteration solution
            [w, dr]= dpsolver.policy(@(V) spp.bellman_spp(mp, V, t, car), v0, mp);
          otherwise 
            error('spp.solmeth, must be eithr "poly" or "policy".'); 
        end 

        abar_spp{t,car}=min(find(dr==1))-1; %-1 because social planner has new cars in state space
        if ~any(dr(1:end-1)==1) 
          abar_spp{t,car}=mp.abar_j0{car}; 
          disp('Not optimal to replace at any age in grid choosing abar_spp{t,car}=mp.abar_j0{car} -> increase spp.abar_max)');
        end

        w=w(1:abar_spp{t,car}+1);
        
        p_spp{t, car}=mp.pnew{car}-(w(1)-w)/mp.mum{t};
        w_spp{t, car}=w;
        if spp.verbosity>0
          fprintf('solved spp for cartype %i consumer type=%i abar_spp=%g\n',car,t,abar_spp{t,car});
        end

      end % end of loop over consumer types
    end % end of loop over car types
  end % end of solve_spp

  function [V, P, dV, Eu]=bellman_spp(mp, V, t, car)
    % bellman_spp: Bellman equation for social planner
    % Input 
    %   V:    Initial guess on value function 
    %         row vector with abar+1 elements V(0),..., V(abar)
    %   mp:   Structure with model parameters 
    %   t:   index for consumer (t<=mp.ntypes)
    %   car:  index for car type (car<=mp.ncartypes)
    %
    % Output
    %   V:  Updated value after one bellman iteration (abar + 1 x 1)

    % Indexing: Values associated with a is stored in V with indexes a+1 
    % (due to base 1 indexing in Matlab)      
    abar=numel(V)-1;
    a=(0:abar-1)';    % ages, a=(0, 1, ,.. abar)
    ia=a+1;

    %accident probabilities and utility
    apr=trmodel.acc_prob_j(mp, (0:abar)',car,abar);
    urep = trmodel.u_car(mp, 0 ,t, car)- mp.mum{t}*(mp.pnew{car}-mp.pscrap{car});
    ukeep = trmodel.u_car(mp, a,t, car);
    vrep= urep +  mp.bet*( (1-apr(1))*V(2) + apr(1)*V(end) );

    vkeep = ukeep  +   mp.bet*((1-apr(ia)).*V(2:end) + apr(ia).*V(end));

    V(1:end-1)=max(vrep, vkeep);
    V(end)= vrep;

    if nargout>1  % Policy function (indicator for replacement)
      P=[sparse(vrep>vkeep);1];  % always replace oldest car
    end

    if nargout>2  %Frechet derivative of bellman operator (needed for dpsolver.policy and dpsolver.poly)
      tpm=sparse(diag((1-P(ia)).*(1-apr(ia)), 1));     % next period car age one year older if not replaced and not in accident
      tpm(:,end)=tpm(:,end)+P.*apr(1) + (1-P).*apr ;   % next period car is clunker if replaced and in accident or not replace car and in accident
      tpm(:,2)=tpm(:,2)+P.*(1-apr(1));                 % next period car is one year old if replaced znd not in accident
      dV = mp.bet.*tpm;
    end

    if nargout>3  % expect utility given policy rule (needed for dpsolver.policy)
      Eu=P.*urep + (1-P).*[ukeep; urep];
    end
  end % end bellman_spp

  function [P, a] = hom_prices(mp, abar, tp, car);
    % This function calculates equilibrium prices in the homogeneous model using the closed-form solution.
    % Syntax: 
    %   [P, a] = prices(mp, abar, tp, car) if abar is a scalar
    %   [P, a] = prices(mp, abar) if abar is a cell array
    %
    % Input 
    %   mp:   Structure with model parameters 
    %   abar: Scalar or mp.ntypes x mp.ncartypes cell arrays scrap date for each consumer type and car type 
    %         (scalar on only if price for a single type is computed and optional inputs tp and car is supplied)
    %   tp:   (optional input) index for consumer (tp<=mp.ntypes)
    %   car:  (optional input) index for car type (car<=mp.ncartypes)
    %
    % Output
    %   P     If called with two inputs, hom_prices(mp, abar), then P is a mp.ntypes x mp.ncartypes 
    %         cell array that hold equilibrium price vectors (abar{tau, j)-1 x 1) for each consumer tau and car j.   
    %         NOTE: The price vector P{tau, j} does NOT new price and scrap price. 
    %
    %         If called with four inputs (mp, abar, tau, j), then P is a price vector for the consumer, car combination tau, j
    %         Prices are found as as linear system X*P=Y with P=(P(1),...,P(abar-1)) and X and Y defined as in the paper

    if nargin ==2; 
      abar_cell=abar;
      CARS=1:mp.ncartypes; 
      CONSUMERS=1:mp.ntypes;
    elseif nargin ==4;
      CARS=car; 
      CONSUMERS=tp; 
      abar_cell{tp, car}=abar;
    end 
      
    for car=CARS; 
      pnew=mp.pnew{car};
      pscrap=mp.pscrap{car};

      for tau=CONSUMERS;
        abar=abar_cell{tau, car};  
        mum=mp.mum{tau};

        x=diag((trmodel.acc_prob_j(mp, (0:abar-2)',car,abar)-1)*mp.bet-1) ...
        +diag(ones(1,abar-2),-1) ...
        +diag((1-trmodel.acc_prob_j(mp, (1:abar-2)',car,abar))*mp.bet,1);

        y=trmodel.u_car(mp, [0:abar-2]',tau,car) / mum ...
        -trmodel.u_car(mp, [1:abar-1]',tau,car) / mum ...
        +mp.bet*pscrap*[trmodel.acc_prob_j(mp, (0:abar-2)',car,abar)] ...
        -mp.bet*pscrap*[trmodel.acc_prob_j(mp, (1:abar-2)',car,abar);1];
        y(1)=y(1)-pnew;

        P{tau, car}=x\y;
        a{tau, car}=(1:abar-1); 
      end % end of tau
    end % end of car
    if nargin ==4; 
      P=P{tau, car};
      a=a{tau, car};
    end
  end % end of hom_prices

  function [] = price_graph(mp, tp, car, fignr)
    % price_graph: solve model at parameters mp and plot homogeneous consumer equilibrium prices
    % syntax:
    %     price_graph(mp, tp, car, fignr):
    %         make price_graph for type (tp, car) in figure with number fignr 
    %
    %     price_graph(mp, tp, car):         
    %         make price_graph for type (tp, car) in new figure
    %
    %     price_graph(mp)
    %         make price_graph for type (1, 1) in new figure

    % solve social planners problem (value function, abar and shadow prices)
    [abar_spp]=spp.solve_spp(mp); 
    
    if nargin <3
      tp=1;
    end
    if nargin <4
      car=1;
    end
    % plot Homogeneous consumer equilibrium prices
    if nargin <5;
      fig=figure;
    else
      fig=figure(fignr);
    end
    set(fig,'defaulttextinterpreter','latex');
    % plot price functions for scrap dates <= optimal scrap date
    for abar=abar_spp{tp, car}:-1:abar_spp{tp, car}-5;
      [p, a] = spp.hom_prices(mp, abar, tp, car);
      hold on;
      if abar==abar_spp{tp, car}
        plot((0:abar)',[mp.pnew{car}; p{tp,car}; mp.pscrap{car}], 'LineWidth',3, 'Color', 'r');
        hold on; 
        plot([0 abar_spp{tp, car}+5],[mp.pscrap{car} mp.pscrap{car}], 'LineWidth',3, 'Color', 'k');
      else
        plot((0:abar)',[mp.pnew{car}; p{tp,car}; mp.pscrap{car}], 'LineWidth',1);
      end
      axis('tight');
      xlabel('Age of car');
      ylabel('Equilibrium price of car');
      xlim([0 max(abar_spp{tp, car})+5])
      ylim([-40 mp.pnew{car}*1.1]);
    end
    title(sprintf('Equilibrium (and non-equilibrium) solutions (with $\\bar{a} = \\le $ optimal scrap age)', abar_spp{tp, car}));
    legend({sprintf('Equilibrium price, abar=%i, type=%i, car=%i', abar_spp{tp, car}, tp, car), sprintf('Scrap price (car=%i)', car)})


    % plot price functions for scrap dates >= optimal scrap date
    if nargin <5;
      fig=figure;
    else
      fig=figure(fignr+1);
    end
    set(fig,'defaulttextinterpreter','latex');
    for abar=abar_spp{tp, car}:abar_spp{tp, car}+5;
      [p, a] = spp.hom_prices(mp, abar, tp, car);
      hold on;
      if abar==abar_spp{tp, car}
        plot((0:abar)',[mp.pnew{car}; p{tp,car}; mp.pscrap{car}], 'LineWidth',3, 'Color', 'r');
        hold on; 
        plot([0 abar_spp{tp, car}+5],[mp.pscrap{car} mp.pscrap{car}], 'LineWidth',3, 'Color', 'k');
      else
        plot((0:abar)',[mp.pnew{car}; p{tp,car}; mp.pscrap{car}], 'LineWidth',1);
      end
      axis('tight');
      xlabel('Age of car');
      ylabel('Equilibrium price of car');
      title(sprintf('Equilibrium (and non-equilibrium) solutions (with $\\bar{a} >$ optimal scrap age)', abar_spp{tp, car}));
      ylim([-40 mp.pnew{car}*1.1]);
      xlim([0 max(abar_spp{tp, car})+5])


    end
    legend({sprintf('Equilibrium price, abar=%i, type=%i, car=%i', abar_spp{tp, car}, tp, car), sprintf('Scrap price (car=%i)', car)})
  end % end of price graph

end %methods
end %class
