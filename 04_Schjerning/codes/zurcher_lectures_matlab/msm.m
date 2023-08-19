classdef msm
  % MSM class: Estimates Rust's engine replacement model Rust(Ecta, 1987) 
  %            using method of simulated moments
  % By Fedor Iskhakov, Bertel Schjerning, and John Rust
  methods (Static)

    function [cri,sim_moments,sim_moments_indiv]=objective(data_moments, W, mp, ap, pnames, theta)
      % Computes the objective function of the MSM estimator 
      % Outputs: cri = criterion value
      %     sim_moments = 1 by mp.n vector of moments averaged over all mp.N buses
      %     sim_moments_indiv = N by mp.n matrix of moment differences individually for each bus
      
      persistent ev0; % save ev between function calls to simplify consecutive calls
                     % running <clear msm.objective> to clear saved ev
      if isempty(ev0)
        ev0=0;
      end

      % update model parameters
      mp=vec2struct(theta, pnames, mp);

      % Transition matrix for mileage
      P = zurcher.statetransition(mp);

      % Pay-off
      u = zurcher.u(mp);

      % solve the model
      bellman0= @(V) zurcher.bellman(V, mp, u, P);
      [ev0,pk0,~,~]=dpsolver.poly(bellman0, ev0, ap, mp.beta);

      if mp.useSolutionForSims
        % use the actual stationary distribution instead of simulated moment
        sim_moments = zurcher.eqb(mp,P,pk0)*mp.momentScale;
        sim_moments(end) = []; %drop the last moment as fraction sum up to one
      else
        %do actual simulations
        seed = 12345;
        sims=zurcher.simdata(mp.nsims,mp.Tsims,mp,P,pk0,seed);
        tabs = tabulate(sims.x); %tabulate x
        sim_moments = tabs(:,3)'/100; %simulated fractions
        sim_moments = [sim_moments zeros(1,mp.n-numel(sim_moments))]; %add zeros for the remaining mileage
        sim_moments = sim_moments *mp.momentScale; %scale
        sim_moments(end) = []; %drop the last moment as fraction sum up to one

        %compute the individual moments
        if nargout>2
          g = grpstats(sims,{'id','x'});
          sim_moments_indiv = full(sparse(g.id,g.x,g.GroupCount/mp.Tsims));
          if size(sim_moments_indiv,2)<mp.n
            sim_moments_indiv=[sim_moments_indiv,zeros(mp.nsims,mp.n-size(sim_moments_indiv,2))];
          end
          sim_moments_indiv = sim_moments_indiv*mp.momentScale; %scale
          sim_moments_indiv(:,end) = []; %drop the last moment as fraction sum up to one

          % moment conditions, treat data_moments as identical, with 1 sim per 1 bus in the data
          sim_moments_indiv = sim_moments_indiv - repmat(data_moments,mp.nsims,1); 
        end
      end

      % moment conditions
      mom = (sim_moments - data_moments);
      % criterion with weighting matrix W
      cri = mom * W * mom';
    end

    function [results, pnames, theta_hat, Avar] = estim(data, mp)
      % Estimation using MSM method

      % TODO
    end

  end %methods
end %classdef