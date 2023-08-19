classdef bbl
  % BBL class: Estimates Rust's engine replacement model Rust(Ecta, 1987) 
  %            using a variant of Bajari, Benkard, Levin (Ecta, 2007) method
  % By Fedor Iskhakov, Bertel Schjerning, and John Rust
  methods (Static)

    function res=objective(pol0, pol1, H, P, mp, pnames, theta)
      % Computes the objective function of the BBL estimator 
      % (squared one-sided deviation of the moment inequalities)
      % pol0 = policy estimated at first stage (CCP)
      % pol1 = perturbed policy
      % H = distribution over the set of inequalities (H(x) in BBL2007)
      % P = transition probabilities

      % update model parameters
      mp=vec2struct(theta, pnames, mp);
      J = numel(pol1); %number of perturbed policies

      % Use Hotz-Miller inversion to get integrated values for the two policies
      V0 = npl.phi(mp,pol0,P);
      res=0;
      for j=1:J
        V1 = npl.phi(mp,pol1{j},P);
        % BBL moment inequality
        g = V0-V1;
        gg = (min(g,0)).^2;
        res = res + H * gg / J; %uniform distribution over perturbations
      end
    end

    function [results, pnames, theta_hat, Avar] = estim(data, mp)
      % Estimation using BBL method

      % TODO
    end

  end %methods
end %classdef