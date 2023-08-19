classdef stats
    % this class computes aggregare statistics from the auto trade equilirbrium
    properties (Constant)
        % External cost of Carbon
        social_cost_carbon_dkk = 50 * 5.8 / 1000; % 5.8 was approximately the exchange rate through 2008
        social_cost_carbon_dkk_high = 250 * 5.8 / 1000; % recommended by the environmental councils
        
        % driving externalities , per km
        % change them here to impact the total_externalities and the
        % social_surplus_total measures calculated by compute_outcomes().
        gas_share = 0.80; % share of cars that are gasoline (vs. diesel)
        extern_cost_per_km =          ...
            .0478  ... noise
            + .2095  ... accidents
            + .3368  ... congestion
            + .0097  ... infrastructure
            + .011083*stats.gas_share ...  % air pollution, gas
            + .044565*(1.0-stats.gas_share); % air pollution, diesel
        kg_co2_per_litre_fuel=2.392;
    end
    
    methods (Static)
        
        function out = compute_outcomes(mp, s, sol)
            % Computes a range of outcomes and puts them in a struct
            % with the following fields:
            % consumer_surplus
            % consumer_surplus_tau: [ntypes×1]
            % revenue_car_tau: [ntypes×ncartypes]
            % revenue_fuel_tau: [ntypes×1]
            % vkt_all: [ntypes×ns]
            % vkt_tau
            % fuel_tau
            % co2_tau
            % revenue_fuel
            % revenue_car
            % total_revenue
            % total_co2
            % total_fuel
            % total_vkt
            % extern_co2
            % extern_vkt
            % total_extern
            % social_surplus_total
            %

            out = struct(); 
            
            out.marketshares = sol.marketshares; 
            out.nocar_share = sum(sol.marketshares(:, end)); 
            
            aa = s.is.age(1:end-1); 
            qq = sol.q(1:end-1); 
            out.Ecarage = aa * qq; 
            
            [welfare, welfare_tau] = stats.welfare(mp, sol);
            out.consumer_surplus = (1-mp.bet) * welfare;
            out.consumer_surplus_tau = (1-mp.bet) * welfare_tau;
            
            % compute tax revenue
            [out.revenue_car_tau, out.revenue_fuel_tau, out.markup_car_tau] = stats.revenue(mp, s, sol);
            
            % compute driving, fuel use and co2
            [vkt, fuel, co2] = stats.driving(mp, s);
            out.vkt_all = vkt; % convenient to save
            out.vkt_tau   = zeros(mp.ntypes, 1);
            out.fuel_tau  = zeros(mp.ntypes, 1);
            out.co2_tau   = zeros(mp.ntypes, 1);
            for tau=1:mp.ntypes
                out.vkt_tau(tau)  =  vkt(tau, :) * sol.h_tau{tau};
                out.fuel_tau(tau) =  fuel(tau,:) * sol.h_tau{tau};
                out.co2_tau(tau)  =  co2(tau, :) * sol.h_tau{tau};
            end
            
            % Aggregate outcomes: simple sum because the outcomes are
            % h-weighted.
            out.revenue_fuel  = sum(out.revenue_fuel_tau);
            out.revenue_car   = sum(out.revenue_car_tau);
            out.markup_car    = sum(out.markup_car_tau);
            out.total_revenue = out.revenue_fuel + out.revenue_car;
            out.total_co2     = sum(out.co2_tau);
            out.total_fuel    = sum(out.fuel_tau);
            out.total_vkt     = sum(out.vkt_tau);
            
            % Externalities
            out.extern_co2    = out.total_co2  * stats.social_cost_carbon_dkk ;
            out.extern_vkt    = out.total_vkt  * stats.extern_cost_per_km; 
            out.total_extern  = out.extern_co2 + out.extern_vkt;
            
            % total net welfare
            W = out.consumer_surplus + out.total_revenue;
            out.social_surplus_total_incl_profits = W - out.total_extern - out.markup_car;
            out.social_surplus_total              = W - out.total_extern;
            out.social_surplus_ex_co2             = W - out.extern_vkt; % simply excludes the monetized value of CO2 emissions
            out.social_surplus_total_highco2      = W - out.extern_vkt - out.total_co2 * stats.social_cost_carbon_dkk_high;  % computes the social cost of CO2 at a higher price
        end
        
        function name = get_outcome_name(var)
            % get_outcome_name(var): Returns a human readable name for the
            % variable "var", which is assumed to be a field in the outcome
            % struct computed by "compute_outcomes()". 
            
            names = struct(); 
            
            names.social_surplus_total = 'Social surplus'; 
            names.social_surplus_ex_co2 = 'Social surplus (excl. CO2)'; 
            
            names.total_co2 = 'CO2 emissions';
            names.extern_co2 = 'CO2 externalities (monetary value)'; 
            names.extern_vkt = 'Driving externalities (excl. CO2)';  
            names.total_extern = 'Driving externalities'; 
            
            names.total_revenue = 'Tax revenue';
            names.revenue_fuel = 'Tax revenue from fuel'; 
            names.revenue_car = 'Tax revenue from registrations'; 
            
            names.consumer_surplus = 'Consumer surplus'; 
            names.total_vkt = 'Driving'; 
            
            names.total_fuel = 'Fuel consumption'; 
            
            if isfield(names, var)
                name = names.(var); 
            else
                error('Unable to find a name for the variable %s, please add in the code above.', var);
            end
            
        end
        
        function [vkt, fuel, co2] = driving(mp, s)
            % vkt, fuel and co2 are ntypes*ns
            
            % initialize
            vkt  = nan(mp.ntypes, s.ns);
            fuel = nan(mp.ntypes, s.ns);
            
            for j=1:mp.ncartypes
                idx = s.ipt.car{j};
                age_vec = s.ipt.age(idx);
                for tau=1:mp.ntypes
                    vkt(tau,idx) = trmodel.driving(mp,age_vec,tau,j);
                    fuel(tau,idx) = vkt(tau,idx) / mp.fe{j};
                end
            end
            
            % no car: no driving
            vkt(:, s.ipt.nocar) = 0;
            fuel(:, s.ipt.nocar) = 0;
            
            co2 = fuel * stats.kg_co2_per_litre_fuel;
        end
        
        function [tax_car, markup_car] = get_tax_and_markup(icar, mp)
            if ~isfield(mp, 'passthrough') || mp.passthrough.rate == 1.0 
                % the default mode for the code: no passthrough add-on
                tax_car=mp.pnew{icar}-mp.pnew_notax{icar}; %tax per new car
                markup_car = 0.0; 
            else
                % passthrough enabled: the car price paid by the consumer
                % is composed of three separate parts: 
                %   taxes, 
                %   manufacturer baseline price 
                %   a cange in the manufacturer's markup (relative to
                %   baseline) 
                p0 = mp.pnew_notax{icar};
                p1 = trmodel.pcar_after_passthrough(p0, mp, icar); 
                p2 = mp.pnew{icar};
                tax_car =  p2-p1;   % difference between what the consumer pays and the firm gets 
                markup_car = p1-p0; % difference between the firm's price under full passthrough and the firm's price at imperfect passthrough 
            end
        end
        
        function [revenue_car_tau, revenue_fuel_tau, markup_car_tau] = revenue(mp, s, sol)
            % This function computes the tax revenue from new car tax and from fuel tax
            %
            % OUTPUTS:
            %   revenue_car_tau: (ntypes*1) vector of car tax revenue from
            %       registration taxes paid by consumer type tau=1,...,ntypes.
            %   revenue_fuel_tau: as above, but for fuel taxes.
            %
            %   Note: both are weighted by tw implicitly.
            
            revenue_car_tau_j  = nan(mp.ntypes, mp.ncartypes);
            markup_car_tau_j   = nan(mp.ntypes, mp.ncartypes); 
                
            % New car sales tax
            for j=1:mp.ncartypes
                [tax_car, markup_car] = stats.get_tax_and_markup(j, mp); 
                for tau=1:mp.ntypes
                    revenue_car_tau_j(tau, j)=tax_car*sol.h_tau{tau}(s.ipt.new{j}); %only new cars in holdings
                    markup_car_tau_j(tau, j)=markup_car*sol.h_tau{tau}(s.ipt.new{j}); %only new cars in holdings
                end
            end
            
            % sum over cars 
            revenue_car_tau = sum(revenue_car_tau_j, 2); 
            markup_car_tau = sum(markup_car_tau_j, 2); 
            
            % Fuel tax
            fuel_tax_per_1000_liter = (mp.p_fuel - mp.p_fuel_notax)*1000;
            [vkt, fuel, co2] = stats.driving(mp, s);
            tax_fuel = fuel * fuel_tax_per_1000_liter; % ntypes*ns
            
            revenue_fuel_tau=nan(mp.ntypes, 1); % ntypes*1: revenue from each HH type
            for tau=1:mp.ntypes
                revenue_fuel_tau(tau)=tax_fuel(tau, :)*sol.h_tau{tau};
            end
            
        end      
        
        function [welfare, welfare_tau] = welfare(mp, sol, print_output)
            % stats.welfare: returns the consumer-specific expected value, where
            % expectations are wrt. q (equilibrium quantities).
            %
            % SYNTAX 1: [welfare, welfare_tau] = stats.welfare(mp, sol, print_output)
            % SYNTAX 2: [welfare, welfare_tau] = stats.welfare(mp, sol)
            %
            
            if nargin < 3
                print_output = false;
            end
            
            welfare_tau    = nan(mp.ntypes, 1);
            for tau=1:mp.ntypes
                % Note: sum(q_tau{t}) = 1.0
                welfare_tau(tau) = sol.ev_tau{tau}' * sol.q_tau{tau} / mp.mum{tau};
            end
            
            % aggregate welfare
            welfare = welfare_tau' * mp.tw;
            
            if print_output
                fprintf('PerCapita welfare: %8.4g\n', welfare);
                fprintf('Consumer-specific welfare\n');
                for tau=1:mp.ntypes
                    fprintf(' %23s: %8.4g\n', mp.lbl_types{tau}, welfare_tau(tau));
                end
            end
        end
        
    end % end of methods
end % end of class

