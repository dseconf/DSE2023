classdef tabdat
  % a collection of functions for working with the tabulated data 
  
methods (Static)

  function [h_tau_data, h_tau_data_timeavg, h_tau_nocar] = create_h_tau(dat)
    % OUTPUT: 
    %   h_tau_data: (ntypes*T cell) Fraction of all households that
    %   year having that is value 
                
    assert(all(ismember({'count', 'year', 'is', 'tau'}, dat.Properties.VariableNames)), 'Table must have the following variables: "count", "year", "tau", "is". ');
    
    taus = unique(dat.tau); 
    ntypes = numel(taus); 
    years = unique(dat.year); 
    T = numel(years); 
    
    G = groupsummary(dat(:, {'count', 'year', 's_car_type', 'is', 'tau'}), {'year', 'is', 'tau', 's_car_type'}, 'sum'); 
    h_tau_data = cell(ntypes, T);
    for t=1:T
      for i_tau = 1:ntypes
        this_G = G(G.tau == i_tau & G.year == years(t), :); 
        this_G.sum_count = this_G.sum_count ./ sum(this_G.sum_count); 
        
        % no-car option goes after everything else 
        I_nocar = this_G.s_car_type == -1; 
        I_car   = this_G.s_car_type ~= -1; 
        reorder = [this_G.sum_count(I_car); this_G.sum_count(I_nocar)]; 
        
        % store 
        h_tau_data{i_tau, t} = reorder; 

      end
    end
    
    if nargout > 1
      h_tau_data_timeavg = cell(ntypes, 1);
      
      for i_tau=1:ntypes
        h_tau_data_timeavg{i_tau} = h_tau_data{i_tau, 1}/T;
        % avg. over years
        for t=2:T
          h_tau_data_timeavg{i_tau} = h_tau_data_timeavg{i_tau} + h_tau_data{i_tau, t}/T;
        end
        %h_tau_data_timeavg{i_tau}(s.ipt.age == 0) = nan; % delete incoming cars that are 0 y/o
        
      end
    end
    
    if nargout > 2
      h_tau_nocar = nan(ntypes, T);
      for i_tau=1:ntypes
        for t=1:T
          h_tau_nocar(i_tau, t) = h_tau_data{i_tau, t}(end); 
        end
      end
    end
  end
  
  function [q_tau_data, q_tau_data_timeavg, q_tau_nocar] = create_q_tau(dat)
    % computes the car distributions, q
    % OUTPUT: 
    %   q_tau_data: (ntypes*T cell) Fraction of all households that
    %   year having that is value 
                
    assert(all(ismember({'count', 'year', 'is', 'tau'}, dat.Properties.VariableNames)), 'Table must have the following variables: "count", "year", "tau", "is". ');
    
    taus = unique(dat.tau); 
    ntypes = numel(taus); 
    years = unique(dat.year); 
    T = numel(years); 
    
    G = groupsummary(dat(:, {'count', 'year', 's_car_type','is', 'tau'}), {'year', 'is', 'tau', 's_car_type'}, 'sum'); 
    q_tau_data = cell(ntypes, T);
    for t=1:T
      for i_tau = 1:ntypes
        this_G = G(G.tau == i_tau & G.year == years(t), :); 
        this_G.sum_count = this_G.sum_count ./ sum(this_G.sum_count); 
        
        % no-car option goes after everything else 
        I_nocar = this_G.s_car_type == -1; 
        I_car   = this_G.s_car_type ~= -1; 
        reorder = [this_G.sum_count(I_car); this_G.sum_count(I_nocar)]; 
        
        % store 
        q_tau_data{i_tau, t} = reorder; 

      end
    end
    
    if nargout > 1
      q_tau_data_timeavg = cell(ntypes, 1);
        
      for i_tau=1:ntypes
        q_tau_data_timeavg{i_tau} = q_tau_data{i_tau, 1}/T;
        % avg. over years
        for t=2:T
          q_tau_data_timeavg{i_tau} = q_tau_data_timeavg{i_tau} + q_tau_data{i_tau, t}/T;
        end
        %h_tau_data_timeavg{i_tau}(s.ipt.age == 0) = nan; % delete incoming cars that are 0 y/o          
      end
    end
    
    if nargout > 2
      q_tau_nocar = nan(ntypes, T);
      for i_tau=1:ntypes
        for t=1:T
          q_tau_nocar(i_tau, t) = q_tau_data{i_tau, t}(end); 
        end
      end
    end
  end
  
  function pr_scrap = get_scrap_data(mp, dat)
    % get_scrap_data(mp, dat)
    %
    % INPUTS:
    %   mp: model parametesr
    %   dat: table with data
    %
    % OUTPUT
    %   pr_scrap: (numcar*numage)-by-ntypes matrix
    %       columns indexed by s.is
    %       rows by tau (HH type)
    s = trmodel.index(mp, mp.abar_j0);
    
    % identify observations where an ownership spell ends
    % each surviving row in d pertains to an incoming car that was not kept
    Ihascar = dat.is ~= s.ipt.nocar;
    Inokeep = dat.id ~= s.id.keep;
    Inomiss = not(isnan(dat.scrap));
    d = dat(Ihascar & Inokeep & Inomiss, :);
    
    % # of scrappages in each cell (Coun = #obs, scrap = share of
    % terminating ownerships resulting in scrappage)
    d.num_scrap= d.count .* d.scrap;
    
    % sum counts wtihin units of (tau, is)
    d = groupsummary(d(:, {'tau', 'is', 'num_scrap', 'count'}), {'tau','is'}, 'sum');
    
    pr_scrap = nan(s.is.nocar-1, mp.ntypes);
    for itau=1:mp.ntypes
      d_ = d(d.tau == itau, :);
      ii = d_.is; % not all rows will be filled out since sometimes, scrap is missing or no HH of this type terminates ownership of that car
      pr_scrap(ii, itau) = d_.sum_num_scrap ./ d_.sum_count;
    end            
  end
end
end
