%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef graphs
    
properties (Constant)
  fontsize = 18; % common across figures (unless explicitly overridden) 
  fontsize_3_wide = 22; % when there are 3 plots next to each other on a page 
  linewidth = 2; 
  fontname = 'times'; 
end %properties

methods (Static, Access=public) %methods callable from the outside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [] = outcomes(mp, s, plots, sol_m, sol_d, pltype, sgtitles)
  % SYNATX : graphs.outcomes(mp, s, plots, sol_m, sol_d, sgtitles)
  %
  % INPUTS: 
  %   pltype: e.g. {'Model', 'Data'}
  %   sgtitles: [] (empty) for no header; otherwise cell array of titles
  if nargin<6 || isempty(pltype)
    pltype={'Model'};
  end
  if nargin<7
    % just print 'Model' and 'Data' as the sgtitle headers 
    sgtitles = pltype; 
  else 
    assert(iscell(sgtitles) || isempty(sgtitles), 'Input "sgtitles" must be a cell array (of length 1 or 2) or empty (for no header)'); 
  end


  if isempty(plots)

      plots={'prices', % equilibrium prices
             'post_trade_dist', % aggregate distribution (post trade)
             'agg_holdings', % aggregate distribution (post trade)
             'holdings', % aggregate distribution by car types (post trade)
             'keep', % keep probabilities
             'scrap', % scrap probabilities
             'taxes' % distribution of tax revenues
             'choice_clunker' % choosing clunker
             'choice_nocar' % choosing no car
      };
  end


  for i=1:numel(plots)

    % find maximum age in the data (for x axes - in case model was solved with more car ages than are in the data)
    if ismember('Data', pltype)
      ncartypes_d = size(sol_d.marketshares, 2)-1; % 1 for outside 
      abar_d = (numel(sol_d.q)-1) / ncartypes_d; 
      if abar_d ~= round(abar_d)
        mp.max_car_age = ceil(abar_d); 
        warning('Confusing maximum car age implied by sol_d, proceding with max_car_age = %8.4f', mp.max_car_age)
      end
    else
      mp.max_car_age = max([mp.abar_j0{:}]); 
    end

    for j=1:numel(pltype)
      switch pltype{j}
        case 'Model'
          sol=sol_m;
        case 'Data'
          sol=sol_d;
        case 'Data-Model'
          sol=struct_diff(sol_d, sol_m);
        otherwise
          disp('pltype must be "Model", "Data", or "Data-Model"')
          return
      end

      switch plots{i}
        case 'prices'
            graphs.prices(mp, s, sol);
        case 'agg_holdings' 
            graphs.agg_holdings(mp, s, sol);
        case 'post_trade_dist' 
            graphs.post_trade_dist(mp, s, sol);
        case 'holdings' 
            graphs.holdings(mp, s, sol);
        case 'quantity' 
          graphs.quantity(mp, s, sol);
       case 'keep' 
            graphs.ccp(mp, s, sol, 'keep')
       case 'scrap' 
            graphs.ccp_scrap(mp, s, sol, pltype{j})
       case 'choice_clunker' 
            graphs.ccp(mp, s, sol, 'choice_clunker')
       case 'choice_nocar' 
            graphs.ccp(mp, s, sol, 'choice_nocar')
       case 'taxes'
            graphs.taxes(mp, s, sol);
       case 'ev_gain'
          if strcmp(pltype{j}, 'Model') 
            graphs.ev_gain(mp);
          end
          otherwise
          error('plot does not exist')
      end
      if ~isempty(sgtitles) && ~strcmp(plots{i}, 'prices') % no header for prices as we have no corresponding data to show 
        sgtitle(sgtitles{j}); 
      end
    end
  end
end %function

function f=myfigure(varargsin)
    % INPUTS: 
    %   varargsin: [optional] cell array of figure properties 
    % OUTPUT: 
    %   f: figure handle 
  if nargin>0
    f=builtin('figure',varargsin);
  else
    f=builtin('figure');
  end
  
  set(gcf,'Color',[1 1 1]);

  % make figures take different positions on screen
  nf=numel(findobj('type','figure'));
  pos=get(gcf,'Position'); %[left bottom width height]
  step=100;
  set(gcf,'Position',[(nf-1)*step 1000 pos(3:4)]);
end

function f = prices(mp, s, sol)
	% SYNTAX: graphs.prices(mp, s, so)
  %
  % OUTPUT: 
  %   f: figure handle 
  if ~isfield(sol, 'p')
    warning('graphs.prices: prices not available: sol struct does not have field "p": exiting') 
    return
  end     

  f=graphs.myfigure;
  
  hold on;
  for j=1:mp.ncartypes
    pj=[mp.pnew{j}; sol.p(s.ip{j}); mp.pscrap{j}];
    h=plot([0; (1:s.abar_j{j})'],pj,'Linewidth', graphs.linewidth);
    if mp.ncartypes>1
      set(h,'DisplayName',mp.lbl_cartypes{j});
    else
      set(h,'DisplayName','Heterogeneous');
    end
  end
  
  legend('Location','Northeast');
  ylabel('Price, 1000 DKK'); 
  xlabel('Car age');
  axis('tight')
  hold off;
  set(gca, 'Ygrid', 'on'); 
  
  graphs.set_fig_layout_post(f); 
  
end

function f = post_trade_dist(mp, s, sol) 
    % syntax graphs.post_trade_dist(mp, s, sol)
  % DISTRIBUTION PLOT: ALL ON ONE AXES
  %
  % OUTPUT: 
  %   f: figure handle 

  h_tau=sol.h_tau;

  % allow for multiple types of car on one graph 
  % (because outside option is in the same distribution!)
  graphs.myfigure(); 
  clf
  set(gcf,'Color',[1 1 1]);
  colormap(summer);
  hold on;

  dty=[]; 
    for j=1:mp.ncartypes
    tmp=[];
        for tau=1:mp.ntypes
        tmp(:,tau)=[100*h_tau{tau}(s.ipt.car{j}); nan];
        end
        dty=[dty; tmp];
  end

  % no car state - from h_tau
  qnc=0;
  for t=1:mp.ntypes
    qnocar(t)=100*h_tau{t}(s.ipt.nocar); %outside option
  end
  qnc=sum(qnocar);

  %split the no car fraction into many columns of height maxy
  maxy=max(sum(dty(1:end-1,:),2)); %max fraction over types of consumer and car
  maxy=ceil(maxy*1000)/1000;  %round up at 3rd diget 
  maxy=1;  % 1 pct.
  nck=floor(qnc/maxy); %number of columns of hight maxy
  nckfrac=maxy*nck/qnc; %fraction of no car prob mass in repeated columns
  dty_nocar=repmat(qnocar*nckfrac/nck,nck,1);
  dty_nocar(end+1,:)=qnocar*(1-nckfrac);

  % add to y values
  dty=[dty;dty_nocar];

  %x ticks and labels
  dtx=[1:numel(h_tau{1})+nck+mp.ncartypes]; %all distributions + nocar columns + separators
  % location of tick marks
  tickstep=4;
  xticks=dtx(1:tickstep:numel(h_tau{1}));  % cars 
  xticks(end+1) = numel(h_tau{1})+mp.ncartypes+floor(nck/2); % no car
  xticks=[10 35 65 90 120]; 

  % tick labels
  ticklabels={};
  for j=1:mp.ncartypes
    lbl_used=compose('%d',[tickstep:tickstep:s.abar_j{j}]);
        %ticklabels= {ticklabels{:}, '0', lbl_used{:}};
        ticklabels= {ticklabels{:}, mp.lbl_cartypes{j}};
  end
  ticklabels{end+1}='No car';

  xlim([dtx(1)-.5,dtx(end)+.5]);
  set(gca,'XTick',xticks,'XTickLabel',ticklabels,'Ygrid','on','TickLength',[0,0]);
  h1=bar(dtx,dty,'LineWidth',1,'BarWidth',1,'BarLayout','stacked');
  l1=legend(flip(h1), flip(mp.lbl_types));

  set(l1,'EdgeColor',[1 1 1],'Location','northeastoutside');
  xlabel('Car age');
  % title(sprintf('{\\bf Holdings distributions (start of period) \\sigma=%g, T=%g \\rho=%g}',mp.sigma,mp.transcost,mp.ptranscost));
  ylabel('Fraction of population (%)');
  axis('tight');
  graphs.set_fig_layout_post;

  hold off;
  
end % post_trade_dist() 


function f = prices_with_hom(mp, s, sol)
    assert(mp.ntypes == 2 && mp.ncartypes == 1, 'Only implemented for 2*1 economy');
    
    % prices in solution 
    pp = [mp.pnew{1}; sol.p]; 
    
    % compute equilibria for each of the 1-type economies 
    pp_onetype = graphs.compute_hom_equilibria(mp);
    
    f = graphs.myfigure();
    aa = (0:mp.abar_j0{1}-1)';
    
    plot(aa, [pp, pp_onetype], 'LineWidth', graphs.linewidth);
    leg_ = {'Two-type equilibrium', ['One-type: ', mp.lbl_types{1}], ['One-type: ', mp.lbl_types{2}]};
    legend(leg_);
    
    ylabel('Price, 1000 DKK'); 
    xlabel('Car age'); 
    set(gca,'Ygrid','on','box','off'); 
    
    graphs.set_fig_layout_post(f); 
 
end

function f = agg_holdings(mp, s, sol) 
	% syntax graphs.holdings(mp, s, sol)
  % DISTRIBUTION PLOT: ALL ON ONE AXES
  %
  % OUTPUT: 
  %   f: figure handle 

  h_tau=sol.h_tau;

  % allow for multiple types of car on one graph 
  % (because outside option is in the same distribution!)
  f=graphs.myfigure(); 
  hold on;

  dty=[]; 
	for j=1:mp.ncartypes
    tmp=[];
		for tau=1:mp.ntypes
	  	tmp(:,tau)=[h_tau{tau}(s.ipt.car{j}); nan];
		end
		dty=[dty; tmp];
  end

  % no car state - from h_tau
  qnc=0;
  for t=1:mp.ntypes
    qnocar(t)=h_tau{t}(s.ipt.nocar); %outside option
  end
  qnc=sum(qnocar);

  %split the no car fraction into many columns of height maxy
 	maxy=max(sum(dty(1:end-1,:),2)); %max fraction over types of consumer and car
  maxy=ceil(maxy*1000)/1000;  %round up at 3rd diget 
  nck=floor(qnc/maxy); %number of columns of hight maxy
  nckfrac=maxy*nck/qnc; %fraction of no car prob mass in repeated columns
  dty_nocar=repmat(qnocar*nckfrac/nck,nck,1);
  dty_nocar(end+1,:)=qnocar*(1-nckfrac);

  % add to y values
  dty=[dty;dty_nocar];

  %x ticks and labels
  dtx=[1:numel(h_tau{1})+nck+mp.ncartypes]; %all distributions + nocar columns + separators
  % location of tick marks
  tickstep=4;
  xticks=dtx(1:tickstep:numel(h_tau{1}));  % cars 
  xticks(end+1) = numel(h_tau{1})+mp.ncartypes+floor(nck/2); % no car
  %xticks=[10 35 65 90 120]; 

  % tick labels
  
  ticklabels={};
  if mp.ncartypes > 1
  for j=1:mp.ncartypes
  	lbl_used=compose('%d',[tickstep:tickstep:s.abar_j{j}]);
        ticklabels= {ticklabels{:}, mp.lbl_cartypes{j}};
  end
  else % just one car 
      ticklabels = compose('%d',[0:tickstep:s.abar_j{1}]);
      
      % hotfix: delete the last tick to avoid clash with 'No car'
      ticklabels{end} = ''; 
  end
  ticklabels{end+1}='No car';
  
  
  xlim([dtx(1)-.5,dtx(end)+.5]);
  set(gca,'Ygrid','on','TickLength',[0,0],'XTick',xticks,'XTickLabel',ticklabels); 
  bar(dtx,dty,'LineWidth',1,'BarWidth',1,'BarLayout','stacked');
  l1=legend(mp.lbl_types{:});
  set(l1,'Location','northeastoutside');
  xlabel('Car age');
  % title(sprintf('{\\bf Holdings distributions (start of period) \\sigma=%g, T=%g \\rho=%g}',mp.sigma,mp.transcost,mp.ptranscost));
  ylabel('Fraction of the population');
  axis('tight');
  graphs.set_fig_layout_post;
  hold off;
  
  graphs.set_fig_layout_post(f); 
  
end % agg_holdings() 

function set_fig_layout_post(f, IS3D, fontsize)
    % INPUTS: 
    %   f: figure handle
    %   IS3D: (boolean, default=false) Switch, set =true if the graph is 3d 
    %   fontsize: (scalar) size of the font in the entire figure 
    %
    % example: 
    % graphs.set_fig_layout_post(gcf, false, 25)
    
    if nargin < 2 
        IS3D = false; 
    end
    if nargin < 3 
        fontsize = graphs.fontsize; 
    end
    
    % global settings 
    set(gca, 'fontname', graphs.fontname, 'box', 'off'); 
    
    if IS3D
        % no special 3d graph settings (yet)
    else
        set(gca, 'Ygrid', 'on'); 
    end
    
    set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
    
end

function h = waitbar(msg, it, nit, h, wait_times)

  if nargin == 1 || it == 1 
    % first iteration
    h = waitbar(0, msg); 
    return 
  else
    % subsequent
    avg_time_per_it = mean(wait_times(1:(it-1))); 
    remaining_it = nit - it; 
    remaining_time = avg_time_per_it * remaining_it; 

    msg_full = sprintf('%s: %d/%d (ETA: %5.2f min)', msg, it, nit, remaining_time/60.); 

    waitbar(it/nit, h, msg_full); 
  end
  
end


function [] = quantity(mp, s, sol, REMOVENOCAR, NORMALIZEBYCATEGORY)
    % dist_by_car_type: Car quantity distribution (State) - BY CAR TYPE
    % SYNATX : graphs.quantity(mp, s, sol)

    q_tau=sol.q_tau;
    if nargin < 4
        REMOVENOCAR = true; 
    end
 
    if nargin < 5
        NORMALIZEBYCATEGORY = true;
    end
    graphs.myfigure();
    for i=1:mp.ncartypes
        subplot(2,2,i); 
        set(gcf,'Color',[1 1 1]);
        colormap(summer);
        hold on;
        dtx=s.is.age(s.is.car{i});
        dtx=[dtx'; max(dtx)+1];
        xlim([dtx(1)-.5,dtx(end)+.5])
        tickstep=2;
        xticks=dtx(1:tickstep:numel(dtx));  % cars
        ticklabels=compose('%d',xticks);
        
        if REMOVENOCAR
            ticklabels{end}=[];
        else
            ticklabels{end}='No car';
        end
        
        set(gca,'XTick',xticks,'XTickLabel',ticklabels,'Ygrid','on');
        dty=[];
        for t=1:mp.ntypes
            dty(:,t)=[q_tau{t}(s.is.car{i}); q_tau{t}(s.is.nocar)]*100;
        end

        if NORMALIZEBYCATEGORY
            dty = dty ./ sum(dty, 2); 
        end
        
        if REMOVENOCAR 
            dty(end, :) = [];
            dtx(end, :) = []; 
        end
        
        if mp.ntypes==numel(mp.lbl_types)
            leglabels=mp.lbl_types;
        else
            leglabels = sprintfc('Household %d', 1:mp.ntypes); 
        end
        
        bar(dtx,dty,'LineWidth',1,...
            'BarWidth',1,...
            'BarLayout','stacked');
        if i == 1 
            l1=legend(leglabels);
            set(l1,'EdgeColor',[1 1 1],'Location','northeastoutside');
        end
        title(sprintf('Holdings: %s', mp.lbl_cartypes{i}));
        
        xlabel('Age of car');
        if NORMALIZEBYCATEGORY
            ylabel('Share of households in group');
        else
            ylabel('Fraction of population');
        end
        axis('tight');
        hold off;
    end
end

function [dta, mp] = holdings(mp, s, sol, REMOVENOCAR, NORMALIZEBYCATEGORY)
    % dist_by_car_type: Car distribution post trade - BY CAR TYPE
    % SYNATX : graphs.holdings(mp, s, h_tau)
    if nargin < 4
        REMOVENOCAR = true; 
    end
 
    if nargin < 5
        NORMALIZEBYCATEGORY = true;
    end

    h_tau=sol.h_tau;
    graphs.myfigure();
    for i=1:mp.ncartypes
        subplot(2,2,i); 
        hold on;
        dtx=s.ipt.age(s.ipt.car{i});
        dtx=[dtx'; max(dtx)+1];
        xlim([dtx(1)-.5,dtx(end)+.5])
        tickstep=2;
        xticks=dtx(1:tickstep:numel(dtx));  % cars
        ticklabels=compose('%d',xticks);
        
        if REMOVENOCAR
            ticklabels{end}=[];
        else
            ticklabels{end}='No car';
        end
        
        set(gca,'XTick',xticks,'XTickLabel',ticklabels,'Ygrid','on');
        dty=[];
        for t=1:mp.ntypes
            dty(:,t)=[h_tau{t}(s.ipt.car{i}); h_tau{t}(s.ipt.nocar)]*100;
        end

        if NORMALIZEBYCATEGORY
            dty = dty ./ sum(dty, 2); 
        end
        
        if REMOVENOCAR 
            dty(end, :) = [];
            dtx(end, :) = []; 
        end
        
        if mp.ntypes==numel(mp.lbl_types)
            leglabels=mp.lbl_types;
        else
            leglabels = sprintfc('Household %d', 1:mp.ntypes); 
        end
        
        bar(dtx,dty,'LineWidth',1,...
            'BarWidth',1,...
            'BarLayout','stacked');
        if i == 1 
            l1=legend(leglabels,'Location','southeast');
        end
        set(l1,'EdgeColor',[1 1 1]);
        % title(sprintf('Holdings: %s', mp.lbl_cartypes{i}));
        title(mp.lbl_cartypes{i});
        xlabel('Age of car');
        if NORMALIZEBYCATEGORY
            ylabel('Fraction of households in group');
        else
            ylabel('Fraction of population');
        end
        axis('tight');
        graphs.set_fig_layout_post;
        hold off;
    end
end

function [revenue_car, revenue_fuel, total_rev] = taxes(mp, s, sol)
  h_tau=sol.h_tau;
  [revenue_car, revenue_fuel] = stats.revenue(mp, s, sol);
  total_rev=sum(revenue_car(:))+sum(revenue_fuel(:));

  %to be stacked by consumer type
  ydata = [(revenue_car/total_rev)'; nan(1, mp.ntypes); (revenue_fuel/total_rev)']; 
  % ydata = [(revenue_car)'; nan(1, mp.ntypes); (revenue_fuel)']; 
  xdata = [1:2*mp.ncartypes+1];

  xticklabels = cell(1, 2*mp.ncartypes + 1); 
  xticklabels(1:mp.ncartypes) = mp.lbl_cartypes; 
  xticklabels(end-mp.ncartypes+1:end) = xticklabels(1:mp.ncartypes); 

  graphs.myfigure;
  clf
  hold on
  set(gca,'XTick',xdata,'XTickLabel',xticklabels,'Ygrid','on',...
          'TickLength',[0,0],'XLim',[xdata(1)-1,xdata(end)+1]);
  bar(xdata,ydata,'LineWidth',1,...
              'BarWidth',0.95,...
              'BarLayout','stacked');
  yl=get(gca,'YLim');
  set(gca,'YLim',[0 ceil(yl(2)*10)/10]);
  l1=legend(mp.lbl_types{:});
  set(l1,'EdgeColor',[1 1 1],'Location','Best');
  title(sprintf('Distribution of tax revenue over sources, consumer and car types'));
  ylabel('Fraction of total tax revenue');
  xlabel(sprintf('New car tax %s Fuel tax',repmat(' ',1,40)));
  hold off;

end

function [] = fit_scrap_dist(mp, s, sol_d, sol_m)

  graphs.myfigure(); 
  sgtitle('Distribution of scrappage') 

  ccp_keep_m=nan(s.ns, mp.ntypes); 
  ccp_keep_d=nan(s.ns, mp.ntypes); 
  for tau=1:mp.ntypes
    ccp_keep_m(:,tau)=sol_m.ccp_tau{tau}(:,s.id.keep);
    ccp_keep_d(:,tau)=sol_d.ccp_tau{tau}(:,s.id.keep);
  end

  scrap_d=([sol_d.q_tau{:}].*[sol_d.ccp_scrap_tau{:}].*(1-ccp_keep_d))*mp.tw;
  scrap_m=([sol_m.q_tau{:}].*([sol_m.ccp_scrap_tau{:}].*(1-ccp_keep_m).*(1-sol_m.F.accprob) + sol_m.F.accprob))*mp.tw;

  for j=1:mp.ncartypes
      subplot(2,2,j)
      idx=s.is.car_ex_clunker{j};
      if isfield(mp, 'max_car_age') 
        idx = s.is.age <= mp.max_car_age; 
      end

      bar([scrap_d(idx) scrap_m(idx)])
      title(mp.lbl_cartypes{j});
      if j==mp.ncartypes
          legend('Data', 'Model');
      end
  end
end


function [] = ccp_scrap(mp, s, sol, pltype, maxcarage)
  ccp_tau=sol.ccp_tau; 
  ccp_scrap_tau=sol.ccp_scrap_tau;
  if nargin<4
    pltype='Data';
  end

  str_title0= 'Scrap probability';
  
  dims=[ceil(sqrt(mp.ncartypes)), round(sqrt(mp.ncartypes))]; 
  graphs.myfigure(); 
  tickstep=2;

  for icar=1:mp.ncartypes
    subplot(dims(1),dims(2), icar);
    str_title= sprintf('%s (%s)',str_title0, mp.lbl_cartypes{icar});
    hold on;
    abarmax=s.abar_j{icar}; 
    for t=1:mp.ntypes
      idx=s.is.car_ex_clunker{icar};
      hold on;
      % probability of scrapping the existing car by consumers
      dty=nan(1, abarmax);
      dtx=nan(1, abarmax);
      if strcmp(pltype, 'Data')
        % dty(1:numel(idx))=ccp_scrap_tau{t}(idx); 
        dty(1:numel(idx))=(1-ccp_tau{t}(idx, s.id.keep)).*ccp_scrap_tau{t}(idx); 
      elseif strcmp(pltype, 'Model')

        % (1-phatkeep(c,a))*[acc(c,a-1)+(1-acc(c,a-1))*P(voluntary scrap of car (c,a) at time t))] 
        dty(1:numel(idx))= (1-ccp_tau{t}(idx, s.id.keep)).*(sol.F.accprob(idx)+(1-sol.F.accprob(idx)).*ccp_scrap_tau{t}(idx));
        %dty(1:numel(idx))=ccp_scrap_tau{t}(idx).*(1-sol.F.accprob(idx)) + sol.F.accprob(idx); 
      else
        dips('graphs.ccp_scrap: pltype must be "Model" or "Data"')
        return        
      end
      dtx(1:numel(idx))=s.is.age(idx); % car ages 1 to abar-1 

      h=plot(dtx,dty,'LineWidth',2,'Marker','o','MarkerSize',5,'MarkerFaceColor','auto');

      if (~strcmp(mp.lbl_cartypes{icar},''))
      set(h,'DisplayName',sprintf('%s',mp.lbl_types{t}));
      else
      set(h,'DisplayName',sprintf('%s, %s',mp.lbl_types{t}, mp.lbl_cartypes{icar}));
      end
    end
    xlim([0,abarmax-1]);

    title(str_title);
    if icar==1
      l1=legend('Location','Northwest');
    end
    ylabel(sprintf('Probability of scrap (%s)', mp.lbl_cartypes{icar}));
    xlabel(sprintf('Age of car (%s)', mp.lbl_cartypes{icar}));
    hold off;

  end
end

function [] = accident_rates(mp,s,sol_m)

  graphs.myfigure(); 
  hold on;
  linetype=cell(mp.ncartypes,1);
  for ct=1:mp.ncartypes
    if (ct == 1)
      linetype{ct}='k--';
    elseif (ct == 2)
      linetype{ct}='b--';
    elseif (ct == 3)
      linetype{ct}='k-';
    else
      linetype{ct}='b-';
    end
 
    h=plot(s.is.age(s.is.car{ct}(1:end-1))',sol_m.F.accprob(s.is.car{ct}(1:end-1)),linetype{ct},'Linewidth',2);
    set(h,'DisplayName',sprintf('%s',mp.lbl_cartypes{ct}));
 
  end
  l1=legend('Location','Northwest');
  title('Estimated accident probabilities by age and car type');
  xlabel('Age of car');
  ylabel('Annual probabilty car is totaled in an accident');
  set(gca,'Ygrid','on');
  hold off;


end

function [] = ccp_scrap_compare(mp, s, sol_d, sol_m)
  ccp_tau_d=sol_d.ccp_tau; 
  ccp_scrap_tau_d=sol_d.ccp_scrap_tau;
  ccp_tau_m=sol_m.ccp_tau; 
  ccp_scrap_tau_m=sol_m.ccp_scrap_tau;

  str_title0= 'Scrap probability';
  
  dims=[ceil(sqrt(mp.ncartypes)), round(sqrt(mp.ncartypes))]; 
  graphs.myfigure(); 
  %tickstep=2;

  for icar=1:mp.ncartypes
    str_title= sprintf('%s (%s)',str_title0, mp.lbl_cartypes{icar});
    str_title= mp.lbl_cartypes{icar};
    abarmax=s.abar_j{icar}; 
    abarmax=21;
    idx=s.is.car_ex_clunker{icar}(1:21);
    dtx(1:numel(idx))=s.is.age(idx)'; % car ages 1 to abar-1 
    dty_d=zeros(abarmax, 1);
    dty_m=zeros(abarmax, 1);

    for t=1:mp.ntypes
        dty_d=dty_d+(1-ccp_tau_d{t}(idx, s.id.keep)).*ccp_scrap_tau_d{t}(idx)*mp.tw(t); 
        dty_m=dty_m+(1-ccp_tau_m{t}(idx, s.id.keep)).*(sol_m.F.accprob(idx)+(1-sol_m.F.accprob(idx)).*ccp_scrap_tau_m{t}(idx))*mp.tw(t);
    end

    subplot(dims(1),dims(2), icar);

    hL= plot(dtx,dty_d, 'r-', dtx,dty_m, 'k-');
    set(hL, {'markerfacecolor'},{'r';'k'}, {'Marker'},{'o'; '+'}, ...
        'LineWidth',2, 'MarkerSize',5);

    if icar==1;
     l1=legend('Data', 'Model');
     set(l1,'EdgeColor',[1 1 1],'Location','Northwest');
    end

    title(str_title);
    ylabel('Scrap Probability');
    xlabel('Age of car');
    axis([0 abarmax 0 .4]);
    %set(gca,'YGrid','on');
    %legend('Model','Data','Location','Northwest');
  end
end

function [] = ccp(mp, s, sol, choice)

  ccp_tau=sol.ccp_tau;

  dims=[ceil(sqrt(mp.ncartypes)), round(sqrt(mp.ncartypes))]; 
  graphs.myfigure(); 
  tickstep=2;

  for icar=1:mp.ncartypes;
    subplot(dims(1),dims(2), icar);
    set(gcf,'Color',[1 1 1]);
    colormap(summer);
    abarmax=s.abar_j{icar}; 
    str_x=''; str_y=''; 
    str_title=mp.lbl_cartypes{icar}; 
    for t=1:mp.ntypes;
      hold on;
      switch choice
        case 'keep'   % probability of keeping the existing car by consumers    
          str_title=mp.lbl_cartypes{icar}; 
          str_y='Probability of keeping'; 
          str_x='Age of current car';
          dty=ccp_tau{t}([s.is.car_ex_clunker{icar}],s.id.keep); %first column, excluding the abar car which can not be kept
          dtx=s.is.age(s.is.car_ex_clunker{icar}); % car ages 1 to abar-1

        case 'choice_clunker'   % probability of purchasing car when owning a clunker     
          str_y='Probability of purchaing a car (same type)'; 
          str_x=sprintf('Age of purchased car'); 
          str_title=sprintf('Consumers owning clunkers (%s)', mp.lbl_cartypes{icar}); 
          dty=ccp_tau{t}([s.is.clunker{icar}],s.id.trade{icar}); %first column, excluding the abar car which can not be kept
          dtx=s.id.age(s.id.trade{icar}); % car ages 0 to abar-1 

        case 'choice_nocar' % probability of purchasing car when owning no car
          str_y='Probability of purchasing a car'; 
          str_x=sprintf('Age of purchased car'); 
          str_title=sprintf('%s', mp.lbl_cartypes{icar}); 
          dty=ccp_tau{t}([s.is.nocar],s.id.trade{icar}); %first column, excluding the abar car which can not be kept
          dtx=s.id.age(s.id.trade{icar}); % car ages 0 to abar-1 
        otherwise
          warning('Choice must be "keep", "choice_clunker", or "choice_clunker"')
          return
      end

      h=plot(dtx,dty,'LineWidth',2,'Marker','o','MarkerSize',5,'MarkerFaceColor','auto');

      if (~strcmp(mp.lbl_cartypes{icar},''))
        set(h,'DisplayName',sprintf('%s',mp.lbl_types{t}));
      else
        set(h,'DisplayName',sprintf('%s, %s',mp.lbl_types{t}, mp.lbl_cartypes{icar}));
      end
      graphs.set_fig_layout_post;

    end
    xlim([dtx(1)-.5,dtx(end)+.5]);
    xticks=dtx(1:tickstep:numel(dtx));  
    ticklabels=compose('%d',xticks);
    set(gca,'XTick',xticks,'XTickLabel',ticklabels,'YGrid','on');
    title(str_title); xlabel(str_x); ylabel(str_y);
    if icar==mp.ncartypes;
      l1=legend('Location','Northeast');
    end
    hold off;
  end
end

function [] = ev_gain(mp)
  warning('graphs.ev_gain no longer functional')
  return
  mp.es=0;
  mp_es=mp; 
  mp_es.es=1;
  
  [abar_spp, p_spp, w_spp]=spp.solve_spp(mp);
  [s_0, F_0, price_j_0, q_tau_0, q_0, ev_tau_0]=equilibrium.solve(mp);
  [s_es, F_es, price_j_es, q_tau_es, q_es, ev_tau_es]=equilibrium.solve(mp_es);
  [s_max, F_max, price_j_max, q_tau_max, q_max, ev_tau_max]=maximal_equilibrium.solve(mp);
  
  
  % PLOT of gains in expected value
  for j=1:mp.ncartypes;
      for t=1:mp.ntypes;
          graphs.myfigure;
          clf
          hold on;
          plot((1:s_0.abar_j{j})',ev_tau_0{t}(s_0.is.car{j}),'Linewidth',2);
          plot((1:s_es.abar_j{j})',ev_tau_es{t}(s_es.is.car{j}),'Linewidth',2);
          plot((1:s_max.abar_j{j})',ev_tau_max{t}(s_max.is.car{j}),'Linewidth',2);
          legend('No scrapping','Endogenous Scrap','Maximal','Location','Northeast');
          xlabel('Age of car');
          ylabel('Expected utility');
          title({sprintf('Expected utility for owners of car type %i, consumer type %i',j,t), ...
              sprintf('EV outside good: ') , ...
              sprintf('no scrapping: %g', ev_tau_0{t}(s_0.is.nocar)), ...
              sprintf('endgenous scrap: %g', ev_tau_es{t}(s_es.is.nocar)), ...
              sprintf('maximal equilibrium: %g', ev_tau_max{t}(s_max.is.nocar)), ...
              });
          axis('tight');
          hold off;
      end
  end
end

function fit_scrap0(mp, dat, ccp_scrap_tau)
    
    if nargin < 3 
        [s, F, price_j, q_tau, q, ev_tau, ccp_tau, ccp_scrap_tau, ctp_tau, ed, h_tau, marketshares]=equilibrium.solve(mp);
    else
        s = trmodel.index(mp, mp.abar_j0); 
    end
    
    % data: in ncars*nHH
    yy_d = tabdat.get_scrap_data(mp, dat); 

    % model 
    yy_m = nan(s.ns-1, mp.ntypes); 
    for tau=1:mp.ntypes
        for j=1:mp.ncartypes
            jj = s.is.car_ex_clunker{j}; 
            yy_m(jj,tau) = ccp_scrap_tau{tau}(jj, 1);
        end
    end

    figure(); 
    for j=1:mp.ncartypes
        %xx = (1:mp.abar_j0{j})';
        ii = s.is.car_ex_clunker{j}; % indices 
        xx = s.is.age(ii); 
        subplot(2,mp.ncartypes,j);
        plot(xx, yy_d(ii, :)); title(['Data: ' mp.lbl_cartypes{j}]); ylabel('Non-kept cars that are scrapped');
        subplot(2,mp.ncartypes,j + mp.ncartypes);
        plot(xx, yy_m(ii, :)); title(['Model: ' mp.lbl_cartypes{j}]); ylabel('Pr(scrap)');
        xlabel('Age of incoming car');
        
        if j == 1
            legend(mp.lbl_types);
        end
    end
end

function data_nocar_over_time(mp, s, dta)
    % SYNTAX: graphs.data_nocar_over_time(mp, s, dat)
    % Share of no-car ownership
    years = unique(dta.year);
    T = numel(years);
    share_nocar_dat = nan(mp.ntypes, T);
    
    % no-car

    for t=1:numel(years)
        for tau = 1:mp.ntypes
          tab=data.tabulate(dta(dta.tau == tau & dta.year == years(t),:), {'is'});
          share_nocar_dta(tau, t) = tab.share(tab.is==s.is.nocar);
        end
    end
    
    % show no-car holdings over time
    if numel(years) > 1
        set(0,'defaultAxesFontSize',14);
        figure()
        [~, ii_] = sort(share_nocar_dta(:,end));
        plot(years, share_nocar_dta(ii_, :), '-o');
        legend(mp.lbl_types{ii_(end:-1:1)}); axis('tight')
        ylabel('Share of households with no car');
    end
end

function fit_driving(mp, dta, sol_m, outcomes_dta, outcomes_mod)
    % show predicted driving
    % TODO: implement driving in sol_m and sol_d structs. 
    
    assert(nargin >= 3, 'Input sol_m required');
    
    if nargin < 4
        DOPRINT=true;
        s = trmodel.index(mp, mp.abar_j0); 
        outcomes_dta = graphs.driving_all_data(dta, mp, DOPRINT);
        outcomes_mod = stats.compute_outcomes(mp, s, sol_m);
    end
    
    % --- data --- 
    figure();
    abar = max([mp.abar_j0{:}]);
    aa = (0:abar-1)';
    for j=1:mp.ncartypes
        vv = squeeze(outcomes_dta.vkt(:,j,:));
        subplot(2,2,j), plot(aa, vv, '-o');
        title(mp.lbl_cartypes{j});
        if j == 1
            legend(mp.lbl_types);
        end
        ylabel('Driving (1000 km/year)'); xlabel('Car age');
    end
    sgtitle('Data');
    
    % --- model ---
    figure();
    for j=1:mp.ncartypes
        aa = (0:mp.abar_j0{j}-1)'; 
        
        % vkt_all is (ntypes*ns); 
        % select (ntypes*abar_j) matrix of predicted driving 
        vv = outcomes_mod.vkt_all(:, s.ipt.car{j});
        
        % show 
        subplot(2,2,j), plot(aa, vv, '-o');
        title(mp.lbl_cartypes{j});
        if j == 1
            legend(mp.lbl_types);
        end
        ylabel('Driving (1000 km/year)'); xlabel('Car age');
    end
    sgtitle('Model');
end

function vkt_data = driving_all_data(dta, mp, DOPRINT)
  % drivin_all_data: helper function for graphs.fit_driving()
  % 
  % OUTPUT: 
  %   vkt_data: struct with two members
  %       vkt: ntypes*ncartypes*abar matrix of avg. driving 
  %       N: 
  
  if nargin < 3
      DOPRINT = true; 
  end
  
  vkt_data = struct();
  abar = max([mp.abar_j0{:}]);
  aa = (0:(abar-1))';
  vkt_data.vkt = nan(mp.ntypes, mp.ncartypes, abar);
  vkt_data.count = nan(mp.ntypes, mp.ncartypes, abar);
  
  h = waitbar(0, 'Computing vkt for all (tau,j,car age)');
  
  I0 = ~isnan(dta.vkt);
  for tau=1:mp.ntypes
      waitbar(tau/mp.ntypes, h, sprintf('Computing vkt for all (tau,j,car age), tau %d/%d', tau, mp.ntypes));
      I1 = dta.tau == tau & I0;
      dta_tau = dta(I1, :); 
      for j=1:mp.ncartypes
          
          % find tau-households driving a type j car
          I2 = (dta_tau.pt_car_type == j);
          
          for ia=1:abar
              age = aa(ia);
              
              I = I2 & dta_tau.pt_car_age == age;
              
              % how many are there?
              N = sum(dta_tau(I, :).count);
              vkt_data.count(tau, j, ia) = N; % store it
              
              % weighted avg. driving over the different car ages and different
              % ways of ending up driving a type j car (keeping a type j car,
              % replacing a type k with a type j car, etc... many different ways
              % of getting to drive a type j car)
              vkt_data.vkt(tau, j, ia) = (dta_tau(I,:).vkt' * dta_tau(I,:).count) / N;
          end
      end
  end
  
  close(h);
  N = nansum(vkt_data.count, 'all'); 
  avg_vkt = nansum(vkt_data.vkt .* vkt_data.count, 'all') / N;
  
  if DOPRINT
      fprintf('Overall avg. vkt in the data = %5.2f km/day\n', avg_vkt);
  end

end

function fit_keep0(mp, dat, ccp_tau)
    if nargin < 3
        % solve model here 
        [s, F, price_j, q_tau, q, ev_tau, ccp_tau, ccp_scrap_tau, ctp_tau, ed, h_tau, marketshares]=equilibrium.solve(mp);
    else 
        s = trmodel.index(mp, mp.abar_j0);
    end
    
    set(groot, 'defaultaxesfontsize', 14); 
    
    figure(); hold on; 
    abars = [mp.abar_j0{:}]; 
    yy_dat = nan(max(abars), mp.ntypes, mp.ncartypes); 
    yy_mod = yy_dat; 
    
    fprintf('Computing empirical keep probs '); 
    for j=1:mp.ncartypes
        for tau=1:mp.ntypes 
            I1 = s.is.car{j}; 
            I2 = s.id.keep; 
            yy_mod(:,tau,j) = ccp_tau{tau}(I1, I2); 
            
            for a=2:mp.abar_j0{j} % from 2: you cannot (should not) have the 1st age (0) incoming 
                I = dat.s_car_type == j & dat.tau == tau & dat.s_car_age == a; 
                if any(I) 
                    Ncell = sum(dat(I, :).count); 
                    Nkeep = sum(dat(I & dat.id == s.id.keep, :).count); 
                    yy_dat(a,tau,j) = Nkeep / Ncell; 
                end
            end
        end
        
        aa = (1:abars(j))'; 
        subplot(2,mp.ncartypes,j),              plot(aa, yy_dat(:,:,j)); 
        ylabel('Observed share who keep this car)'); title(['Data: ' mp.lbl_cartypes{j}]); 
        subplot(2,mp.ncartypes,j+mp.ncartypes), plot(aa, yy_mod(:,:,j)); 
        ylabel('Pr(keep|incoming car)'); title(['Model: ' mp.lbl_cartypes{j}]); 
        if j == 1 
            legend(mp.lbl_types); 
            title('Fit: keeping the existing car'); 
        end
        fprintf('%d/%d ', j, mp.ncartypes); 
    end
    hold off; 
    
    fprintf('Done! \n'); 
end

function fit_market_shares_agg(mp, s, sol_m, sol_d)
    % graphs.fit_market_shares(mp, s, sol_m, sol_d)  
    % market shares aggregated over consumer types
  
    s_model=sol_m.marketshares;
    s_data=sol_d.marketshares;

    ymax_d=ceil(100*max(max(sol_d.marketshares)))/100; 
    ymax_m=ceil(100*max(max(sol_m.marketshares)))/100; 

    % by HH types 
    graphs.myfigure;
    xlab = {mp.lbl_cartypes{:}, 'no car'}'; 

    ym = sum(sol_m.marketshares, 1); 
    yd = sum(sol_d.marketshares, 1); 
    bar([yd', ym']); 
    
    title('Aggregated market shares');
    set(gca, 'xticklabels', xlab);
    ylabel('Market share (%)');
    legend('Data', 'Model','Location','Northwest');
    graphs.set_fig_layout_post;

end

function fit_market_shares(mp, s, sol_m, sol_d)
    % graphs.fit_market_shares(mp, s, sol_m, sol_d)    

    ymax_d=ceil(max(max(sol_d.marketshares*100))); 
    ymax_m=ceil(max(max(sol_m.marketshares*100))); 

    % by HH types 
    graphs.myfigure;
    xlab = {mp.lbl_cartypes{:}, 'no car'}'; 
    for tau=1:mp.ntypes
        subplot(2,4,tau), bar(100*[sol_d.marketshares(tau, :)', sol_m.marketshares(tau, :)']); title(mp.lbl_types(tau));
        set(gca, 'xticklabels', xlab);
        xtickangle(90);
        ylabel('Market share (%)');
        ylim([0,max(ymax_m,ymax_d)]);
        graphs.set_fig_layout_post;
    end
    legend('Data', 'Model');

    % nicify 
    set(gca, 'xticklabels', xlab);
    ylabel('Market share');
    xlabel('Car type');
    legend('Data', 'Model','Location','Northwest');
    

end


function fit_market_shares_45(mp, s, sol_m, sol_d)
    % graphs.fit_market_shares(mp, s, sol_m, sol_d)    
    s_model=sol_m.marketshares*100;
    s_data=sol_d.marketshares*100;

    ymax_d=ceil(max(max(s_model))); 
    ymax_m=ceil(max(max(s_data))); 
    ymax=max(ymax_m, ymax_d);

    graphs.myfigure;
    xlab = {mp.lbl_cartypes{:}, 'no car'}'; 
    pl=plot(s_data, s_model, 'o');



    col={'r';'k';'g';'b';'r'};
    col={'#0072BD';'#D95319';'#EDB120';'#7E2F8E';'k'; '#77AC30';'#4DBEEE';'#A2142F'};
    col={'#EDB120';'#77AC30';'k';'g';'r'; '#0072BD';'#77AC30';'#4DBEEE';'#A2142F'};

    set(pl, 'MarkerSize',8, {'markerfacecolor'}, col(1:mp.ncartypes+1));
    xlabel('Actual market share (%)'); 
    ylabel('Predicted market share (%)'); 
    hold on
    plot(0:0.01:ymax, 0:0.01:ymax, '--k');
    legend(xlab, 'Location', 'SouthEast', 'EdgeColor',[1 1 1]);
    graphs.set_fig_layout_post;
    set(gca, 'fontname', graphs.fontname, 'Ygrid', 'off','Xgrid', 'off', 'box', 'on'); 
    xlim([0,ymax])
    ylim([0,ymax])
    axis square


end

function fit_purchases(mp,s, sol_m, sol_d, plottype, idx)
    % syntax graphs.fit_purchases(mp,s, sol_m, sol_d, plottype)

    graphs.myfigure;

    if nargin<5
      plottype='by car type';
    end

    if nargin<6
      idx=[s.is.car{:} s.is.nocar];
    end
    
    ccp_m=0; ccp_d=0; q_d=0;q_m=0; 
    for j=1:mp.ncartypes
        if strcmp(plottype, 'by car type')
            ccp_m=0; ccp_d=0; 
        end
        dtx=s.id.age(s.id.trade{j})';

        for tau=1:mp.ntypes

            q_tau_m= sol_m.q_tau{tau}(idx)*mp.tw(tau); 
            q_tau_d= sol_d.q_tau{tau}(idx)*mp.tw(tau); 

            if j==1
                q_m=q_m + q_tau_m; 
                q_d=q_d + q_tau_d; 
            end

            ccp_d=ccp_d+nansum(q_tau_d.*sol_d.ccp_tau{tau}(idx,s.id.trade{j}), 1)'; 
            ccp_m=ccp_m+nansum(q_tau_m.*sol_m.ccp_tau{tau}(idx,s.id.trade{j}), 1)'; 
        end

        if strcmp(plottype, 'by car type')
            subplot(2,2,j), bar(dtx, 100*[ccp_d/nansum(q_d) ccp_m/nansum(q_m)]); 
            title(mp.lbl_cartypes(j))
            xlabel('Age of purchased car')
            if j==mp.ncartypes
                legend({'Data', 'Model'})
            end
        elseif strcmp(plottype, 'all cars')
            if j==mp.ncartypes
                bar(dtx, 100*[ccp_d/nansum(q_d) ccp_m/nansum(q_m)]);
                xlabel('Age of purchased car')
                legend({'Data', 'Model'})
            end
        end
        xlim([-1,22]);
        ylabel('Choice probability (%)') ;
        graphs.set_fig_layout_post;
    end
end

function fit_scrap(mp,s, sol_m, sol_d, plottype)
    % syntax graphs.fit_purchases(mp,s, sol_m, sol_d, plottype)

    abarmax=23;
    graphs.myfigure;

    if nargin<5
      plottype='by car type';
    end

    dty_d=0; dty_m=0; q_m=0; q_d=0;
    for icar=1:mp.ncartypes
        idx=s.is.car_ex_clunker{icar}(1:21);
        dtx(1:numel(idx))=s.is.age(idx)'; % car ages 1 to abar-1 
    
         if strcmp(plottype, 'by car type')
              dty_d=0; dty_m=0; q_m=0; q_d=0;
         end

         for tau=1:mp.ntypes

              q_tau_j_m= sol_m.q_tau{tau}(idx).*mp.tw(tau); 
              q_tau_j_d= sol_d.q_tau{tau}(idx).*mp.tw(tau); 

              q_m= q_m+ q_tau_j_m; 
              q_d= q_d+ q_tau_j_d; 

              ccp_tau_j_m=(1-sol_m.ccp_tau{tau}(idx, s.id.keep)).*(sol_m.F.accprob(idx)+(1-sol_m.F.accprob(idx)).*sol_m.ccp_scrap_tau{tau}(idx));
              ccp_tau_j_d=(1-sol_d.ccp_tau{tau}(idx, s.id.keep)).*sol_d.ccp_scrap_tau{tau}(idx); 

              dty_d=dty_d+ccp_tau_j_d.*q_tau_j_d; 
              dty_m=dty_m+ccp_tau_j_m.*q_tau_j_m;
        end

        if strcmp(plottype, 'by car type')
            dims=[ceil(sqrt(mp.ncartypes)), round(sqrt(mp.ncartypes))]; 

            dty_m=100*dty_m./q_m;
            dty_d=100*dty_d./q_d;

            subplot(dims(1),dims(2), icar);

            hL= plot(dtx,dty_d, '-', dtx,dty_m, '-');
            set(hL, {'Marker'},{'o'; '+'}, 'LineWidth',2, 'MarkerSize',5);

            if icar==1; % only plot legend in first panel
                l1=legend('Data', 'Model');
                set(l1,'EdgeColor',[1 1 1],'Location','Northwest');
            end

            title(mp.lbl_cartypes{icar});
            if icar/2~=round(icar/2)
                ylabel('Scrap prob. (%)');
            end
            if icar>2
                xlabel('Age of car');
            end
            axis([0 abarmax 0 40]);
            graphs.set_fig_layout_post(); 

        end
    end

  if ~strcmp(plottype, 'by car type');
     dty_m=100*dty_m./q_m;
     dty_d=100*dty_d./q_d;

     hL= plot(dtx,dty_d, '-', dtx,dty_m, '-');
         set(hL, {'Marker'},{'o'; '+'}, ...
             'LineWidth',2, 'MarkerSize',5);
     l1=legend('Data', 'Model');
     set(l1,'EdgeColor',[1 1 1],'Location','Northwest');
     ylabel('Scrap Probability');
     xlabel('Age of car');
     axis([0 abarmax 0 40]);
     graphs.set_fig_layout_post(); 

  end
end


function fit_keep(mp,s, sol_m, sol_d, plottype)
    % syntax graphs.fit_purchases(mp,s, sol_m, sol_d, plottype)

    abarmax=21;
    graphs.myfigure;

    if nargin<5
      plottype='by car type';
    end

    dty_d=0; dty_m=0; q_m=0; q_d=0;
    for icar=1:mp.ncartypes
        idx=s.is.car_ex_clunker{icar}(1:21);
        dtx(1:numel(idx))=s.is.age(idx)'; % car ages 1 to abar-1 
    
         if strcmp(plottype, 'by car type')
              dty_d=0; dty_m=0; q_m=0; q_d=0;
         end

         for tau=1:mp.ntypes

              q_tau_j_m= sol_m.q_tau{tau}(idx).*mp.tw(tau); 
              q_tau_j_d= sol_d.q_tau{tau}(idx).*mp.tw(tau); 

              q_m= q_m+ q_tau_j_m; 
              q_d= q_d+ q_tau_j_d; 

              ccp_tau_j_m=sol_m.ccp_tau{tau}(idx, s.id.keep);
              ccp_tau_j_d=sol_d.ccp_tau{tau}(idx, s.id.keep); 

              dty_d=dty_d+ccp_tau_j_d.*q_tau_j_d; 
              dty_m=dty_m+ccp_tau_j_m.*q_tau_j_m;
        end

        if strcmp(plottype, 'by car type')
            dims=[ceil(sqrt(mp.ncartypes)), round(sqrt(mp.ncartypes))]; 

            dty_m=100*dty_m./q_m;
            dty_d=100*dty_d./q_d;

            subplot(dims(1),dims(2), icar);

            hL= plot(dtx,dty_d, '-', dtx,dty_m, '-');
            set(hL, {'Marker'},{'o'; '+'}, 'LineWidth',2, 'MarkerSize',5);

            if icar==1; % only plot legend in firs'Data',t panel 
                l1=legend('Data', 'Model');               
                set(l1,'EdgeColor',[1 1 1],'Location','Northeast');
            end

            title(mp.lbl_cartypes{icar});
            if icar/2~=round(icar/2)
                ylabel('Keep prob. (%)');
            end
            if icar>2
                xlabel('Age of car');
            end
            axis([0 abarmax 40 100]);
            graphs.set_fig_layout_post(); 

        end
    end

  if ~strcmp(plottype, 'by car type');
     dty_m=100*dty_m./q_m;
     dty_d=100*dty_d./q_d;

     hL= plot(dtx,dty_d, '-', dtx,dty_m, '-');
         set(hL, {'Marker'},{'+'; 'o'}, ...
             'LineWidth',2, 'MarkerSize',5);
     l1=legend('Data','Model');
     set(l1,'EdgeColor',[1 1 1],'Location','Northeast');
     ylabel('Keep probability (%)');
     xlabel('Age of car');
     axis([0 abarmax 40 100]);
     graphs.set_fig_layout_post(); 

  end
end


function fit_scrap_by_car_type(mp, s, sol_m, sol_d, plottype)
    % graphs.fit_scrap(mp, s, sol_m, sol_d) 

    if nargin<5
      plottype='conditional';
    end

    if strcmp(plottype, 'conditional')
      str_title0= 'Scrap probability conditional on trading';
    else
      str_title0= 'Unconditional scrap probability';
    end

    dims=[ceil(sqrt(mp.ncartypes)), round(sqrt(mp.ncartypes))]; 
    graphs.myfigure(); 
    tickstep=2;

    ccp_tau=sol.ccp_tau; 
    ccp_scrap_tau=sol.ccp_scrap_tau;
  
    for icar=1:mp.ncartypes;
      subplot(dims(1),dims(2), icar);
      str_title= sprintf('%s (%s)',str_title0, mp.lbl_cartypes{icar});
      
      set(gcf,'Color',[1 1 1]);
      colormap(summer);
      hold on;
      abarmax=s.abar_j{icar}; 
      
      for t=1:mp.ntypes;
      hold on;
      % probability of scrapping the existing car by consumers
      dty=nan(1, abarmax);
      dtx=nan(1, abarmax);
      if strcmp(plottype, 'conditional')
        dty(1:numel(s.is.car_ex_clunker{icar}))=ccp_scrap_tau{t}([s.is.car_ex_clunker{icar}]); 
      else
        dty(1:numel(s.is.car_ex_clunker{icar}))=ccp_scrap_tau{t}([s.is.car_ex_clunker{icar}]) ... 
                .* (1-ccp_tau{t}([s.is.car_ex_clunker{icar}],s.id.keep)); 
      end
      dtx(1:numel(s.is.car_ex_clunker{icar}))=s.is.age(s.is.car_ex_clunker{icar}); % car ages 1 to abar-1 
      % h=stairs([dtx,dtx(end)+1]-.5,[dty,dty(end)],'LineWidth',2);
      h=plot(dtx,dty,'LineWidth',2,'Marker','o','MarkerSize',5,'MarkerFaceColor','auto');

      if (~strcmp(mp.lbl_cartypes{icar},''))
      set(h,'DisplayName',sprintf('%s',mp.lbl_types{t}));
      else
      set(h,'DisplayName',sprintf('%s, %s',mp.lbl_types{t}, mp.lbl_cartypes{icar}));
      end
    end
    xlim([0,abarmax-1]);
    % tickstep=2;
    % xticks=dtx(1:tickstep:numel(dtx));  
    % ticklabels=compose('%d',xticks);
    title(str_title);
    l1=legend('Location','Southwest');
    ylabel(sprintf('Probability of scrap (%s)', mp.lbl_cartypes{icar}));
    xlabel(sprintf('Age of car (%s)', mp.lbl_cartypes{icar}));
    % set(gca,'XTick',xticks,'XTickLabel',ticklabels,'YGrid','on');
    set(gca,'YGrid','on');
    hold off;
  end


    % aggregated 
    graphs.myfigure;
    ym = sum(sol_m.marketshares, 1); 
    yd = sum(sol_d.marketshares, 1); 
    bar([ym', yd']); 
    
    % nicify 
    title('Aggregated market shares');
    set(gca, 'xticklabels', xlab);
    xtickangle(90);
    ylabel('Market share');
    legend('Data', 'Model');
end

function car_hist(plots, mp, s, dta,sol)
  % data.car_hist: Visualize car distribution using histograms
  %
  % SYNTAX graphs.car_hist(plots, mp, s, dta, sol)
  %
  % INPUTS: 
  %     plots:  cell array with plotnames
  %             available plots:
  %             plots={'carage_by_type',
  %                    'carcage',
  %                    'carcage_by_year',
  %                    'purchases',
  %                    'purchases_by_year',
  %                    'scrap','scrap_by_type'};
  %             if plots=[] all available plots are made
  %
  %     mp:     structure with model parameters
  %     s:      a struct with the following elements:
  %     dta:    state and choice data with variables:  

  numCarAges=numel(unique(dta.s_car_age(dta.s_car_type > -1,:)));
  years=unique(dta.year); 
  T=numel(years);

  if isempty(plots)
    plots={'carage_by_type', 'carcage', 'carcage_by_year', 'purchases', 'purchases_by_year', 'scrap','scrap_by_type'};
  end

  for i=1:numel(plots)
    graphs.myfigure;
    switch plots{i}
      case 'carage_by_type'
        nC=numel(unique(dta.s_car_type(dta.s_car_type>0)));
        for j=1:nC   
          out=data.tabulate([dta], {'s_car_age'}, {}, (dta.s_car_type==j));
          subplot(2,2,j), bar(out.s_car_age, out.share); xlabel('car age'); 
          ylabel('frequency'), title(mp.lbl_cartypes{j})
        end
      case 'carcage'
        out=data.tabulate([dta], {'s_car_age'}, {}, (dta.s_car_type>-1)); % select population who owns a car
        bar(out.s_car_age, out.share); xlabel('car age'); ylabel('frequency'), title('car age distribution');
      case 'carcage_by_year'
        out = nan(numCarAges, T); 
        for t=1:T
          out_t=data.tabulate([dta], {'s_car_age'}, {}, dta.year == years(t) & dta.s_car_age >= 0); % select population who owns a car
          out(:, t) = (out_t.share)' ; 
        end
        surf(years, 1:numCarAges-1, out(2:end, :)); 
        view(-240, 33) 
        ylabel('Car age'); xlabel('Year'); title(sprintf('Car stock, %d - %d', years(1), years(end)))
      case 'scrap'
        out=data.mean(dta, @(data) data.d_scrap, {'s_car_age'}, dta.s_car_age>0)
        bar(out.s_car_age, out.mean); xlabel('car age'); ylabel('frequency'), title('scrappage');
      case 'scrap_by_type'
        for j=1:mp.ncartypes  
          out=data.mean(dta, @(data) data.d_scrap, {'s_car_age'}, (dta.s_car_type==j));
          AX{j}=subplot(2,2,j); 
          bar(out.s_car_age, out.mean); 
          xlabel('car age'), ylabel('scrap, %'), title(mp.lbl_cartypes{j});
          allYLim(j,:) = ylim;
        end
        for j=1:mp.ncartypes 
          set(AX{j}, 'YLim', [min(allYLim(:,1)), max(allYLim(:,2))]);
        end
      case 'scrap-nokeep'
        out=data.mean(dta, @(data) data.d_scrap, {'s_car_age'}, dta.id>1)
        bar(out.s_car_age, out.mean); xlabel('car age'); ylabel('frequency'), title('scrappage (conditional on not keeping)');
      case 'purchases'
        out=data.tabulate([dta], {'d_car_age'}, {}, (dta.d_car_age>=0  & dta.id~=s.id.keep)); % select population who buys a car 
        bar(out.d_car_age, out.share); xlabel('car age'); ylabel('frequency'), title('car purchases');
      case 'purchases_by_year'
        out = nan(numCarAges+1, T); 
        for t=1:T
          out_t=data.tabulate([dta], {'d_car_age'}, {}, dta.year == years(t) & dta.d_car_age>=0 & dta.id~=s.id.keep); % select population who owns a car
          out(:, t) = (out_t.share)' ; 
        end
        bar3(out(1:end, :)); 
        set(gca,'XTick',1:numel(years),'XTickLabel',compose('%d',[years]));
        set(gca,'YTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        view(-240, 33) 
        ylabel('Age of purchased car'); xlabel('Year'); title('car purchases')
      case 'purchases_by_is'
        out = nan(numCarAges+1, 20); 
        for t=1:20
          out_t=data.tabulate([dta], {'d_car_age'}, {}, dta.s_car_age == t  & dta.d_car_age>=0 & dta.id~=s.id.keep); % select population who owns a car
          out(:, t) = (out_t.share)' ; 
        end
        bar3(out(1:end, :)); 
        set(gca,'XTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        set(gca,'YTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        view(135, 45) 
        ylabel('Age of purchased car'); xlabel('Age of incoming car'); title('Car purchases -- Data')
        if (nargin > 3)
        graphs.myfigure;
        out = nan(numCarAges+1, 20); 
        for t=1:20
          wcprobs=zeros(numCarAges+1,1);
          qa=0;
          for tau=1:mp.ntypes
            for ct=1:mp.ncartypes
              cprobs=zeros(numCarAges+1,1);
              for ctp=1:mp.ncartypes
                 cprobs=cprobs+sol.ccp_tau{tau}(s.is.car{ct}(t),s.id.trade{ctp}(1:numCarAges+1))';
              end
              wcprobs=wcprobs+cprobs*sol.q_tau{tau}(s.is.car{ct}(t))*mp.tw(tau);
              qa=qa+sol.q_tau{tau}(s.is.car{ct}(t))*mp.tw(tau);  
            end
          end
          wcprobs=wcprobs/qa;
          wcprobs=wcprobs/sum(wcprobs);
          out(:,t)=wcprobs;
        end
        bar3(out(1:end, :)); 
        set(gca,'XTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        set(gca,'YTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        view(135, 45) 
        ylabel('Age of purchased car'); xlabel('Age of incoming car'); title('Car purchases -- Model')

        graphs.myfigure;
        tau=1;
        ct=3;
        out = nan(numCarAges+1, 20); 
        for t=1:20
          out_t=data.tabulate([dta], {'d_car_age'}, {}, dta.tau==tau & dta.is == s.is.car{ct}(t) & dta.s_car_type==ct & dta.d_car_type==ct & dta.d_car_age>=0 & dta.id~=s.id.keep); % select population who owns a car
          ss=size(out_t.share,1);
          out(1:ss, t) = (out_t.share); 
        end
        bar3(out(1:end, :)); 
        set(gca,'XTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        set(gca,'YTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        view(135, 45) 
        ylabel('Age of purchased car'); xlabel('Age of incoming car'); title(sprintf('Car purchases of car type %i by households of type %i holding car type %i  -- Data',ct,tau,ct))
        graphs.myfigure;
        out = nan(numCarAges+1, 20); 
        for t=1:20
          cprobs=sol.ccp_tau{tau}(s.is.car{ct},s.id.trade{ct}(1:numCarAges+1))';
          wcprobs=wcprobs/sum(wcprobs);
          out(:,t)=wcprobs;
        end
        bar3(out(1:end, :)); 
        set(gca,'XTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        set(gca,'YTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        view(135, 45) 
        ylabel('Age of purchased car'); xlabel('Age of incoming car'); title(sprintf('Car purchases of car type %i by households of type %i holding car type %i  -- Model',ct,tau,ct))
        end
      case 'posttrade_by_year'
        out = nan(numCarAges+1, T); 
        for t=1:T
          out_t=data.tabulate([dta], {'d_car_age'}, {}, dta.year == years(t) & dta.d_car_age>=0); % select population who owns a car
          out(:, t) = (out_t.share)' ; 
        end
        bar3(out(1:end, :)); 
        view(75, 40) 
          % tick labels
        set(gca,'XTick',1:numel(years),'XTickLabel',compose('%d',[years]));
        set(gca,'YTick',1:2:numel(out_t.d_car_age),'YTickLabel',compose('%d',[out_t.d_car_age(1):2:out_t.d_car_age(end)]));
        ylabel('Age of car'); xlabel('Year'); title('post trade holdings distribution')
      otherwise
        fprintf('requested plot does not exist')
    end
  end
end

function pp = compute_hom_equilibria(mp)
    % For a 2-type economy, this solves the two 1-type equilibria and
    % returns just the two abar price-vectors. 
    %
    % OUTPUT: 
    %   pp: (ncarages*ntypes) matrix of equilibrium prices (including the
    %   exogenous new car price) for a homogenous equilibrium consisting
    %   only of consumers of type tau (in columns). 

    mp_ = mp; % backup for safe-keeping 
    
    % hard-coded list of fields that have to be updated. 
    % NOTE: there could be more type-specific fields, just 
    % not in what we are using in the paper currently 
    vv = {'u_0', 'u_a', 'mum', 'u_even', 'u_a_sq', 'lbl_types'};
    
    assert(strcmp(mp.modeltype, 'reducedform'), 'Only implemented for modeltype == "reducedform"'); 
    assert(mp.ntypes == 2, 'Only implemented for examples with two household types');
    assert(mp.ncartypes == 1, 'Only implemented for one-car settings');
    
    % preallocate output 
    pp = nan(mp.abar_j0{1}, mp.ntypes); 
    
    for tau = 1:2
        
        % copy model parameters 
        mp = mp_; 
        mp.ntypes = 1;
        mp.tw = [1];
        
        % overwrite type-specific coefficients 
        for i=1:numel(vv)
            v = vv{i};
            mp.(v) = mp_.(v)(tau,:);
        end
        
        % solve model 
        mp = trmodel.update_mp(mp);
        sol_homo = equilibrium.solve(mp);

        % store equilibrium prices from the solution 
        pp(:, tau) = [mp.pnew{1}; sol_homo.p];
    end

end % compute_hom_equilibria()

end % end of methods
end % end of data
