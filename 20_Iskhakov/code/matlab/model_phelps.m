% Matlab class to implement the Phelps model of consumption and savings
% with income shocks and credit constraint
% Written by Fedor Iskhakov, University of New South Wales, 2015
% See Iskhakov, Jorgensen, Rust and Schjerning "The Endogenous Grid Method for Discrete-Continuous Dynamic Choice Models with (or without) Taste Shocks" (QE, 2017)

classdef model_phelps < handle
% This class defines a Deaton consumption-savings model

properties (Access=public)
	%Default parameters of the model
	label		= 'Phelps''s model'; %name of this model
	Tbar		= 25		; %number of periods (fist period is t=1) 
	ngridm	= 100		; %number of grid points over assets
	mmax		= 75		; %maximum level of assets
	expn		=	10		; %number of quadrature points used in calculation of expectations
	nsims		= 10		; %number of simulations
	r 			= 0.05	; %interest rate
	df 			= 0.95	; %discount factor
	sigma   = 0.05	; %sigma parameter in logNormal shock to returns
	init    =[0 30] ; %interval of the initial wealth
end %properties

properties (SetAccess=private, GetAccess=public)
	%Entities computed inside of the model
	policy=polyline;
	value=polyline;
	sims=struct();
end %properties

methods (Access=public)

	%Definition of the model
	%Vectorize!!!
	function u=util(me,consumption) %utility
		u=log(consumption);
	end %util
	function mu=mutil(me,consumption) %marginal utility
		mu=1./consumption;
	end %mutil
	function cons=imutil(me,mutil) %inverse marginal utility
		cons=1./mutil;
	end %imutil
	function w1=budget(me,it,savings,shocks) 
		%wealth in period t+1, where it=t
		%inputs: savings = 1x(ngridm) row vector of savings
		%				 shocks = (expn)x1 column vector of shocks
		%output: w1 = (expn)x(ngridm) matrix of all possible next period wealths
		w1=(ones(size(shocks,1),1)*savings(1,:)) ...
			 .* (exp(me.r+shocks(:,1)-(me.sigma^2)/2)*ones(1,size(savings,2))) ...
			  + 1e-8; %need tiny positive number to avoid NaNs in solution
	end %budget
	function mw1=mbudget(me,it,savings,shocks) 
		%derivative of wealth in t+1 w.r.t. savings
		%inputs and outputs as above
		mw1=exp(me.r+shocks(:,1)-(me.sigma^2)/2)*ones(1,size(savings,2));
	end %mbudget

	%Solver EGM
	function solve_egm(me,compute_value)
		%solve the model with EGM algorithm
		if nargin<2
			compute_value=false;%do not compute value by default (not needed for EGM)
		end
		me.policy=polyline;%delete previous solution
		me.value=polyline;%delete previous solution
		[quadp quadw]=quadpoints(me.expn,0,1);	%quadrature points
		quadstnorm=norminv(quadp,0,1);	%normally distributed
		savingsgrid=linspace(0,me.mmax,me.ngridm); %grid over savings (should start with 0)
    %main EGM loop
    for it=me.Tbar:-1:1
	    if it==me.Tbar
        %terminal period
        me.policy(it)=polyline([0 me.mmax],[0 me.mmax],'Policy function in period T');
        if compute_value
	        me.value(it)=polyline([0 me.mmax],[0 NaN],'Value function in period T');
        	%vf(1)=0 is important, otherwise vf can not be computed in terminal period
        end
	    else
        %not the terminal period
        wk1=me.budget(it,savingsgrid,quadstnorm*me.sigma);%next period wealth matrix
        cons1=me.policy(it+1).interpolate(wk1);%next period consumption
        mwk1=me.mbudget(it,savingsgrid,quadstnorm*me.sigma);%next period marginal wealth
        rhs=quadw'*(me.mutil(cons1).*mwk1);%RHS of Euler equation
        cons0=me.imutil(me.df*rhs);%current period consumption
        me.policy(it)=polyline(savingsgrid+cons0,cons0,sprintf('Policy in periond %d',it));
      	me.policy(it)=me.policy(it).inset(0,0,0);%connect the dots
        if compute_value
	        ev=quadw'*me.value_function(it+1,wk1);%expected value function
	        me.value(it)=polyline(savingsgrid+cons0,me.util(cons0)+me.df*ev,sprintf('Value function in periond %d',it));
	        %add special first point: ev of saving zero
	        me.value(it)=me.value(it).inset(0,ev(1),0);%last zero to make first point
	      end
	    end % if(terminal period)
    end %it
	end

	%Solver VFI
	function solve_vfi(me)
		%solve the model with Value functions
		me.policy=polyline;%delete previous solution
		me.value=polyline;%delete previous solution
		[quadp quadw]=quadpoints(me.expn,0,1);	%quadrature points
		quadstnorm=norminv(quadp,0,1);	%normally distributed
		wealthgrid=exp(linspace(log(1e-4),log(me.mmax),me.ngridm));%log-grid over wealth
		opt=optimset('largescale','off','HessUpdate','bfgs', ...
								 'Display','none','TolFun',1e-8,'TolX',1e-8);
		function res=bellman(it,wealth,consumption)
			%the maximand of Bellman equation
			if consumption<wealth
				wk1=me.budget(it,wealth-consumption,quadstnorm*me.sigma);%next period wealth column-vector
				ev=quadw'*me.value(it+1).interpolate(wk1);%expectation over tomorrow's value
				res=me.util(consumption)+me.df*ev;
				res=-res;%for minimization
			else
				res=NaN;
			end
		end
    %main VFI loop
    for it=me.Tbar:-1:1
			fprintf('%2d ',it);
	    if it==me.Tbar
        %terminal period
        me.policy(it)=polyline(wealthgrid,wealthgrid,'Policy function in period T');
        me.value(it)=polyline(wealthgrid,me.util(wealthgrid),'Value function in period T');
        fprintf('%s\n',repmat('.',1,me.ngridm));
	    else
        %not the terminal period
        for i=1:me.ngridm
        	%loop over state space
        	if i==1
        		x0=wealthgrid(i)/2;
        	else
        		x0=me.policy(it).x(i-1);%use last found optimum as starting value
        	end
        	[xo,fo]=fminunc(@(c)bellman(it,wealthgrid(i),c),x0,opt);
        	me.policy(it)=me.policy(it).inset(wealthgrid(i),xo,i-1);%add point
        	me.value(it)=me.value(it).inset(wealthgrid(i),-fo,i-1);%minus due to minimization
        	fprintf('.');
        end %state space
        fprintf('\n');
	    end %if(terminal period)
    end %it
	end

	%Solver Euler
	function solve_euler(me)
		%solve the model by solving Euler euqations
		error 'Not implemented: homework'
	end

	%Simulator
	function sim(me,seed)
		%simulate from the model
		%input: init = initial wealth
		%			  seed = seed for random number generator 
		%							(to run identical or varying simulations)
		if me.policy(1).len<1
			error 'The model should be solved first'
		end
		%fix the stream of random numbers
		if nargin<2
			rng(7134,'twister');
		else
			rng(seed,'twister');
		end
		%allocate
		me.sims=struct('wealth0',nan(me.nsims,me.Tbar), ...
									 'wealth1',nan(me.nsims,me.Tbar), ...
									 'consumption',nan(me.nsims,me.Tbar), ...
									 'shock',nan(me.nsims,me.Tbar), ...
									 'return',nan(me.nsims,me.Tbar));
		%simulate choices and states
    for it=1:me.Tbar
    	if it==1
    		%draw initial wealth form uniform distribution on given interval
				me.sims.wealth0(:,it)=me.init(1)+rand(me.nsims,1)*(me.init(2)-me.init(1));
				me.sims.shock(:,it)=nan(me.nsims,1);
				me.sims.return(:,it)=nan(me.nsims,1);
			else
				%transition
				me.sims.shock(:,it)=norminv(rand(me.nsims,1),0,1)*me.sigma;%normal scaled
				me.sims.wealth0(:,it)=diag(me.budget(it-1,me.sims.wealth1(:,it-1)',me.sims.shock(:,it)));%match savings and shocks one-to-one
				me.sims.return(:,it)=exp(me.r+me.sims.shock(:,it)-(me.sigma^2)/2);
			end
			%choice
			me.sims.consumption(:,it)=me.policy(it).interpolate(me.sims.wealth0(:,it));
			%end of period wealth
			me.sims.wealth1(:,it)=me.sims.wealth0(:,it)-me.sims.consumption(:,it);
    end
	end

	%Plotting
	function ax=plot(me,what2plot)
		%plot computed entities: policy, value, "sims smth"
		if nargin<2
			what2plot='solution';
		end
		what2plot=strsplit0(lower(what2plot),' ');%explode by space, convert to cell array
		if sum(ismember(what2plot,{'policy','solution','pol'}))>0
			%plotting the policy functions
			if me.policy(1).len<1
				error 'Nothing to plot'
			end
			fig1 = figure('Color','white');
			ax = axes('Parent',fig1,'FontSize',14);
			me.policy.plot(ax,'Marker','none','LineWidth',1);
			hold(ax,'all');
			set(ax,'XLim',[0 me.mmax],'YGrid','on','XGrid','on');
			box(ax,'on');
			xlabel(ax,'Wealth','FontSize',14);
			title(ax,sprintf('%s: %s',me.label,'optimal consumption rules'),'FontSize',14);
		elseif sum(ismember(what2plot,{'value','value_function','val','valfunc','vf'}))>0
			%plotting value functions
			if me.value(1).len<2
				error 'Nothing to plot'
			end
			fig1 = figure('Color','white');
			ax = axes('Parent',fig1,'FontSize',14);
			%replace the analytical region with polylines
			k=100;%points to be added in analytical region
			[~,data]=me.value.chop(1); %through first points from all value functions
			for it=1:me.Tbar
				pt=exp(linspace(log(eps^.5),log(me.value(it).x(2)),k));%log-grid
				data(it)=data(it).grow(polyline(pt,me.value_function(it,pt)),true); %true for adding in the front
				data(it).label=sprintf('t=%d',it);
			end
			data.plot(ax,'Marker','none','LineWidth',1);
			hold(ax,'all');
			set(ax,'XLim',[0 me.mmax],'YGrid','on','XGrid','on');
			box(ax,'on');
			xlabel(ax,'Wealth','FontSize',14);
			title(ax,sprintf('%s: %s',me.label,'value functions'),'FontSize',14);
		elseif sum(ismember(what2plot,{'simulations','simulation','sim','sims'}))>0
			flds=fields(me.sims)';
			map=ismember(flds,what2plot);
			map=map | sum(map)==0; %plot everything if nothing is chosen
			for fld=flds(map)
				if numel(me.sims.(fld{1}))<1
					error 'Nothing to plot'
				end
				fig1 = figure('Color','white');
				ax = axes('Parent',fig1,'FontSize',14);
				h=plot([1:me.Tbar]',me.sims.(fld{1})');
				set(h,'Color','k','LineWidth',.5);
				hold(ax,'all');
				set(ax,'XLim',[1 me.Tbar],'YGrid','of','XGrid','of');
				%prevent Y scale from getting too detailed
				corr=@(x,tol) [mean(x)-max(max(x)-min(x),tol)/2,mean(x)+max(max(x)-min(x),tol)/2];
				set(ax,'YLim',corr(get(ax,'YLim'),1e-8));
				box(ax,'off');
				xlabel(ax,'Age','FontSize',14);
				title(ax,sprintf('%s: simulated %s',me.label,fld{1}),'FontSize',14);
			end
		else
			error 'Didn''t understand what to plot..'
		end
	end

	%Calculator of value functions that uses the analytical part in credit constrained region
	function res=value_function(me,it,x)
		%interpolates value function at period t=it using analytical part
		if me.value(it).len<2
			error(sprintf('Can not compute value function at period %d because it only has %d points',it,me.value(it).len))
		end
		res=nan(size(x)); %output of the same size as x
		mask=x<me.value(it).x(2); %all points in credit constrained region
		mask=mask | it==me.Tbar; %in the terminal period all points are in the constrained region
		res(mask)=me.util(x(mask))+me.df*me.value(it).y(1); %the first value in me.value is EV from zero savings!
		res(~mask)=me.value(it).interpolate(x(~mask));		
	end

end %methods
end %of classdef

