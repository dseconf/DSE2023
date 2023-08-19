% run_errorbound: Code to make plot of error_bound and log(error_bound) against iteration count
clear;
close all

% switches
do_nk=0; % 1 for Newton iterations, 0 for successive approximations

% Model parameters 
mp.bellman_type='iv';
mp.n=90;
mp=zurcher.setup(mp);

legends={};
betavec=[0.95 0.99 0.999 0.9999 0.99999];
for ibeta=1:numel(betavec)
	if ~do_nk % force dpsolver to do 10000 successive approximations
		mp.ap.sa_min=10000;
	end

	V0=0;
	mp.beta=betavec(ibeta);

  % belman equation
	bellman= @(V) zurcher.bellman(V, mp);
	[V, pk0, dV, iterinfo(ibeta)]=dpsolver.poly(bellman, V0, mp.ap, mp.beta);	
	legends{ibeta}=sprintf('beta=%g', betavec(ibeta)) ;

	tol=iterinfo(ibeta).sa.tol;
	if do_nk
			tol=[tol; iterinfo(ibeta).nk.tol];	
			legends_log_err{ibeta}=legends{ibeta};
	else
		y=log(tol); 
		y=y(~isinf(y));
		x=[ones(numel(y),1) (1:numel(y))'];
		b=inv(x'*x)*x'*y;
	  legends_log_err{ibeta}=sprintf('slope=%10.5f, beta=%-10g ', b(2),betavec(ibeta));
	end

	colorOrder = get(gca, 'ColorOrder');

	figure(1)
	hold on
	plot(1:numel(tol), (tol), 'Color', colorOrder(ibeta,:), 'LineWidth', 1.5);

	figure(2)
	hold on
	plot(1:numel(tol), log(tol), 'Color', colorOrder(ibeta,:),'LineWidth', 1.5);
	grid
end; 

figure(1)
set(gca,'FontSize',20)
xlabel('Iteration count');
ylabel('Error bound');
legend(legends, 'Location', 'SouthEast')
title('Error bound vs iteration count');
xlim([(0+mp.ap.sa_min*do_nk-2) mp.ap.sa_min-1+5*do_nk])

figure(2)
set(gca,'FontSize',20)
xlabel('Iteration count');
ylabel('log(Error bound)');
legend(legends_log_err, 'Location', 'SouthEast')
title('Log of error bound vs iteration count');
xlim([(0+mp.ap.sa_min*do_nk-2) mp.ap.sa_min-1+5*do_nk])
