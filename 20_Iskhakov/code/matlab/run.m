% Run script to demonstrate the properties of models with consumption and savings
% as well as polyline class. To be used as starting point for working with this models code.
% Written by Fedor Iskhakov, Australian National University, 2016

close all
clear
clear classes % update classes in the memory
addpath('utils');

%% Polyline demo

for j=1:10
  p(j)=polyline(sort(rand(10,1)),rand(10,1));
end
p.plot;
env=p.upper_envelope(true);
h=env.plot(gca);
set(h,'Color','red','LineWidth',2);

p1=polyline(rand(10,1),rand(10,1));
p1.plot;
env1=p1.secondary_envelope;
h=env1.plot(gca);
set(h,'Color','red','LineWidth',2);

p.demo;
methods(p);

fprintf('Press any key to continue..');pause;fprintf('\n\n');


%% Flat simulated consumption path using Deaton model

m1=model_deaton;
m1.label='Deaton model solved with VFI';
m1.df=1/(1+m1.r);
m1.sigma=0;
m1.init=[30 35];
fprintf('Solving %s with value function iterations:\n',m1.label)
tic;
m1.solve_vfi;
t=toc;
fprintf('VFI solver done in %s\n',ht(t))

m2=model_deaton;
m2.label='Deaton model solved with EGM';
m2.df=1/(1+m2.r);
m2.sigma=0;
m2.init=[30 35];
nn=100;
fprintf('Solving %s with EGM %d times: ',m2.label,nn)
tic;
for i=1:nn
  m2.solve_egm;
end
t=toc;
fprintf('average run time is %s\n',ht(t/nn))

m1.nsims=100;
m2.nsims=100;
m1.sim;
m2.sim;
error1=mean((max(m1.sims.consumption')-min(m1.sims.consumption'))./mean(m1.sims.consumption'));
error2=mean((max(m2.sims.consumption')-min(m2.sims.consumption'))./mean(m2.sims.consumption'));
fprintf('Average relative error in consumption (min-max/mean over lifecycle, averaged over individuals):\n');
fprintf('VFI %1.15f\n',error1);
fprintf('EGM %1.15f\n',error2);
m2.plot('sim consumption');

m1.nsims=2;
m2.nsims=2;
m1.sim;
m2.sim;
m1.plot('sim consumption');
m2.plot('sim consumption');

fprintf('\nPress any key to continue..');pause;fprintf('\n\n');



%% Nice simulation graphics using Phelps model
m3=model_phelps;
m3.nsims=100;
m3.solve_egm;
m3.sim;
m3.plot('sim');
fprintf('Simulation plots for Phelps model produced\n')

fprintf('\nPress any key to continue..');pause;fprintf('\n\n');

%% Flat simulated consumption path using Retirement model
m4=model_retirement;
m4.ngridm=500;
m4.df=1/(1+m4.r); %flat consumption hopefully
m4.sigma=0;
m4.lambda=eps; %no EV taste shocks
m4.nsims=5;
m4.init=[5 20];
tic
m4.solve_dcegm;
t=toc;
fprintf('Retirement model solved with DC-EGM in %s\n',ht(t));
m4.plot('policy')
m4.plot('value')
m4.plot('prob_work')
m4.sim;
m4.plot('sim consumption');

fprintf('\nPress any key to continue..');pause;fprintf('\n\n');

%% Nice simulation graphics using retirement model
m5=model_retirement;
m5.ngridm=500;
m5.df=1/(1+m5.r); %flat consumption hopefully
m5.sigma=0.35;
m5.lambda=0.2; %some EV taste shocks
m5.nsims=50;
m5.init=[5 20];
tic
m5.solve_dcegm;
t=toc;
fprintf('Retirement model solved with DC-EGM in %s\n',ht(t));
m5.plot('policy');
m5.plot('value');
m5.plot('prob_work');
m5.sim;
m5.plot('sim');
fprintf('Simulation plots for retirement model produced\n')


