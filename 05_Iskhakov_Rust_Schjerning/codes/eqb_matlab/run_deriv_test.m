% run_test.deriv
clear all
addpath('matlabinclude');
addpath('autotrade');
mp=setparams.default(); % standard model parameters

% parameter adjustments
mp.es=1;
mp.u_og={0};
mp.sigma_s=1;
mp.ntypes=2; 
mp.tw=ones(mp.ntypes,1)/mp.ntypes;
mp.ncartypes=2;
mp.abar_j0={25, 25};
mp.nocheck_aggderivs=0; % set this to 1 to avoid doing checks of aggregate holdings distribution, prices, excess demand, etc

% parameters to compute derivatives wrt
mp.pnames={'u_a', 'u_0', 'mum', 'transcost', 'ptranscost', 'tc_sale', 'tc_sale_age', 'tc_sale_even',  'ptc_sale', 'acc_0', 'acc_a'}; 
%mp.pnames={'bet','mum', 'acc_0','acc_a','acc_even'};

% specification and stating values (elemnts in mp.p0 must conform with mp.pnames)
mp.p0.u_0={60; 60}; 	   % consumer spefic baseline utility of cars
mp.p0.u_a={-2.1,-2.1}; % car spefic slope coefficients on car agíng
mp.p0.mum={1;1};	   % consumer spefic marginal utility of money
mp.p0.pnew_notax=mp.pnew_notax; % car specific new car prices before taxes
mp.p0.tc_sale=3;
mp.p0.ptc_sale=.05;
mp.p0.tc_sale_age=1;
mp.p0.tc_sale_even=2;
mp.p0.acc_0=mp.acc_0;
mp.p0.acc_a=mp.acc_a;
mp.p0.acc_even=mp.acc_even;
mp.p0.acc_even{1}=-.03;
mp.p0.acc_even{2}=.03;
mp.p0.bet=0.95;
mp.p0.acc_0{2}=-2.9;
mp.p0.acc_a{2}=.01;
mp.p0.transcost=mp.transcost;
mp.p0.ptranscost=mp.ptranscost;

[pvec0,mp] = estim.mp2pvec(mp,mp.pnames);

mp=estim.update_mp(pvec0,mp);

tests.run(mp);





