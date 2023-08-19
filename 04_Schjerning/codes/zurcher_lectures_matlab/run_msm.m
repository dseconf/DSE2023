% run_bbl: Demonstrate MSM method using Rust's engine replacement model Rust(Ecta, 1987) 

% This demonstration is loosely based on MSM estimation method, but:
% 1. uses actual solution of the model as data moments, instead of estimating these in the first step
% 2. uses a particular scheme to randomize policy when forming the moment inequalities
%    (instead of perturbing the replacement threshold point, it adds EV(1) shocks to all
%     decision specific values and computes the perturbed policy by taking the max)

clear all
clc
% close all

useSolutionForSims = 0; %use moments computed from solution instead of simulations
updateWeightingMatrix = 1; %use updated weighting matrix
nPicGrid = 20; %number of grid points in the surface plot in each dimension

%% Set up the model
mpopt.bellman_type='ev'; %use expected value functions
mp=zurcher.setup(mpopt);

%additional setting for MSM
mp.useSolutionForSims=useSolutionForSims; 
mp.nsims = 1000; %number of buses to simulate
mp.Tsims = 1000; %number of time periods to simulate
mp.momentScale = 1e2; %scale of the moments

% Bellman operator (in integrated values space)
bellman= @(V) zurcher.bellman(V, mp);

% Solve the model using poly-algorithm (use combination of SA and NK)
V0=0; % initial guess on fixed point W
[ev,pk,~,~]=dpsolver.poly(bellman, V0, mp.ap, mp.beta); 

mp.ev_true=ev;
mp.pk_true=pk;

% solve for the stationary distribution over state space
P = zurcher.statetransition(mp);
data_moments = zurcher.eqb(mp,P,pk)*mp.momentScale; %scale
data_moments(end) = []; %drop the last moment as fraction sum up to one
% data_moments is "computed fractions of data in mileage bins"
% we use actual stationary distribution computed from the solution of 
% the model to represent data, so this is the best case scenario for MSM

% Weighting matrix
if useSolutionForSims
  % using solution of the model instead of simulations
  W=eye(mp.n-1);
  Wlabel='computed moments';
elseif ~updateWeightingMatrix
  % simulations with identity weighting matrix
  W=eye(mp.n-1);
  Wlabel='identity weighting matrix';
else
  % updated (second step) weighting matrix
  W=eye(mp.n-1); % first step weighting
  % use true values as consistent estimates from the first step
  est1c = mp.c; 
  est1RC = mp.RC;
  % compute the moment differences individually for each bus with W=eye
  [~,mom,mom_indiv]=msm.objective(data_moments, W, mp, mp.ap, {'c','RC'}, [est1c, est1RC]');
  C = cov(mom_indiv); %covariance matrix of the moment conditions
  cC=cond(C); %condition number
  if cC>1e5
    fprintf('Condition number of covariance matrix is %1.3e, switching to diagonal W\n',cC) 
    W = diag(1./var(mom_indiv)); % if convariance matrix is badly conditioned, use diagonal
  else
    W = inv(cov(mom_indiv));
  end
  Wlabel='updated weighting matrix';
end

if ~useSolutionForSims
  % Plot the data vs simulated moment at true parameter
  [~,sims_moments]=msm.objective(data_moments, W, mp, mp.ap, {'c','RC'}, [mp.c, mp.RC]');
  fig=figure('Color',[1 1 1]);
  ax=axes('Parent',fig);
  hold(ax);
  plot(ax,mp.grid(1:end-1),data_moments,'Linewidth',2,'color','black');
  stairs(ax,mp.grid(1:end-1),sims_moments,'Linewidth',2,'color','red');
  set(ax,'XLim',[min(mp.grid(1:end-1)),max(mp.grid(1:end-1))]);
  title(sprintf('Solution-based and simulated moments at the true parameter (%s)',Wlabel));
  ylabel('Fraction of fleet');
  xlabel('Mileage');
  drawnow %draw before going to heavier computation below
end

% Make a MSM criterion plot
pnames = {'c','RC'};
cgrid=linspace(2,3,nPicGrid)';
RCgrid=linspace(11,12,nPicGrid)';
msmsurf=zeros(nPicGrid,nPicGrid);

% loop over c and RC values
for i=1:nPicGrid;
  fprintf('%2d ',i);
  for j=1:nPicGrid;
    theta = [cgrid(i), RCgrid(j)]';
    msmsurf(j,i)=msm.objective(data_moments, W, mp, mp.ap, pnames, theta);
    fprintf('.');
  end
  fprintf('\n');  
end

fig=figure('Color',[1 1 1]);
ax=axes('Parent',fig);
surf(ax,cgrid,RCgrid,msmsurf);
axis('tight');
title(sprintf('Method of simulated moments objective (%s)',Wlabel));
ylabel('RC parameter');
xlabel('c parameter');
hold on;
plot3(ax,[mp.c,mp.c],[mp.RC,mp.RC],ax.get('ZLim'),'-m','Linewidth',2);
hold off;

% find optimal parameters
% x=fminsearch(@(theta) msm.objective(data_moments, W, mp, mp.ap, pnames, theta), [0;0])
% hold on;
% plot3(ax,[x(1),x(1)],[x(2),x(2)],ax.get('ZLim'),'-r','Linewidth',2);
