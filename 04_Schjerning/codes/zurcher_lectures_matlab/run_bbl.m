% run_bbl: Demonstrate BBL method using Rust's engine replacement model Rust(Ecta, 1987) 

% This demonstration is loosely based on BBL estimation method, but:
% 1. uses actual solution of the model at true parameter values for the policy 
%    functions, instead of estimating these in the first step
% 2. uses a particular scheme to randomize policy when forming the moment inequalities
%    (instead of perturbing the replacement threshold point as in BBL, it adds shocks 
%     to all decision specific values and computes the perturbed policy)

clear all
clc
% close all

%% BBL is sensitive to exact implementation: we several types of perturbations
PolicyPerturbationType = 3; %1, 2 or 3
PolicyPerturbationSigma = 0.0005; %amount of perturbation
nPicGrid = 20; %number of grid points in the surface plot in each dimension

%% Set up the model
mpopt.integrated=1; %IMPORANT to use integrated value functions
mp=zurcher.setup(mpopt);
% Transition matrix for mileage
P = zurcher.statetransition(mp);
% Pay off
u = zurcher.u(mp);
% Bellman operator (in integrated values space)
bellman= @(V) zurcher.bellman(V, mp, u, P);

% Solve the model using poly-algorithm (use combination of SA and NK)
V0=0; % initial guess on fixed point W
[~,ccp,~,~]=dpsolver.poly(bellman, V0, mp.ap, mp.beta); 
% ccp is "first stage estimated" probability of keeping

% Solve for the stationary distribution over state space
pp = zurcher.eqb(mp,P,ccp);

% Perturb the "first stage estimated" policy
if PolicyPerturbationType==1
  npert=1;
  pert{1} = ccp(mp.n:-1:1);
elseif PolicyPerturbationType==2
  npert=20; %number of perturbations
  rng(123,'twister');
  Vccp = npl.phi(mp,ccp,P); %values from CCPs
  for j=1:npert
    Vccp1 = Vccp .* (1+randn(mp.n,1)*PolicyPerturbationSigma);% add noise 
    pert{j} = npl.lambda(Vccp1,mp,P); %back to policy
  end
elseif PolicyPerturbationType==3
  npert=40; %number of perturbations
  for j=1:npert
    pert{j} = ccp; %copy the ccp policy
    jj=50+2*j; %threshold indexes  
    if mod(j,2)==0
      %step up from threshold half of the times
      pert{j}(jj+1:end) = ccp(jj+1:end) + (1-ccp(jj))/2;
    else
      %step down from threshold half of the times
      pert{j}(jj+1:end) = ccp(jj+1:end) - (1-ccp(jj))/2;
    end 
  end
else
  error 'Unknown PolicyPerturbationType'
end

% Plot the ccp and the perturbed choice probabilities
fig=figure('Color',[1 1 1]);
ax=axes('Parent',fig);
hold(ax);
for j=1:npert
  plot(ax,mp.grid,[ccp,pert{j}],'Linewidth',0.5,'color','black');
end
plot(ax,mp.grid,ccp,'Linewidth',2,'color','red');
set(ax,'XLim',[min(mp.grid),max(mp.grid)],'YLim',[0,1]);
title('Estimated and perturbed policy functions');
ylabel('Probability of keeping');
xlabel('State space');
drawnow %draw before going to heavier computation below

% Check BBL criterion
check=bbl.objective(ccp, pert, pp, P, mp, {'c','RC'}, [mp.c, mp.RC]');
fprintf('BBL criterion check: %1.5e (should not be different from zero)\n',check);

% Make a BBL criterion plot
pnames = {'c','RC'};
cgrid=linspace(1,5,nPicGrid)';
RCgrid=linspace(5,15,nPicGrid)';
misurf=zeros(nPicGrid,nPicGrid);

% loop over c and RC values
for i=1:nPicGrid;
  fprintf('%2d ',i);
  for j=1:nPicGrid;
    theta = [cgrid(i), RCgrid(j)]';
    misurf(j,i)=bbl.objective(ccp, pert, pp, P, mp, pnames, theta);
    fprintf('.');
  end
  fprintf('\n');  
end

fig=figure('Color',[1 1 1]);
ax=axes('Parent',fig);
surf(ax,cgrid,RCgrid,misurf);
axis('tight');
title('Moment inequality estimation objective, bus problem');
ylabel('RC parameter');
xlabel('c parameter');
hold on;
plot3(ax,[mp.c,mp.c],[mp.RC,mp.RC],ax.get('ZLim'),'-m','Linewidth',2);
hold off;
