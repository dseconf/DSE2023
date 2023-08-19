% run_fxp: Solve Rust's engine replacement model Rust(Ecta, 1987) 
clear all
% clc

% algorithm switch 
algorithm = 'sa';   % sa or poly
algorithm = 'poly'; % sa or poly

% Read default parameters in to struct mp 
mp.n=90;

%% Parameters for solution algorithm (mp.ap used in dpsolver.m)
mp.ap.printfxp=2;	% (0= no output), (compressed output), (2= detailed iteraion output)
mp.ap.sa_max=500000;	% Set minimum number of contraction steps (successive approximations)
mp.ap.sa_min=5;	% Set minimum number of contraction steps (successive approximations)

mp.beta=0.9999;
mp.bellman_type='ev';  	% bellman in expected value ('ev') or ('iv') integraded value function space  

mp=zurcher.setup(mp);

% Transition matrix for mileage
[P] = zurcher.statetransition(mp);
[u] = zurcher.u(mp);

V0=0; % Initial guess on fixed point V0
% bellman equation
bellman= @(V) zurcher.bellman(V, mp, u, P);    		
switch algorithm
	case 'sa' % solve by successive approximations (SA)
		[W, iterinfo]=dpsolver.sa(bellman, V0, mp.ap);	
	case 'poly' % solve using poly-algorithm (use combination of SA and NK)
		[W, pk, dV, iterinfo]=dpsolver.poly(bellman, V0, mp.ap, mp.beta);	
	otherwise
 		error('Algorithm must be ''sa'' or ''poly''');
end

% recall that bellman either integrated value, V or Expected value ev 
if strcmp(mp.bellman_type, 'iv');
	ev=P{1}*W;
else
	ev=W;
end


