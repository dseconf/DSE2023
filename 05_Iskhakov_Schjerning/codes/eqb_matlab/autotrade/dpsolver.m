% This class is generic solver for Bellman equations
% in dynamic programming models using John's poly-algorithm
% (sa and newton are parts of this class as well)
% 
% July 2020
%

classdef dpsolver
  	methods (Static)
    	function ap=setup(apopt)
    		% SOLVE.SETUP: setup of algorithm solution parameters, ap,  used in solve
			%  syntax: ap=setup(apopt)
			%
			%  INPUT:  apopt (optional): If apopt is specified, default parameters will overwritten with elements in apopt.
			%
			%  OUTPUT: apopt: algorithm parameter structure
			%
			% See also:
			%   dpsolver.sa, dpsolver.nk

			% default values of ap
			ap.sa_max=			20;			% Maximum number of contraction steps
			ap.sa_min=			2; 			% Minimum number of contraction steps
			ap.sa_tol=			1.0e-3;	% Absolute tolerance before (in dpsolver.poly: tolerance before switching to N-K algorithm)
			ap.max_fxpiter= 	35;          % Maximum number of times to switch between Newton-Kantorovich iterations and contraction iterations.
			ap.pi_max=			40;			% Maximum number of Newton-Kantorovich steps
			ap.pi_tol=			1.0e-13;	% Final exit tolerance in fixed point algorithm, measured in units of numerical precision
			ap.tol_ratio=		1.0e-03;	% Relative tolerance before switching to N-K algorithm  
											% when discount factor is supplied as input in dpsolver.poly
			ap.printfxp=		0; 			% Print iteration info for fixed point algorithm
	    									% ap.printfxp=0 (No printing), ap.printfxp=1 (Summary info), ap.printfxp>1 (Detailed info)
			ap.check_deriv=	0;				% Check analytical derivatives of bellman against numerical derivatives 
			ap.check_deriv_tol=	1e-6;		% Absolute tolerance before for difference between analytical and numerical derivatives of bellman operator

 
	    	if nargin>0
	        	pfields=fieldnames(apopt);
	        	for i=1:numel(pfields);
	            	ap.(pfields{i})=apopt.(pfields{i});
	        	end
	    	end
		end % end of SOLVE.setup

		function [df_nm, df_an, df_err] =check_deriv(fun, x0, ap)
			[f_nm, ~, df_an]=fun(x0); 
			df_nm=gradp(fun,x0); 
			df_err=df_an-df_nm;
		end


		function [V, P, dV, iter]=poly(bellman, V0, ap, bet)
			% dpsolver.poly: Solve for fixed point using a combination of Succesive Approximations (SA) and Newton-Kantorovich (NK) iterations 
			%
			%  syntax:	[V, P, dV, iter]=dpsolver.poly(bellman, V0, ap, bet)
			%
			%  INPUT:
			%     bellman:  [V, P, dV]=Bellman equation with fixed point, V
			%	 							Matlab function on the form [V,dV]=bellman(V)
			%								where V is the value function (m x 1), P is the policy function 
			%								and dV is the (m x m) Frechet derivative of the Bellman operator]
			%     V=  			Initial guess value function,V.
			%								[m x 1 matrix]
			%
			%			ap:				Algorithm paramaters. See dpsolver.setup
			%
			%			bet:		  Discout factor. Enters rule for stopping SA and switching to NK iterations. 
			%								SA should stop prematurely when relative tolerance is close to bet. 
			%
			%  OUTPUT:
			%     V:         m x 1 matrix. Fixed point, V
			%			P:				 Policy function at fixed point
			%			dV:				 Frechet derivative of the Bellman operator

			%-----------------------------------------------------------------------------------------------------------------------------

			% Set default settings for fixed point algorithm, for ap's not given in input
			% (overwrites defaults if ap is given as input)
			if isempty(ap)
				ap=dpsolver.setup;
			else
				ap=dpsolver.setup(ap);
			end

			solutiontime=tic;
			for k=1:ap.max_fxpiter; %poly-algorithm loop (switching between SA and N-K and back)

			    % SECTION A: CONTRACTION ITERATIONS
			    if ap.printfxp>0
			        fprintf('\n');
			        fprintf('Begin contraction iterations (for the %d. time)\n',k);
			    end;
			    if nargin>3
						[V0, iter(k).sa]=dpsolver.sa(bellman, V0, ap, bet);			    	
			    else
						[V0, iter(k).sa]=dpsolver.sa(bellman, V0, ap);
  				end

			    % SECTION 2: NEWTON-KANTOROVICH ITERATIONS
			    if ap.printfxp>0
			        fprintf('\n');
			        fprintf('Begin Newton-Kantorovich iterations (for the %d. time)\n',k);
			    end
			    [V0, P, dV, iter(k).nk]=dpsolver.nk(bellman, V0, ap);	

			    if iter(k).nk.converged
			    	if ap.printfxp>0
			             fprintf('Convergence achieved!\n\n');
			             fprintf('Elapsed time: %3.5f (seconds)\n',toc(solutiontime));
			         end
			         break; %out of poly-algorithm loop
			     else
			         if k>=ap.max_fxpiter
			             warning('dpsolver.poly: Non-convergence! Maximum number of iterations exceeded without convergence!');
			             break; %out of poly-algorithm loop with no convergence
			         end
			    end
			end
			V=V0;
    	end % end of SOLVE.poly 

		function [V, iter]=sa(bellman, V0, ap, bet)
			% dpsolver.sa: Solve for fixed point using successive approximations
			%
			%  syntax:	[V, P, dV, iter] = dpsolver.sa(bellman, V0, ap, bet):
			%
			%  INPUT:
			%     bellman:  V=Bellman equation with fixed point, V
			%	 							Matlab function on the form V=bellman(V)
			%  							where V is the value function (m x 1) 
			%     V0:  			Initial guess value function,V.
			%								m x 1 matrix
			%			ap:				Algorithm paramaters. See dpsolver.setup
			%
			%			bet:		  Discout factor. Enters rule for stopping SA and switching to NK iterations. 
			%								SA should stop prematurely when relative tolerance is close to bet.  
			%
			%  OUTPUT:
			%     V:      m x 1 matrix. Approximation of fixed point, V


			% Set default settings for fixed point algorithm, for ap's not given in input
			% (overwrites defaults if ap is given as input)
			ap=dpsolver.setup(ap);

			solutiontime=tic;
			iter.tol=nan(ap.sa_max,1);
			iter.rtol=nan(ap.sa_max,1);
			iter.converged=false;
			for i=1:ap.sa_max; 
				V=bellman(V0);
				iter.tol(i)=max(max(abs(V-V0)));
				iter.rtol(i)=iter.tol(i)./iter.tol(max(i-1,1));
				V0=V; % accept SA step and prepare for new iteration
		   	% Stopping criteria
		   	if nargin >3
		    	if (i>=ap.sa_min) && (abs(bet-iter.rtol(i)) < ap.tol_ratio)
		      		iter.message='SA stopped prematurely due to relative tolerance. Start NK iterations';
		      		break
			  	end
			end

		    % Rule 2: 
        	adj=ceil(log10(abs(max(max(V0)))));
        	ltol=ap.sa_tol*10^adj;  % Adjust final tolerance
        	ltol=ap.sa_tol;
				if (i>=ap.sa_min) && (iter.tol(i) < ltol)
					iter.message=sprintf('SA converged after %d iterations, tolerance: %10g\n',i, iter.tol(i));
					iter.converged=true;
					break
				end
			end;

			iter.n=i;
			iter.tol=iter.tol(1:i);
			iter.rtol=iter.rtol(1:i);

			iter.time=toc(solutiontime);
			dpsolver.print(iter, ap);
		end % end of SOLVE.sa

		function [V, P, dV, iter]=policy(bellman, V0, ap)
			% dpsolver.nk: Solve for fixed point using Newton-Kantorovich iterations 
			%
			%  syntax:	[V, P, dV, iter]=dpsolver.nk(bellman, V0, ap):
			%
			%  INPUT:
			%     bellman:  [V, P, dV]=Bellman equation with fixed point, V
			%	 							Matlab function on the form [V,dV]=bellman(V)
			%								where V is the value function (m x 1), P is the policy function 
			%								and dV is the (m x m) Frechet derivative of the Bellman operator]
			%     V=  			Initial guess value function,V.
			%								[m x 1 matrix]
			%
			%  OUTPUT:
			%     V:         m x 1 matrix. Fixed point, V
			%			P:				 Policy function at fixed point
			%			dV:				 Frechet derivative of the Bellman operator


			% Set default settings for fixed point algorithm, for ap's not given in input
			% (overwrites defaults if ap is given as input)
			ap=dpsolver.setup(ap);

			solutiontime=tic;
			iter.tol=nan(ap.pi_max,1);
			iter.rtol=nan(ap.pi_max,1);
			iter.converged=false;

			m=numel(V0);			
		  	for i=1:ap.pi_max; %do at most pi_max N-K steps
	        	%Do policy iteration step
	        	[V1, P, dV, Eu]=bellman(V0); % also return value and policy function
	        	F=speye(m) - dV; % using dV from last call to bellman
	        	V=inv(F)*Eu; % policy-iteration


	        	% % do additional SA iteration for stability and accurate measure of error bound
	        	% V0=bellman(V); 
				% tolerance 
				iter.tol(i)=max(max(abs(V-V0)));
				V0=V;

		        %adjusting the N-K tolerance to the magnitude of ev
		        adj=ceil(log10(abs(max(max(V0)))));
		        ltol=ap.pi_tol*10^adj;  % Adjust final tolerance
		        ltol=ap.pi_tol;  % tolerance

		        if (iter.tol(i) < ltol);
		             % Convergence achieved
		             iter.message=sprintf('Policy converged after %d iterations, tolerance: %10g\n',i, iter.tol(i));
		             iter.converged=true;
		             break
		        end
			end %Next N-K iteration

			iter.time=toc(solutiontime);
			iter.n=i;
			iter.tol=iter.tol(1:i);
			iter.rtol=iter.rtol(1:i);

			dpsolver.print(iter, ap);
		end % end of SOLVE.nk

		function [V, P, dV, iter]=nk(bellman, V0, ap)
			% dpsolver.nk: Solve for fixed point using Newton-Kantorovich iterations 
			%
			%  syntax:	[V, P, dV, iter]=dpsolver.nk(bellman, V0, ap):
			%
			%  INPUT:
			%     bellman:  [V, P, dV]=Bellman equation with fixed point, V
			%	 							Matlab function on the form [V,dV]=bellman(V)
			%								where V is the value function (m x 1), P is the policy function 
			%								and dV is the (m x m) Frechet derivative of the Bellman operator]
			%     V=  			Initial guess value function,V.
			%								[m x 1 matrix]
			%
			%  OUTPUT:
			%     V:         m x 1 matrix. Fixed point, V
			%			P:				 Policy function at fixed point
			%			dV:				 Frechet derivative of the Bellman operator


			% Set default settings for fixed point algorithm, for ap's not given in input
			% (overwrites defaults if ap is given as input)
			ap=dpsolver.setup(ap);


			solutiontime=tic;
			iter.tol=nan(ap.pi_max,1);
			iter.rtol=nan(ap.pi_max,1);
			iter.converged=false;

			m=numel(V0);			
		  	for i=1:ap.pi_max; %do at most pi_max N-K steps
	        	%Do N-K step
	        	[V1, P, dV]=bellman(V0); % also return value and policy function
	        	if ap.check_deriv
	        		[df_nm, df_an, df_err] =dpsolver.check_deriv(bellman, V0, ap);
	        		if max(abs(df_err))>ap.check_deriv_tol;
	        			warning('Max abs difference between analytical and numerical derivatives of bellman operator is larger than tolerance')
	        		end
	        	end
	        	F=speye(m) - dV; % using dV from last call to bellman
	        	V=V0-F\(V0-V1); % NK-iteration

	        	% do additional SA iteration for stability and accurate measure of error bound
	        	V0=bellman(V); 
				
				% tolerance 
				iter.tol(i)=max(max(abs(V-V0)));

		        %adjusting the N-K tolerance to the magnitude of ev
		        adj=ceil(log10(abs(max(max(V0)))));
		        ltol=ap.pi_tol*10^adj;  % Adjust final tolerance
		        ltol=ap.pi_tol;  % tolerance

		        if (iter.tol(i) < ltol);
		             % Convergence achieved
		             iter.message=sprintf('N-K converged after %d iterations, tolerance: %10g\n',i, iter.tol(i));
		             iter.converged=true;
		             break
		        end
			end %Next N-K iteration

			iter.time=toc(solutiontime);
			iter.n=i;
			iter.tol=iter.tol(1:i);
			iter.rtol=iter.rtol(1:i);

			dpsolver.print(iter, ap);
		end % end of SOLVE.nk

		function iter=print(iter, ap)
			if (ap.printfxp>1); 	% print detailed output
				fprintf('iter           tol        tol(j)/tol(j-1) \n');
				for i=1:numel(iter.tol)
					fprintf(' %3.0f   %16.8e %16.8f\n',i, iter.tol(i),iter.rtol(i));
				end	
			end		

			if (ap.printfxp>0) 			% print final output
				if iter.converged==1;
					fprintf('%s\n',iter.message);
 				else
					fprintf('Maximum number of iterations reached, tolerance: %10g\n',iter.tol(end));
				end
				fprintf('Elapsed time: %3.5f (seconds)\n',iter.time)
			end
		end % end SOLVE.print


 	end % end of methods
end % end of SOLVE

