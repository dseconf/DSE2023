% This class is generic solver for non-linear equations
% using newtons method to solve equations on the form F(x)=0 
% 
% [ap]=setup(apopt);
% [x, f, convergence]=newton.solve(F, x0, ap);
% [df_nm, df_an, df_err] =check_deriv(F, x0);
% 
% example [x, f, convergence]=newton.solve(@(x) x.^2 + x -10, 1)
%
% syntax of main function: 
% 	[x, f, convergence]=newton.solve(fun, x0, ap)
%
% July 2020

classdef newton
	methods (Static)

		function [x, f, convergence]=solve(F, x0, ap)
			% newton.solve:  solves newton.solve attempts to solve non-linear equations of the form:
			%          
			% F(x) = 0    where F and x may be vectors   
			%
			%
			%  syntax:	[x, f, convergence]=newton.solve(fun, x0, ap)
			%
			%  INPUT:
			%     F:  			Matlab function on the form [F,dF]=F(x)
			%								where x nx1 vector of variable to solve for and dF is the (m x n) matrix of frist order derivatives]
			%     x0=  			Initial guess, x.
			%								[n x 1 matrix]
			%			ap:				Algorithm paramaters. See newton.setup
			%
			%			bet:		  Discout factor. Enters rule for stopping SA and switching to NK iterations. 
			%								SA should stop prematurely when relative tolerance is close to bet. 
			%
			%  OUTPUT:
			%     x:         n x 1 matrix. Fixed point, V
			%	
			%			F:				 Function at solution
			%
			%     convergence: 
			%					1 of newton.solve converged to a root.
			%					0 of newton.solve did not converge
			%-----------------------------------------------------------------------------------------------------------------------------

			% Set default settings for newton.solve for ap's not given in input
			% (overwrites defaults if ap is given as input)
			if nargin>=3;
				ap=newton.setup(ap);
			else
				ap=newton.setup;
			end

			if any(strcmp(ap.display, {'iter', 'on'}));
			  fprintf('Starting newton.solve\n');
			end

			if ap.check_deriv==1
				if ap.use_deriv==1
					[df_err, df_nm, df_an]=newton.check_deriv(F, x0);
					if any(norm(df_err)>ap.deriv_tol)
						warning('Analytical derivatives do not match at initial value, x0')
					else
						fprint('Analytical and numerical derivatives match at initial value, x0\n')
					end 
				else 
					warning('Must have ap.use_deriv=1 to test derivatives')
				end
			end 

			solutiontime=tic;

	   	convergence=0;
	    tol=1;

	    i=1;
	    while (tol > ap.ftol & i < ap.maxiter)
	    	if 	ap.use_deriv;	% use analytical derivatives
	      	[Fx, dFx]=F(x0);
	     	else % use numerical derivatives
	     		[Fx]=F(x0);
	     		dFx=gradp(F, x0);
	     	end
	      tol=max(abs(Fx));
	      x1=x0-inv(dFx')*Fx;
	      % pp=p-inv(dpp)*pp;

	       if tol < ap.ftol
	       	convergence=1;
	       end
	       if any(strcmp(ap.display, {'iter'}));
	         fprintf('i=%i tol=%g\n',i,tol);
	       end
	       x0=x1;
	       i=i+1;
	    end
	    x=x0; 

			if any(strcmp(ap.display, {'iter', 'on'}));
				if convergence
					fprintf('Convergence achieved!\n');
				else
					fprintf('No convergence! Maximum number of iterations exceeded without convergence!');
				end
				fprintf('Elapsed time: %3.5f (seconds)\n',toc(solutiontime));
			end
	 	end % end of newton.solve 

  	function ap=setup(apopt)
  		% newton.setup: setup of algorithm parameters, ap, used in newton.solve
			%  syntax: ap=setup(apopt)
			%
			%  INPUT:  apopt (optional): If apopt is specified, default parameters will overwritten with elements in apopt.
			%
			%  OUTPUT: apopt: algorithm parameter structure
			%
			% See also:
			%   newton.solve

			% default values of ap
			ap.use_deriv=	1;		  	% use analytical derivatives
			ap.check_deriv=	0;			% Check analytical derivatives of F against numerical derivatives 
			ap.deriv_tol=1e-5;			% Tolerance on norm of error of derivatives; 
			ap.maxiter=		200;		% Maximum number of iterations
			ap.ftol=		1.0e-10;	% Tolerance on F(x)
 			ap.display=		'iter'; 	% Display output from solution algorithm
	    								% 'off' (No printing), 'on' (Summary info), 'iter' (itetation info)
			ap.check_deriv_tol=	1e-6;	% Absolute tolerance before for difference between analytical and numerical derivatives of bellman operator

    	if nargin>0
      	pfields=fieldnames(apopt);
      	for i=1:numel(pfields);
          	ap.(pfields{i})=apopt.(pfields{i});
      	end
    	end
		end % end of newton.setup

		function [df_err, df_nm, df_an] =check_deriv(F, x0)
			[f_nm, df_an]=F(x0); 
			df_nm=gradp(F,x0); 
			df_err=df_an-df_nm;
		end % end of newton.check_deriv

		function [F, dF] =testfun(x)
			F=x.^2 + x - 100;
			dF=2*x + 1;
		end % end of newton.testfun



 	end % end of methods
end % end of newton