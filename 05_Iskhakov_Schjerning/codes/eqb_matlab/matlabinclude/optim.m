% This class iinlcudes several maximization routines
% 
% [ap]=setup(apopt);
% [x, f, convergence]=optim.solve(F, x0, ap);
% [df_nm, df_an, df_err] =check_deriv(F, x0);
% 
% example [x, f, convergence]=optim.solve(@(x) x.^2 + x -10, 1)
%
% syntax of main function: 
%   [x, f, convergence]=optim.solve(fun, x0, ap)
%
% July 2020

classdef optim
  methods (Static)

    function ap=setup(apopt)
      % optim.setup: setup of algorithm parameters, ap, used in optim.solve
      %  syntax: ap=setup(apopt)
      %
      %  INPUT:  apopt (optional): If apopt is specified, default parameters will overwritten with elements in apopt.
      %
      %  OUTPUT: apopt: algorithm parameter structure
      %
      % See also:
      %   optim.solve

      % default values of ap
      ap.hess_update='bfgs'; % 'fd', 'user' 
      ap.use_deriv= 1;        % use analytical derivatives
      ap.check_deriv= 0;      % Check analytical derivatives of F against numerical derivatives 
      ap.deriv_tol=1e-2;      % Tolerance on norm of error of derivatives; 
      ap.maxiter=   200;    % Maximum number of iterations
      ap.xtol=    1.0e-8;  % Tolerance on x1-x0
      ap.ftol=    1.0e-10;  % Tolerance on F(x)
      ap.gtol=    1.0e-6;   % Tolerance on gradient
      ap.display=   'iter';   % Display output from solution algorithm
                      % 'off' (No printing), 'on' (Summary info), 'iter' (itetation info)
      ap.check_deriv_tol= 1e-6; % Absolute tolerance before for difference between analytical and numerical derivatives of bellman operator

      if nargin>0
        pfields=fieldnames(apopt);
        for i=1:numel(pfields);
            ap.(pfields{i})=apopt.(pfields{i});
        end
      end
    end % end of optim.setup

    function [x, fx, convergence]=solve(fun, x0, ap)
      % optim.solve:     solves optim.solve attempts to maximize or minimize non-linear equations
      %          
      % max fun(x) = 0   where F is scalar and x may be a vector   
      %
      %
      % syntax: [x, f, convergence]=optim.solve(fun, x0, ap)
      %
      %  INPUT:
      %    fun:   Matlab function on the form [f,g,h]=fun(x)
      %           where x nx1 vector of variable to solve for and dF is the (m x n) matrix of frist order derivatives]
      %
      %     x0:   Initial guess, x (n x 1 matrix)
      %          
      %     ap:   Algorithm paramaters. See optim.setup
      %
      %
      %  OUTPUT:
      %     x:    local minimizer/maximizer (n x 1 matrix)
      % 
      %     f:    Function at optimum
      %
      %     convergence: 
      %         1 of optim.solve converged to a local optimum
      %         0 of optim.solve did not converge
      %-----------------------------------------------------------------------------------------------------------------------------

      % Set default settings for optim.solve for ap's not given in input
      % (overwrites defaults if ap is given as input)
      if nargin>=3;
        ap=optim.setup(ap);
      else
        ap=optim.setup;
      end
      if any(strcmp(ap.display, {'iter', 'on'}));
        fprintf('Starting optim.solve using %s as hessian udpate\n', ap.hess_update);
      end

      % check derivatives at x0
      if ap.check_deriv==1 & ap.use_deriv==1
        optim.check_deriv(fun, x0, ap);
      end 

      solutiontime=tic;
      convergence=0;
      iter=0;
      
      % STEP 0: initialize f, g and H at x0
      if  ap.use_deriv; % use analytical derivatives
        [fx0, gx0]=fun(x0);
      else % use numerical derivatives
        fx0=fun(x0);
        gx0=gradp(fun, x0);
      end

      if strcmp(ap.hess_update, 'fd')
        if  ap.use_deriv; % use analytical derivatives
          Hx0=gradp(@(x) optim.gradfun(fun, x), x0); 
        else
          Hx0=hessp(fun, x0);
        end
      else
        Hx0=eye(numel(x0)); 
      end

      while (iter < ap.maxiter & convergence==0)

        % STEP 1: compute search direction
        s0=-Hx0\gx0';

        % STEP 2: find step length
        fx1=fun(x0+s0);
        if fx1<=fx0; % accept lambda=1 and move on without line search
          lambda=1;
        else % Failed to decrease objective: line searching....
          opt.Display ='off';
          opt.MaxIter =100;
          lambda=fminbnd(@(lambda) fun(x0+lambda*s0), -1, 4, opt); 
        end
        
        % STEP 3: update parameter vector
        x1=x0+lambda*s0;
        dx=x1-x0;

        % print iteration info
        if any(strcmp(ap.display, {'iter'}));
           optim.print_out(iter, x0, x1, fx0, gx0, lambda)
        end

        % STEP 4: update objective function, gradient and Hessian at x1
        if  ap.use_deriv; % use analytical derivatives
          [fx1, gx1]=fun(x1);
        else % use numerical derivatives
          fx1=fun(x1);
          gx1=gradp(fun, x1);
        end

        if strcmp(ap.hess_update, 'fd')
          if  ap.use_deriv; % use analytical derivatives
            Hx1=gradp(@(x) optim.gradfun(fun, x), x1); 
          else
            Hx1=hessp(fun, x1);
          end
        elseif strcmp(ap.hess_update, 'bfgs')
          y=gx1'-gx0';
          Hx1=Hx0- Hx0*dx*dx'*Hx0/(dx'*Hx0*dx) + y*y'/(y'*dx);
        end

        % STEP 5: stopping rule
        if  max(abs(gx0)) < ap.gtol*(1+abs(fx1));
            convergence=1;
        end

        if max(abs(dx)) < ap.xtol*(1+max(abs(x1)));
          if  max(abs(gx0)) < ap.gtol*(1+abs(fx1));
            convergence=1;
          else
            convergence=-1;
          end
        end
        
        % prepare for next iteration
        x0=x1;
        fx0=fx1;
        gx0=gx1;
        Hx0=Hx1;
        iter=iter+1;
      end
      x=x0; 
      fx=fx0; 
      gx=gx0; 

      if any(strcmp(ap.display, {'iter', 'on'}));
        if convergence==1
          fprintf('Convergence to optimal point achieved\n');
        elseif convergence==-1
          fprintf('Convergence to non-optimal point\n');
        else
          fprintf('Maximum number of iterations exceeded without convergence!\n');
        end
        fprintf('Elapsed time: %3.5f (seconds)\n',toc(solutiontime));
      end
    end % end of optim.solve 

    function []=print_out(iter, x0, x1, fx, gx, step);
      if iter==0
          fprintf('%4s %10s  %13s  %13s  %13s\n', 'iter','step-length', '||x1-x0||','||g(x)||', 'f(x)')
          fprintf('%s\n', repmat('-',100,1))
      end
      fprintf('%4d  %10.2g  %13.6g  %13.2e  %13.6g\n', iter,step, max(abs(x1-x0)),max(abs(gx)),mean(fx))
    end

    function [df_err, df_nm, df_an] =check_deriv(fun, x0, ap)
      [f_nm, df_an]=fun(x0); 
      df_nm=gradp(fun,x0); 
      df_err=df_an-df_nm;
      if any(norm(df_err)>ap.deriv_tol)
        warning('Analytical derivatives do not match at initial value, x0')
      else
        fprintf('Analytical and numerical derivatives match at initial value, x0\n')
      end 
    end % end of optim.check_deriv

    function [g] =gradfun(fun, x0)
      [~, g]=fun(x0);
      g=reshape(g, numel(g), 1);
    end % end of optim.testfun

    function [f, g] =rosenbrock(x,a,b) 
      f=(a-x(1))^2+b*(x(2)-x(1)^2)^2;
      g=nan(1,2);
      g(1)=-2*(a-x(1))-4*x(1)*b*(x(2)-x(1)^2);
      g(2)=2*b*(x(2)-x(1)^2);
    end

  end % end of methods
end % end of newton