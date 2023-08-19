% NPL class: Solves Rust's engine replacement model Rust(Ecta, 1987) 

classdef npl
    methods (Static) 
        function [pk]=Psi(mp, pk0, P)  % Policy iteration operator, pk=Psi(pk0)             
            % Psi is the mapping taking choice probabilities pk0 given current model parameters
            % and Markov transitions, and returning the optimal ccps pk given pk0 and parameters. 
            Vsigma=npl.phi(mp, pk0, P);
            pk=npl.lambda(Vsigma, mp, P);
        end % end of Psi

        function [Vsigma]=phi(mp, pk, P, Finv)  % Policy evaluation operator, Vsigma=phi(pk)
            % This function does most of the work for calculating phi            
            if nargin <=3
                Fu = bsxfun(@times, P{1}, pk)  + bsxfun(@times, P{2}(1,:), 1-pk);
                Finv=inv(speye(mp.n) - mp.beta*Fu);
            end
            eulerc=0.5772156649015328606065120900824024310421;

            c=mp.c*0.001*mp.grid;
            vK= -c          + eulerc - log(pk);
            vR= -c(1)-mp.RC + eulerc - log(1-pk);

            pv = bsxfun(@times, vK, pk) + bsxfun(@times, vR, 1-pk);
            Vsigma= Finv*pv;
        end % end of phi

        function [pK]=lambda(Vsigma, mp, P) % Calculates the choice probabilities at the grid points.          
            if numel(Vsigma)==1
                v=zeros(size(mp.grid,1),3);
            end
            c = mp.c*0.001*mp.grid;
            VK = - c + mp.beta*P{1}*Vsigma;
            VR = - mp.RC - c(1) + mp.beta*P{2}*Vsigma;
            pK = 1./(1+exp(VR-VK));   
        end
        
        function [f,g,h ]=ll(theta, pk0, data, P, mp, Finv) % pseudo log-likelihood
            %% Procedure for calculating the likelihood value, gradients and Hessian approximation
            %  based on AM (2002) and Rust (1987).
            global pk
            % parse parameters
            mp.RC=theta(1);  mp.c=theta(2);

            % The following two lines is equivalent to calling Psi(mp, pk0, P) except that 
            % we utilize that we can hold Finv fixed for a given K, when computing Psi (for speed)
            [Vsigma]=npl.phi(mp, pk0, P, Finv);         % Smoothed value function
            pk=npl.lambda(Vsigma, mp, P);               % Probability of keep


            pKdata=pk(data.x);
            dk=(data.d==0);
            dr=(data.d==1);

            logl=log(pKdata.*dk+(1-pKdata).*dr);
            
            f=-mean(logl);

            NT=size(logl, 1);
        
            if nargout >=2; % compute gradient      
                dc=0.001*mp.grid;

                res=dk-pKdata;                   % model "residual"
                score=zeros(NT,numel(theta));
                dP=bsxfun(@minus, P{2}(1,:), P{1});

                dvdRC=(-1 + mp.beta*dP*Finv*((1-pk).*(-1)));
                dvdc=(dc + mp.beta*dP*Finv*(pk.*(-dc)));
                score(:,1)=res.*dvdRC(data.x); % RC
                score(:,2)=res.*dvdc(data.x); % c

                g=mean(score);   
            end
            if nargout >=3; % Compute BHHH approx to Hessian
                h=score'*score/NT;
            end
        end % end of ll

        function [f,g,h]=ll_logit(theta, y, x) % log-likelihood logit
            % px=exp(x*theta)./(1+exp(x*theta));
            px=1./(1+exp(-x*theta));
            logl=y.*log(px)+(1-y).*log(1-px);
            f=-mean(logl);
        end % end of ll_logit

        function [mp, pk, logl, K]=estim(theta0, pk0, data, P, mp, Kmax);  % NPL algorithm for K-PI estimators
            global pk

            options =  optimset('Algorithm','trust-region','Display','off', 'GradObj','on', 'TolFun',1E-10,'TolX',0,'Hessian','on');
            if nargin==5;
                Kmax=100;
            end

            %% Print estimation results
            fprintf('*************************************************************************\n'); 
            fprintf('Fixed parameters and spaces\n'); 
            fprintf('*************************************************************************\n'); 
            fprintf('Beta           = %10.5f \n',mp.beta);
            fprintf('n              = %10.5f \n',mp.n);
            fprintf('Sample size    = %10.5f \n',size(data.x,1));
            fprintf('\n'); 
            fprintf('*************************************************************************\n'); 
            fprintf('Method: Nested Pseudo Likelihood (NPL)\n'); 
            fprintf('*************************************************************************\n'); 
            fprintf('\n\n');

            fprintf('Parameter Estimates\n');
            fprintf('--------------------------------------------------------------------------\n');
            fprintf('%10s     %14s %14s %14s \n', '', 'RC', 'C' ,'log-like');
            fprintf('--------------------------------------------------------------------------\n');

            tic
            for K = 1:Kmax
                % Step 0)  Pre-compute unconditional transition matrix Fu and the inverse of I-beta*Fu to used in step 1 

                Fu = bsxfun(@times, P{1}, pk0)  + bsxfun(@times, P{2}, 1-pk0);
                Finv=inv(speye(mp.n) - mp.beta*Fu);

                % Step 1)  Maximize the pseudo-likelihood function given step K-1 CCPs
                [theta_npl,FVAL,EXITFLAG,OUTPUT1] = fminunc(@(theta) npl.ll(theta, pk0, data, P, mp, Finv),theta0,options);

                logl=-FVAL;
                NPL_metric = max(abs(theta0-theta_npl));

                % Step 2)  Update CCPs using theta_npl from step 1)
                pk0=pk;              % CCPs step K-1

                theta0 = theta_npl;  % Theta step K-1 
               
                %% Print estimation results
                fprintf('K = %-10d %14.3f %14.3f %14.3f \n' , ...
                K, mean(theta_npl(1)) ,  mean(theta_npl(2)) , mean(logl)*size(data.x,1)); 

                if NPL_metric < 1e-4
                    break
                end
            end
            timenpl=toc;
            fprintf('--------------------------------------------------------------------------\n');
            fprintf('Time to solve engine replacement model: %6.6f seconds\n', timenpl);

            fprintf('\n\n');
            mp.RC=theta_npl(1); mp.c=theta_npl(2);
        end % end of estim

        function [pk_solve] = solve(mp, P, pk0, fig)  % Solves pk=Psi(pk) using successive approximations 
            % fig is optional. If fig is input, npl will plot replacement probability for each iteration      
            tol = 1;
            for i_Psi = 1:100
                pk1 = npl.Psi(mp, pk0, P);

                if nargin>3
                    figure(fig)
                    hold all
                    plot(mp.grid,1-pk1,'-r','LineWidth', 2);
                end
                
                tol = max(abs(pk1-pk0));
                tol_vec(i_Psi) = tol;
                pk0 = pk1;
                if tol < 1e-12
                    break
                end
            end
            pk_solve = pk0;
            if nargin>3
                plot(mp.grid,1-pk_solve,'-b','LineWidth', 2);
                title('Convergence of \Psi');
                xlabel('Milage grid')
                ylabel('Replacement probability')
            end
        end
    end % end of methods
end % end of estim class