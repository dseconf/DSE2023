classdef mpec
  % MPEC class for structural estimation of discrete choice models.
  % By Bertel Schjerning & Patrick Kofod Mogensen
  methods (Static)
    function [results, theta_hat, Avar] = estim(data, mp)

      n_u = numel(struct2vec(mp,mp.pnames_u)); 
      n_P = numel(struct2vec(mp,mp.pnames_P));

      if n_P>0;
          error('MPEC only implemented for twostep=1');
      end

      samplesize=numel(data.d);

      % ************************************
      % STEP 1: ESTIMATE p 
      % ************************************

      % Estimate model using Partial MLE
      tab=tabulate(data.dx1);
      tab=tab(tab(:,3)>0,:);
      p=tab(1:end-1,3)/100;

      % ************************************
      % STEP 2: ESTIMATE structural parameters
      % ************************************

      [J_pattern, H_pattern] = mpec.sparsity_pattern(n_u,n_P,mp.n, n_P+1);

      % set options optimizer options for fmincon
      % options_mpec = optimset('DerivativeCheck','off','Display','on',...
      %     'GradConstr','on','GradObj','on','TolCon',1E-8,'TolFun',1E-8,'TolX',1E-15, ...
      %     'JacobPattern',J_pattern, 'HessPattern', H_pattern, ...
      %     'MaxFunEval', 100000, 'Algorithm','interior-point' ); 

    % set options optimizer options for fmincon
      options_mpec = optimset('DerivativeCheck','off','Display','off',...
          'GradConstr','on','GradObj','on','TolCon',1E-10,'TolFun',1E-8,'TolX',1E-15, ...
          'JacobPattern',J_pattern, 'HessPattern', H_pattern, ...
          'MaxFunEval', 100000, 'Algorithm','interior-point' ); 


      % options_mpec = optimset('DerivativeCheck','off','Display','iter',...
      %     'GradConstr','on','GradObj','on','TolCon',1E-10,'TolFun',1E-8,'TolX',1E-15, ...
      %     'MaxFunEval', 100000); 


      % starting values 
      theta_start_mpec = [struct2vec(mp,[mp.pnames_u mp.pnames_P]); zeros(mp.n,1)];

      % Upper and lower bounds on parameters
      lb = zeros(2+mp.n,1); ub = zeros(2+mp.n,1);
      lb(1) = -0.000000000001; ub(1) = inf; %thetaCost
      lb(2) = -0.000000000001; ub(2) = inf; %RC

      %Put bound on EV; this should not bind, but is a cautionary step to help keep algorithm within bounds
      lb(3:end) = -inf; ub(3:end) = 0; %EV
      lb(3:end) = -10000; ub(3:end) = 0; %EV

      % No linear equality constraints
      Aeq = []; beq = [];

      % Create mp.n dummy variables for observed mileage being equal to a particular value in the grid. 
      % Needed for derivatives of likelihood
      for i = 1:mp.n;
          data.xd(:,i) = (data.x == i);           
      end
      data.xd = sparse(data.xd);  % Utilize sparsity pattern. Matrix is N*T by mp.n and has only N*T nonzero elements. 

      % Define objective functions and constraints
      ll_p_mpec = @(theta) mpec.ll(data, mp, theta); % objective function
      con_p_bellman = @(theta) mpec.con_bellman(data, mp, theta); % Costraint (Bellman equation)     

      outsidetimer = tic;
      [theta_hat,FVAL,EXITFLAG,OUTPUT] = fmincon(ll_p_mpec,theta_start_mpec,[],[],Aeq,beq,lb,ub,con_p_bellman,options_mpec);  
      timetoestimate = toc(outsidetimer);

      % Compute Variance-Covaiance matrix
      Avar=nan(2,2);  
      
      % Collect results
      results=vec2struct(theta_hat,{'RC', 'c'},mp);
      results.cputime=timetoestimate;
      results.converged   =   (EXITFLAG == 1);

      results.MajorIter = OUTPUT.iterations;
      results.funcCount = OUTPUT.funcCount;

      results.llval   =   - ll_p_mpec(theta_hat)*samplesize;
    end

    function [f,g,h ]=ll(data, mp, theta) 
      global ev0;

      % parse parameters
      mp.RC=theta(1);
      mp.c=theta(2);
      if numel(theta)>(2+mp.n); % FIML: estimate state transition parameters
          n_p=numel(mp.p);
          mp.p=theta(3:3+n_p)';
          ev=theta(4+n_p:end);
      else % Two-step estimation. Hold state transition parameters fixed
          n_p=0;
          ev=theta(3:end);
      end

      % Update u, du and P evaluated in grid points
      dc=0.001*mp.grid;     % derivative of cost function 
      c=mp.c*0.001*mp.grid;

      VK=-c+mp.beta*ev;
      VR=-mp.RC-c(1)+mp.beta*ev(1);
      pk=1./(1+exp((VR-VK)));
      lp=pk(data.x);
      dk=(data.d==0);
      dr=(data.d==1);

      logl=log(dk.*lp+(1-lp).*dr);

      % add on log like for mileage process
      if n_p>0;
          logl=logl + log(mp.p(1+ data.dx1)');
      end

      f=-mean(logl);

      if nargout >=2; % compute gradient        
          res=lp-dk; % model "residual"
          g=zeros(2+mp.n,1);
          g(1)=-mean(res);                            %RC
          g(2)=mean(res.*(dc(data.x)-dc(1)));         %c
          g(3)=-mp.beta*mean(res.*(data.xd(:,1)-1));  %ev(1)

          NT=numel(res);
          for i=2:mp.n;        
              g(2 + n_p+ i)=-mp.beta*sum(res.*(data.x==i))/NT;       %ev(2:N)
              g(2 + n_p+ i)=-mp.beta*sum(res(data.x==i))/NT;         %ev(2:N)
              g(2 + n_p+ i)=-mp.beta*sum(res(data.xd(:,i)))/NT;      %ev(2:N) sparse index precalculated
          end
          g=-g;
      end

      if nargout >=3; % Compute hessian
          %Hessian is computed for only RC, c, and ev
          if n_p>0
              error 'RTFM!'
          end

          % Parameters are: {RC, c, ev}
          % Derivative of the diff between V(replace) and V(keep)
          ddiff=zeros(mp.n,2+mp.n);
          ddiff(:,1)=-1;                          % d(diff)/dRC
          ddiff(:,2)=dc-dc(1);                    % d(diff)/dc
          ddiff(2:end,3)=mp.beta*ones(mp.n-1,1);   % d(diff)/ev(1)
          ddiff(2:end,4:end)=-mp.beta*eye(mp.n-1);

          %Compute Hessian of likelihood [2x2]
          H=sparse(2+mp.n,2+mp.n);
          for i=1:2+mp.n
               for j=i:2+mp.n;
                  if ((i<=3)||(i==j))
                       % TODO: skip elements (sparsity), utilize symetry and precompute indices
                       H(i,j)=sum(pk(data.x).*(pk(data.x)-1).*ddiff(data.x,i).*ddiff(data.x,j));
                       H(j,i)=H(i,j);
                  end
              end
          end       
          h=-H/NT; % since we are minimizing the negative of the mean log-likleihood
      end
    end % end of ll

    function [c,ceq,DC,DCeq] = con_bellman(data, mp, theta)

      mp.RC=theta(1);
      mp.c=theta(2);
      if numel(theta)>(2+mp.n)
          error('RTFM! con_bellman not implemented for full mle')
      else
          n_p=0;
          ev=theta(3:end);
      end

      % Update u, du and P evaluated in grid points
      dc=       0.001*mp.grid;

      % Transition matrix for mileage
      P = zurcher.statetransition(mp); 
      u = zurcher.u(mp); 

      if nargout >= 3
        [ev1, pk, dev]=zurcher.bellman_ev(ev, mp, u, P);
      else
        [ev1, pk]=zurcher.bellman_ev(ev, mp, u, P);
      end

      % Define and evaluate nonlinear inequality constraints
      c = [];

      % Define and evaluate nonlinear equality constraints
      ceq = ev1-ev;

      % Define and evaluate the constraint Jacobian (DC, DCeq).   
      if nargout >= 3
          DC= [];
          DCeq=[];

          DCeq=sparse(mp.n,2+mp.n); 
          tmp=(P{1}*pk);
          DCeq(:,1)=- P{1}*(1*(1-pk));
          DCeq(:,2)=-P{1}*((dc-dc(1)).*pk);

          DCeq(:,3:end)=(dev-speye(mp.n));
          DCeq=DCeq';
      end
    end % end of con_bellman

    function [JacobSpaPattern, HessSpaPattern] = sparsity_pattern(nc,np,N, M)
      % nc: number of cost function parameters to be estimated
      % np: number of  free transition matrix parameters to be estimated 
      % N: number of EV parameters (i.e. number of gridpoints) 
      % M: number transition matrix parameters.

      BandDiag=0;
      for i=0:M-1
          a=diag(ones(N,1),i);
          BandDiag=BandDiag+a(1+i:end,1+i:end);
      end
      BandDiag(:,1)=ones(N,1);     
      JacobSpaPattern=[ones(N, nc+np) BandDiag];	

      HessSpaPattern = ones(nc+np+1,nc+np+N);
      HessSpaPattern = [HessSpaPattern; ones(N-1,nc+np+1) eye(N-1)];
    end % end of sparsity_pattern
  end
end

