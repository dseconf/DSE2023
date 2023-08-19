classdef zurcher
% zurcher class: Holds model parts for Rust's engine replacement model Rust(Ecta, 1987) 
% By Fedor Iskhakov, Bertel Schjerning, and John Rust
methods (Static)
  function mp = setup(mpopt)
    % zurcher.setup: Sets default parameters equal to parameters in Table X in Rust (1987).    
    mp.bustypes=[1,2,3,4];    % Vector with chosen bus-types in Rust's data (can be 1,2,3,4) 

    % Spaces
    mp.n=175;                       % Number of grid points     
    mp.max=450;                     % Max of mileage

    % Structural parameters
    mp.p     =[0.0937 0.4475 0.4459 0.0127]';  % Transition probabilities
    mp.RC    = 11.7257;                        % Replacement cost
    mp.c     = 2.45569;                        % Cost parameter
    mp.beta  = 0.9999;                         % Discount factor
    % Type of bellman equation: 
    %   mp.bellman_type='ev': use bellman equation in expected value function space 
    %                         (state and choice specific)
    %   mp.bellman_type='iv'  use bellman equation in integrated value function space 
    %                         (only state specific)
    mp.bellman_type='ev';  

    % parameters to estimate (set to empty to hold parameter fixed and skip estimation)
    mp.pnames_u={'RC', 'c'};  
    mp.pnames_P={'p'};  

    mp.ap=dpsolver.setup;

    if nargin>0
      pfields=fieldnames(mpopt);
      for i=1:numel(pfields);
        mp.(pfields{i})=mpopt.(pfields{i});
      end
    end

    global ev0 V0;  
    ev0=0; V0=0;

    mp.grid= (0:mp.n-1)';           % Grid over mileage
  end % end of zurcher.setup

  function [u, du] = u(mp)
    % zurcher.u: Computes state and choice specific utility and it's derivatives 
    %
    % Inputs
    %  mp       structure with model parameters (see setup)
    %
    % Outputs:
    %  u        mp.n x 2 matrix with state and choice specific pay-off:
    %           u(:,1): utility conditional on keeping engine
    %           u(:,2): utility conditional on replacing engine

    u=zeros(mp.n,2);
    u(:,1)=-0.001*mp.c*mp.grid;              % Pay-off - keeping
    u(:,2)=-mp.RC-0.001*mp.c*mp.grid(1);     % Pay-off - replacement

    if nargout>1
      n_u = numel(struct2vec(mp,mp.pnames_u)); 
      n_P = numel(struct2vec(mp,mp.pnames_P));
      du=zeros(mp.n, n_u+n_P, 2);   % note that derivatives of utility wrt p is zero
      for iP=1:numel(mp.pnames_u)        
        if strcmp(mp.pnames_u(iP), 'c'); 
          du(:,iP,1)= -0.001*mp.grid; 
          du(:,iP,2)= -0.001*mp.grid(1); 
        elseif strcmp(mp.pnames_u(iP), 'RC'); 
          du(:,iP,1)= 0; 
          du(:,iP,2)=-1; 
        end
      end
    end
  end % end of zurcher.u

  function P = statetransition(mp)
    % zurcher.statetransition: Computes state transition matrices conditional on choice. 
    %
    % Inputs
    %  mp       structure with model parameters (see setup)
    %
    % Outputs:
    %  P        2 dimensional cell array with mp.n x mp.n conditional transition matrices:
    %           P{1} is transition matrix conditional on keeping engine
    %           P{2} is Transition matrix conditional on replacing engine
    
    % First compute transition matrix conditional on keeping engine  
    p=[mp.p; (1-sum(mp.p))];
    n=mp.n;
    P{1}=0;
    for i=0:numel(p)-1;
      P{1}=P{1}+sparse(1:n-i,1+i:n,ones(1,n-i)*p(i+1), n,n);
      P{1}(n-i,n)=1-sum(p(1:i));
    end
    P{1}=sparse(P{1});

    % Then transition matrix conditional on replacing engine
    P{2}=sparse(n,n); 
    for i=1:numel(p);
      P{2}(:,i)=p(i);
    end
  end % end of zurcher.statetransition

  function [V1, pk, dbellman_dV]=bellman(V0, mp, u, P)

    % zurcher.bellman:   Procedure to compute Bellman operator
    %   mp.bellman_type='ev': use bellman equation in expected value function space 
    %                         (state and choice specific)
    %   mp.bellman_type='iv'  use bellman equation in integrated value function space 
    %                         (only state specific)
    % To obtain ev=bellamn_ev(ev) from V=bellamn_iv(V) compute ev=P{1}*V

    % Inputs and outputs: see bellman_iv and zurcher.bellman_ev

    if numel(V0)==1;
      V0=zeros(mp.n,1);
    end

    if nargin<3
      u = zurcher.u(mp);
    end

    if nargin<4
      P = zurcher.statetransition(mp);
    end

    if strcmp(mp.bellman_type, 'iv')
      [V1, pk, dbellman_dV]=zurcher.bellman_iv(V0, mp, u, P);
    elseif strcmp(mp.bellman_type, 'ev')
      [V1, pk, dbellman_dV]=zurcher.bellman_ev(V0, mp, u, P);
    else
      error('mp.bellman_type must be iv or ev')
    end
  end % end of zurcher.bellman

  function [V1, pk, dBellman_dV]=bellman_iv(V0, mp, u, P)
    % zurcher.bellman_iv:  integrated Bellman equation
    %
    % Inputs
    %  V0       mp.n x 1 matrix of initial guess on integrated value function
    %  mp       structure with model parameters (see setup)
    %  u        mp.n x 2 matrix of state and choice specific utilities
    %  P        2 dimensional cell array with mp.n x mp.n conditional 
    %           transition matrices for each choice (see statetransition)
    %
    % Outputs:
    %  V1       mp.n x 1 matrix of integrated value function given initial guess V0
    %  pk       mp.n x 1 matrix of choice probabilities (Probability of keep)

    vK= u(:,1) + mp.beta*P{1}*V0;   % Value of keeping
    vR= u(:,2) + mp.beta*P{2}*V0;   % Value of replacing

    % Need to recenter logsum by subtracting max(vK, vR)
    maxV=max(vK, vR);
    V1=(maxV + log(exp(vK-maxV)  +  exp(vR-maxV))); 

    % If requested, also compute choice probability from ev (initial input)
    if nargout>1 
      pk=1./(1+exp((vR-vK)));
    end

    if nargout>2 % compute Frechet derivative
      dBellman_dV=mp.beta*(P{1}.*pk + P{2}.*(1-pk));
    end
  end % end of zurcher.bellman_iv  

  function [ev, pk, dbellman_dev]=bellman_ev(ev0, mp, u, P)
    % zurcher.bellman_EV: bellman equation based on expected values
    %                     (Similar to Rust NFXP_Manual)
    % Inputs
    %  ev0      mp.n x 1 matrix of expected values given initial guess on value function
    %  u        mp.n x 2 matrix of state and choice specific utilities
    %  mp       structure with model parameters (see setup)
    %  P        2 dimensional cell array with mp.n x mp.n conditional 
    %           transition matrices for each choice (see statetransition)
    %
    % Outputs:
    %  ev       mp.n x 1 matrix of expected values given initial guess of ev 
    %  pk       mp.n x 1 matrix of choice probabilities (Probability of keep)

    vK= u(:,1) + mp.beta*ev0;      % Value off keep
    vR= u(:,2) + mp.beta*ev0(1);   % Value of replacing

    % Need to recenter logsum by subtracting max(vK, vR)
    maxV=max(vK, vR);
    V=(maxV + log(exp(vK-maxV)  +  exp(vR-maxV))); 
    ev=P{1}*V;


    % If requested, also compute choice probability
    if nargout>1 
      pk=1./(1+exp((vR-vK)));
    end
    if nargout>2 % compute Frechet derivative
      dbellman_dev=mp.beta*(P{1}.*pk');    
      % Add additional term for derivative wrt Ev(1), 
      % since Ev(1) enter logsum for all states            
      dbellman_dev(:,1)=dbellman_dev(:,1)+mp.beta*P{1}*(1-pk);        
    end
  end % end of zurcher.bellman_ev  
    
  function [f,g,h]=ll(data, mp, theta)
    global V0;
    
    % Update model parameters
    mp=vec2struct(theta, [mp.pnames_u mp.pnames_P], mp);
    mp.p=abs(mp.p); % HACK: force transition probabilities to be positive

    % Update transition matrices 
    P = zurcher.statetransition(mp); 

    % Update pay-off 
    if nargout >=2;  % utility and it's compute derivatives          
      [u, du] = zurcher.u(mp);   
    else % skip derivatives
      u = zurcher.u(mp);
    end         

    % Solve model
    bellman= @(V) zurcher.bellman(V, mp, u, P);
    [V0, pk, dBellman_dV]=dpsolver.poly(bellman, V0, mp.ap, mp.beta);  
    
    % Evaluate likelihood function
    % 1. log likelihood regarding replacement choice
    y_j=[(1-data.d) data.d];         % choice dummies [keep replace]
    px_j=[pk(data.x) 1-pk(data.x)];  % ccps at observed states [keep replace]
    logl=log(sum(y_j.*px_j,2));      % log likelihood

    % 2. add on log likelihood for mileage process
    if numel(mp.pnames_P)>0;
      p=[mp.p; 1-sum(mp.p)];
      logl=logl + log(p(1+ data.dx1));
    end

    % Objective function (negative mean log likelihood)
    f=mean(-logl);

    if nargout >=2;       % compute scores and gradient of likelihood
      score = zurcher.score(data, mp, P, pk, px_j, V0, du, dBellman_dV);
      g=mean(-score,1);   % gradient (of negative of mean log likelihood)
    end

    if nargout >=3;  % compute BHHH Hessian approximation
      h=score'*score/(size(logl, 1)); 
    end 
  end % end of zurcher.ll_integrated 

  function score = score(data, mp, P, pk, px_j, V0, du, dBellman_dV);
      % Compute scores (use chain rule - three steps)

      y_j=[(1-data.d) data.d];                  % choice dummies [keep replace]

      n_u = numel(struct2vec(mp,mp.pnames_u));  % number of utility parameters
      n_P = numel(struct2vec(mp,mp.pnames_P));  % number of transition parameters

      % STEP 1: compute derivative of bellman operator wrt. parameters
      % ...first wrt. utility parameters
      dbellman=pk.*du(:,:,1) + (1-pk).*du(:,:,2); 

      if strcmp(mp.bellman_type, 'ev');
        dbellman=P{1}*dbellman;
      end

      % ...then wrt. state transition parameters
      if n_P>0        
        for iP=1:n_P;
          dbellman(1:mp.n-iP, n_u+iP)= [V0(iP:end-1)] - [V0(n_P+1:mp.n); repmat(V0(end), n_P-iP, 1)];
        end
      end

      % STEP 2: compute derivative of fixed point, V, wrt. all parameters
      dV=(speye(size(dBellman_dV)) - dBellman_dV)\dbellman;  

      % STEP 3: compute derivative of log-likelihood wrt. parameters
      % ...first wrt. utility parameters
      score=0;
      for j=1:size(y_j, 2); 
        dv= du(:,:,j) + mp.beta*P{j}*dV; 
        score = score+ (y_j(:,j)-px_j(:,j)).*dv(data.x,:); 
      end

       % ...then wrt. state transition parameters
      if n_P>0
        p=[mp.p; 1-sum(mp.p)];
        invp=exp(-log(p));
        invp=[sparse(1:n_P,1:n_P,invp(1:n_P),n_P,n_P); -ones(1,n_P)*invp(n_P+1)];
        for iP=1:n_P;
          score(:,n_u+iP)= score(:,n_u+iP) + invp(1+data.dx1,iP);
        end
      end
  end % end of zurcher.score

  function [data] = simdata(N,T,mp, P, pk, seed)
    % zurcher.simdata: simulates data from engine replacement model. 
    %
    % Inputs
    %  N       Number busses to simulate
    %  T       Number of time periods to be simulated for each bus
    %  mp       structure with model parameters (see setup)
    %  P        2 dimensional cell array with mp.n x mp.n conditional transition matrices (see statetransition)
    %  pk       mp.n x 1 matrix of choice probabilities (probability of keep)

    % Initial conditions: all buses start uniformly distributed across statespace

    % Outputs:
    %  data: table with NxT rows and 6 columns 
    %       data.id  :  Bus id 
    %       data.t   :  Time period
    %       data.d   :  Replacement indicator (Replace = 1, Keep=0)
    %       data.x   :  Mileage in period t
    %       data.x1  :  Mileage in period t+1
    %       data.dx1 :  Change in mileage between period t and t+1

    if exist('seed')
      rng(seed,'twister');
    end

    id=repmat((1:N),T,1);
    t=repmat((1:T)',1,N);
    u_init=randi(mp.n,1,N);
    u_dx=rand(T,N);
    u_d=rand(T,N);

    csum_p=cumsum(mp.p);
    dx1=0;
    for i=1:numel(csum_p);
      dx1=dx1+ (u_dx>(csum_p(i)));
    end;

    x=nan(T, N);
    x1=nan(T, N);
    x(1,:)=u_init; % Initial conditions

    for it=1:T;
      d(it,:)=(u_d(it,:)<(1-pk(x(it,:)')'))*1;  % Replace = 1, Keep=0
      x1(it,:)=min(x(it,:).*(1-d(it,:)) + d(it,:) + dx1(it,:), mp.n);
      if it<T;
        x(it+1,:)=x1(it,:);
      end
    end

    data.id=id;
    data.t=t;
    data.d=d;
    data.x=x;
    data.x1=x1;
    data.dx1=dx1;

    pfields=fieldnames(data);
    for i=1:numel(pfields);
      data.(pfields{i})=reshape(data.(pfields{i}), T*N,1);
    end
    data=struct2table(data);
  end % end of zurcher.simdata

  function dta = readbusdata(mp)
    load('busdata1234.mat');

    % Select busses
    data=data(ismember(data(:,2), mp.bustypes),:);
    
    id=data(:,1);       % Bus id
    bustype=data(:,2);  % bustype: 1,2,3,4
    d1=data(:,5);       % Lagged replacement dummy, d_(t-1)
    d=[d1(2:end);0];    % Replacement dummy, d_t
    x=data(:,7);        % Odometer, x_t


    % Discretize odometer data into 1, 2, ..., n
    x=ceil(x.*mp.n/(mp.max*1000));                                 

    % Monthly mileage
    dx1=x(1:end,1)-[0;x(1:end-1,1)];           % Take first difference on odometer                           
    dx1=dx1.*(1-d1)+x.*d1;                     % Make true x_t-x_(t-1) (replace first difference by x_t if replacement dummy lagged is 1)

    % remove observations with missing lagged mileage
    data=[id, bustype, d,x,dx1];                                           % [ id , replace dummy lagged, odometer, "change in odometer" , bustype]
    remove_first_row_index=data(:,1)-[0;data(1:end-1,1)];                  % Take first difference of ID...
    data=data((remove_first_row_index==0),:);                              % ... and only keep lines where ID hasn't changed 

    % Save data structure
    dta.d=data(:,3);           
    dta.x=data(:,4);
    dta.dx1=data(:,5);
    dta=struct2table(dta);
  end % end of zurcher.readbusdata    

  function [pp, pp_K, pp_R] = eqb(mp,P, pk)
    % This program computes the equilibrium distribution of the controlled
    % stochastic process P(i(t)|x(t))p(x(t)|x(t-1),i(t-1)). The equilibrium
    % distribution pp is an (1,n) matrix, where n is the number of discrete
    % mileage cells: pp=pp*P, where P is the transition probability matrix
    % for the (nxn) merged state process, where states with i=1 are 
    % identified with the state x=1. Estimates of the complete (1,2*n)
    % equilibrium distribution are then uncovered from pp. 

    % Outputs    
    % pp: Pr{x} (Equilibrium distribution of mileage)
    % pp_K: Pr{x,i=Keep}
    % pp_R: Pr{x,i=Replace}

    function ed=ergodic(p);
      % ergodic.m: finds the invariant distribution for 
      % an NxN Markov transition probability
      n=size(p,1);
      if n ~= size(p,2);
        fprint('Error: p must be a square matrix\n');
        ed=NaN;
      else
        ap=eye(n)-p';
        ap=[[ap; ones(1,n)], ones(n+1,1)];
        if (rank(ap) < n+1);
          fprintf('Error: transition matrix p is not ergodic\n');
          ed=NaN;
        else
          ed=[ones(n,1); 2];
          ed=ap\ed;
          ed(end)=[];
          ed=reshape(ed,1,[]); %row vector
        end;
      end
    end % end of ergodic

    % state transition matrix 
    pl = P{1}.*pk +P{2}.*(1-pk);

    % Equilibrium distribution of mileage
    pp = ergodic(pl);   % Pr{x} (Equilibrium distribution of mileage)
    pp_K=pp.*(pk');     % Pr{x,i=Keep}
    pp_R=pp.*(1-pk');   % Pr{x,i=Replace}
  end % end of zurcher.eqb

end % end of methods
end % end of zurcher class