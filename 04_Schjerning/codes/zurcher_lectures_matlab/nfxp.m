classdef nfxp
  % NFXP class: Estimates rust's engine repplacement model Rust(Ecta, 1987) 
  % By Fedor Iskhakov, Bertel Schjerning, and John Rust
  methods (Static)
    function [results, theta_hat, Avar] = estim(data, mp)
      samplesize=numel(data.d);
 
      % set options optimizer options for fminunc      
      optim_options= optimset('Algorithm','trust-region', 'Display','off', 'GradObj','on', 'TolFun',1E-8,'TolX',1E-8 ,'Hessian','on');

      % **********************************************************************************
      % STEP 1: ESTIMATE state transition matrix using frequency estimator
      % **********************************************************************************
      tab=tabulate(data.dx1);
      tab=tab(tab(:,3)>0,:);

      % Use first step estimates as starting values for p
      mp.p=tab(1:end-1,3)/100;

      % **********************************************************************************
      % STEP 2: ESTIMATE structural parameters
      % **********************************************************************************

      outsidetimer=tic;
      
      % STEP 2a: Only estimate utility parameters
      mp1=mp;
      mp1.pnames_P={}; 
      theta1=struct2vec(mp1,mp.pnames_u);
      llfun=@(theta) zurcher.ll(data, mp1, theta);  
      [theta_hat,FVAL,EXITFLAG,OUTPUT1]=fminunc(llfun,theta1,optim_options);
      mp1=vec2struct(theta_hat,mp.pnames_u,mp);

      % Estimate RC, c and p
      mp2=mp1;
      mp2.pnames_P=mp.pnames_P; 
      llfun=@(theta) zurcher.ll(data, mp1, theta);  


      if ~isempty(mp.pnames_P);
        llfun=@(theta) zurcher.ll(data, mp2, theta);  
        theta2=struct2vec(mp2,[mp.pnames_u mp.pnames_P]);
        [theta_hat,FVAL,EXITFLAG,OUTPUT2]=fminunc(llfun,theta2,optim_options);
      else
        OUTPUT2.iterations=0;
        OUTPUT2.funcCount=0; 
      end

      timetoestimate=toc(outsidetimer);

      % Compute Variance-Covaiance matrix
      [f,g,h]=llfun(theta_hat);
      Avar=inv(h*samplesize);
      
      % Collect results
      results=vec2struct(theta_hat,{'RC', 'c'},mp);
      results.grad_direc=g*inv(h)*g';
      results.cputime=timetoestimate;
      results.converged   =   (EXITFLAG>=1 &  EXITFLAG<=3);

      results.MajorIter=OUTPUT1.iterations+OUTPUT2.iterations;
      results.funcCount=OUTPUT1.funcCount+OUTPUT2.funcCount;
      results.llval   =   -FVAL*samplesize;
    end % end of nfxp.estim
  end % end of methods
    
    
end % end of estim class