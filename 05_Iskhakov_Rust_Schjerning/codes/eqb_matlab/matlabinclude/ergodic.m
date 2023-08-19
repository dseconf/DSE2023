function [ed,ded]=ergodic(p,dp);
  % ergodic.m: finds the invariant distribution for an NxN Markov transition probability, p, 
  %             and the gradient of the invariant distribution for an NxN Markov transition probability
  %             when the transition probability matrix depends on a vector of variables. 
  %
  %             ergodic can be called as without the second output and input 
  %             if you only want to compute an invariant distribution and not its gradient
  %
  %   INPUTS
  %     p       The transition probability matrix of dimension nxn
  %     dp      A 3 dimensional array of size nxnxk so it is a vector of derivatives of transition matrices 
  %             with respect to a k dimensional parameter vector theta
  %
  %   OUTPUTS
  %     ed      The ergodic distribution, a (1:n) rows vector that satisfies q=q*p where p is an nxn transition probability matrix
  %     ded:    The gradient of ed with respect to a vector of parameters theta where p depends on theta, p(theta)
  %             this will be returned as an nxk matrix where k is the dimension of the theta vector
  %
  %   SYNTAX
  %     [ed]=ergodic(p): 
  %     [ed, ded]=ergodic(p,dp): 
  %
  %  John Rust, Georgetown University, January 2019


     n=size(p,1);

     if (size(p,2)~= n);
        fprintf('Error: p must be a square matrix\n');
        ed=NaN;
     end;

     ap=eye(n)-p';
 
     ap=[ap; ones(1,n)];

     ap=[ap ones(n+1,1)];

     if (rank(ap) < n+1);
        fprintf('Error: transition matrix p is not ergodic\n');
        ed=NaN;
     end;

     ed=ones(n,1);

     ed0=[ed; 2];

     inva=inv(ap);

     ed=inva*ed0;

     ed=ed(1:n);  % this is the ergodic distribution, the first return of this function

     % now compute ded, the gradient of the ergodic distribution with respect to theta

     if nargout>=2

       k=size(dp,3);  % dimension of theta

       ded=zeros(n,k);

       % loop over each component of theta, the vector of variables that p depends on

         da=zeros(n+1,n+1);  
       for i=1:k;

         dpk=dp(:,:,i);  % this is an n x n matrix of derivatives of p with respect to theta(i)

         da(1:n,1:n)=-dpk';
         dedk=-(inva*da*inva)*ed0;
         ded(:,i)=dedk(1:n);

       end;
     end
   end

