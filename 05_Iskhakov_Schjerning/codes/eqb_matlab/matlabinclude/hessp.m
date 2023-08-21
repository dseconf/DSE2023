
%  Purpose:  Computes the matrix of second partial derivatives
%            (Hessian matrix) of a function defined by a procedure.
%
%  Format:   h = hessp(@fhandle,x0);
%
%  Inputs:   @fhandle     pointer to a single-valued function f(x), defined
%                         as a procedure, taking a single Kx1 vector
%                         argument (f:Kx1 -> 1x1).  
%
%            x0           Kx1 vector specifying the point at which the Hessian
%                         of f(x) is to be computed.
%
%  Output:   h            KxK matrix of second derivatives of f with respect
%                         to x at x0. This matrix will be symmetric.
%

function h = hessp(fhandle,x0)

% CHECK FOR COMPLEX INPUT 
realcheck = isreal(x0); 
if realcheck == 0
    error('ERROR: Not implemented for complex matrices.');
else
  x0 = real(x0);
end;

% INITIALIZATIONS
    k = size(x0,1);
    hessian = zeros(k,k);
    grdd = zeros(k,1);
    eps = 6.0554544523933429e-6;
    
    
% COMPUTATION OF STEPSIZE (dh) 
ax0 = abs(x0);
    if x0 ~= 0
        dax0 = x0./ax0;
    else
        dax0 = 1;
    end;
    max0= max([ax0,(1e-2)*ones(size(x0,1),1)],[],2);  
    dh = eps*max0.*dax0;
 
    xdh = x0+dh;
    dh = xdh-x0;
    dh2mat=(ones(k,1)*dh')';
    ee = eye(k).*dh2mat;
    
% COMPUTATION OF f0=f(x0) 
    f0 = sum(fhandle(x0),1);

% COMPUTE FORWARD STEP 
    for i=1:k

        grdd(i,1) = sum(fhandle(x0+ee(:,i)),1);

    end;

% COMPUTE "DOUBLE" FORWARD STEP   
     for i=1:k
        for j=1:k

            hessian(i,j) = sum(fhandle(x0+(ee(:,i)+ee(:,j))),1);
            if i ~= j
                hessian(j,i) = hessian(i,j);
            end;

        end;
    end;

grdd2mat=(ones(k,1)*grdd')';
h = (((hessian - grdd2mat) - grdd2mat')+ f0)./ (dh2mat.*dh2mat');

end
