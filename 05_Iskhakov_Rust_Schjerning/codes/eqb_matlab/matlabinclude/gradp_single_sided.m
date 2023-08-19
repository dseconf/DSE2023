function g = gradp(fhandle,x0)
%  Purpose:    Computes the gradient vector or matrix (Jacobian) of a
%              vector-valued function that has been defined in a procedure.
%              Single-sided (forward difference) gradients are computed.
%
%  Format:     g = gradp(@fhandle,x0);
%
%  Input:      @fhandle    scalar, m-file function pointer to a vector-valued function:
%
%                                          f:Kx1 -> Nx1
%
%              x0           Kx1 vector of points at which to compute gradient.
%
%  Output:     g            NxK matrix containing the gradients of f with respect
%                           to the variable x at x0.

% CHECK FOR COMPLEX INPUT
realcheck = isreal(x0); 
if realcheck == 0
    error('ERROR: Not implemented for complex matrices.');
else
  x0 = real(x0);
end;
if size(x0,1)>1 && size(x0,2)>1; 
    error('ERROR: Gradient evaluation for a matrix input not allowed!'); 
end; 
if size(x0,2)>1; 
    error('Error in gradient evaluation: Only implemented for x0 as a column vector!'); 
end; 
f0 = fhandle(x0);
n = numel(f0);
k = numel(x0);
grdd = zeros(n,k);

% COMPUTATION OF STEPSIZE (dh)
    ax0 = abs(x0);
    if x0~=0
        dax0 = x0./ax0;
    else
        dax0 = 1;
    end;
    max0= max([ax0,(1e-2)*ones(size(x0,1),1)],[],2); 
    dh = (1e-8)*max0.*dax0;
    dh = (1e-6)*max0.*dax0;
    dh = (1e-6)*max0.*dax0;
    
    xdh = x0+dh;
    dh = xdh-x0;     
    
    reshapex0=zeros(k,k);
    for j=1:k
        reshapex0(:,j) = x0(:,1);
    end;
    
    arg=reshapex0;
    for j=1:k
        arg(j,j) = xdh(j,1);
    end;
    
    for i=1:k
        grdd(:,i) = fhandle(arg(:,i));
    end;

    f02mat = (ones(k,1)*f0')';
    dh2mat = (ones(n,1)*dh');

    g = (grdd-f02mat)./(dh2mat);
    

end
