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
% modified by John Rust on 3/20/2021 to do double-sided numerical derivs for higher accuracy

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
%f0 = fhandle(x0);
%n = numel(f0);
k = numel(x0);
%grdd = zeros(n,k);
grddu=[];
grddl=[];

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
    
    xdhu = x0+dh;
    xdhl = x0-dh;
    dh = xdhu-xdhl;     
    
    reshapex0=zeros(k,k);
    for j=1:k
        reshapex0(:,j) = x0(:,1);
    end;
    
    argu=reshapex0;
    argl=reshapex0;
    for j=1:k
        argu(j,j) = xdhu(j,1);
        argl(j,j) = xdhl(j,1);
    end;
fprintf('gradp x0  argu argl \n');
    
    for i=1:k
x0
argu(:,i)
argl(:,i)
        gu=fhandle(argu(:,i));
        grddu=[grddu gu];
        gl=fhandle(argl(:,i));
        grddl=[grddl gl];
    end;

    %f02mat = (ones(k,1)*f0')';
    n=numel(gu);
    dh2mat = (ones(n,1)*dh');

    g = (grddu-grddl)./(dh2mat);
    

end
