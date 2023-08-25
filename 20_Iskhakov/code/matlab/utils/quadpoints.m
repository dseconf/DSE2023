% quadpoints.m: Matlab function to return n  quadrature abscissa and weights 
%               over interval [lbnd,ubnd]. Adapted from Numerical Recipes by
%               John Rust, Yale University, February, 1997, in Gauss
%
%               Harry J. Paarsch, University of Iowa, December 2003

  function [x,w]=quadpoints(n,lbnd,ubnd);
	
%  local x2,x1,x,w,eps,m,xm,xl,i,z,z1,p1,p2,j,p3,pp;

   x2  = ubnd;
   x1  = lbnd;
   x   = zeros(n,1);
   w   = x;
   EPS = 3.e-14;
   m   = floor((n+1)/2);
   xm  = (x2+x1)/2;
   xl  = (x2-x1)/2;
   i = 1; 
   z1 = 1.e99;
   while (i <= m);
     z  = cos(pi*(i-.25)/(n+.5));
     while (abs(z-z1)>EPS); 
       p1 = 1;
       p2 = 0;
       j=1; 
       while (j <= n);
	 p3 = p2;
	 p2 = p1;
	 p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;
         j=j+1; 
       end; 
       pp = n*(z*p1-p2)/(z*z-1);
       z1 = z;
       z  = z1 - p1/pp;
     end;
     x(i)     = xm - xl*z;
     x(n+1-i) = xm + xl*z;
     w(i)     = 2*xl/((1-z*z)*pp*pp);
     w(n+1-i) =w(i);
     i = i+1;
   end; 
