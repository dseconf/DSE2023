  % matgrad: function to compute jacobian of function f=fx(x) that returns matrix of outputs, f, from a vector inputs x.
  %
  %     SYNTAX: [g] = matgrad(@fhandle, x0)
  %     
  %     INPUT   : function handle on the form f=fx(x), 
  %                 where x is p x 1 vector of inputs and y is m x n matrix of outputs f
  %           
  %            x0 : p x 1 vector of inputs where to evaluate gradient
  %     
  %     OUTPUT  g : m x n x p array of gradients. If n=1, g is m x p
  % see also gradp, hessp

  function [g] = matgrad(fx, x0)

    function [fvec] = f_wrapper(fx, x0)
      f=fx(x0);
      fvec=reshape(f, numel(f),1); 
    end

    f=fx(x0);
    szf=size(f);
    if sum(szf(2:end))==1;
      szf=numel(f);
    end 
    fvec=@(x) f_wrapper(fx, x);
    g=gradp(fvec,x0);
    
    g=reshape(g, [szf, numel(x0)]);
  end