function P=logit(v, sigma, i)
	% logit.m: Returns a mxn vector of logit probabilities of the mxn matrix of values v
	% INPUTS: 
	% 	v: 		mxn matrix of values v
	% 	sigma: 	(optional) scalar (scale fctor on iid ev error). If not passed, sigma=1 is assumed.   
	% 	i: 		(optional) scale integer: to only return choice probability of alternative i
	%
	% OUTPUS: 
	% 	logsum: mxb lector of probailities values
	%
	% syntax:  
	% 	P=logit(v, sigma, i):   P is a vector of logit probabilities for alternative i 
	% 	P=logit(v, sigma);      P is a matrix of logit probabilities for all alternatives
	% 	P=logit(v):             P is a matrix of logit probabilities for all alternatives assuming sigma=1
	%
	% See also:
	% 	logsum
	    
	if nargin<2
		sigma=1
	end
	
	if sigma==0;
		[~, idxmax]=nanmax(v,[], 2);
		P=zeros(size(v));
		P(:,idxmax)=1;
		return
	end
	maxv=nanmax(v,[], 2);
	expv=exp(bsxfun(@minus, v, maxv)/sigma);
	if nargin>2
			P=bsxfun(@rdivide,expv(:,i), nansum(expv,2));
	else
		P=bsxfun(@rdivide,expv, nansum(expv,2));
	end
	P(isnan(P))=0;
end % end logit
