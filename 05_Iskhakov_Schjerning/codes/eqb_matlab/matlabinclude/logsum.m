function logsum=logsum(v, sigma)
	% logsum.m: Returns a mx1 vector of log-sum values of the mxn matrix of values v
	%
	% INPUTS: 
	% 	v      	m x n matrix of values v
	% 	sigma   scalar (scale fctor on iid ev error) 
	% 
	% OUTPUS: 
	% 	logsum: mx1 lector of logsum values 
	%
	% SYNTAX:  
	% 	P=logsum(v, sigma);     
	% 	P=logsum(v):
	%
	% See also:
	% 	logit            

	if nargin ==1
		sigma=1
	end

	maxv=nanmax(v,[], 2);
	if sigma==0
		logsum=maxv;
		return
	end
	v=v- maxv;
	logsum=maxv + sigma*log(nansum(exp(v/sigma),2));
end % end logsum
