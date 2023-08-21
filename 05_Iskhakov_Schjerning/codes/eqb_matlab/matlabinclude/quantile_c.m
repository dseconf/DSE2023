function [Y, FX, Xgrid]=quantile_c(X, count_X, P, n_Xgrid)
	% [Y, FX, Xgrid] = quantile_c(X, count_X, P, n_Xgrid) returns quantiles of the values in X with counts C.  
	%    P is a scalar or a vector of cumulative probability values.  
	%	 Y is the same size as P, and Y(i) contains the P(i)-th quantile. 
	%    n_Xgrid (optional) integer scaler that gives number of grid-points where to approximate CDF 
	%	 and interpolate between to obtain quantiles (default is 100). Increase for precision, decrease for speed.   

	if nargin<3;
		n_quantiles=10;
		P=linspace(1/n_quantiles,1, n_quantiles)
	end
	if nargin<4;
		n_Xgrid=100;
	end

	Xgrid=linspace(min(X), max(X), n_Xgrid)';
	FX=nan(n_Xgrid,1);

	sumc=sum(count_X);
	for iy=1:n_Xgrid;
		FX(iy)=sum(count_X(X<=Xgrid(iy),:))/sumc;
	end
	
	[Fu,i_Fu,~] = unique(FX);
	Y=interp1(Fu, Xgrid(i_Fu),  P);
end