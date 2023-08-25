function inv = stdnormal_inv (x)
% STDNORMAL_INV  Quantile function of the standard normal distribution
%  INV = stdnormal_inv(X)
%  For each component of X, compute the quantile (the inverse of the
%  CDF) at X of the standard normal distribution.

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/stdnormal_inv.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright (C) 1995, 1996, 1997, 1998, 2000, 2002, 2005, 2006, 2007
%               Kurt Hornik
% Copyright (C) 2008 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

  if (nargin ~= 1)
    error('stdnormal_inv: you should provide one argument');
  end

  inv = sqrt (2) * erfinv (2 * x - 1);

end

