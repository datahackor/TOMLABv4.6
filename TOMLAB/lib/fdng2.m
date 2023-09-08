%
% Numerical gradient approximation 
%
% The numerical approximation is computed by use of the Matlab spline
% routine (NumDiff=2) or the Spline Toolbox routines 
% csapi (NumDiff=4,SplineTol <  0), csaps (NumDiff=3,splineSmooth) or 
% spaps (NumDiff=4,SplineTol >= 0)
%
% splineTol    = Prob.optParam.splineTol
% splineSmooth = Prob.optParam.splineSmooth
%
% The numerical step size h is given by Prob.optParam.CentralDiff;
% 
% function g = fdng2(x, Prob, g, fx, varargin)
%
% If g is nonempty, estimate any elements of g set to NaN.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written Oct 27, 1998.  Last modified Apr 11, 2004.
%

function g = fdng2(x, Prob, g, fx, varargin)

n  = Prob.N;

if isempty(g)
   ALLg = 1;
   g  = zeros(n,1);
else
   ALLg = 0;
   ix = find(isnan(g));
   if isempty(ix), return, end
end

if nargin < 4, fx = []; end

x  = x(:);

if isempty(fx)
   fx = nlp_f(x, Prob, varargin{:});
end

NumDiff      = abs(Prob.NumDiff);
h            = Prob.optParam.CentralDiff;
splineSmooth = Prob.optParam.splineSmooth;
splineTol    = Prob.optParam.splineTol;

if ALLg
   for dim = 1:n
      Prob.FDVar = dim;
      z          = x(dim);
      x(dim)     = z + h;
      fx_ph      = nlp_f(x,Prob, varargin{:});
      x(dim)     = z - h;
      fx_mh      = nlp_f(x,Prob, varargin{:});
      x(dim)     = z;
   
      XX = [ z-h   z  z+h ];
      YY = [fx_mh fx fx_ph];
   
      if NumDiff == 2       % Use Matlab spline routine
         pp  = spline(XX,YY);
      elseif NumDiff == 3   % Use SPLINE Toolbox routine csaps.m
         pp  = csaps(XX,YY,splineSmooth);
      elseif splineTol >= 0 % Use SPLINE Toolbox routine spaps.m
         pp  = spaps(XX,YY,splineTol);
      else                  % Use SPLINE Toolbox routine csapi.m
         pp  = csapi(XX,YY);
      end
   
      g(dim) = ppval(splineDer(pp, 1), z );
      %g(dim) = fnval(fnder(pp, 1),z);
   end   
else
   for i = 1:length(ix)
      dim        = ix(i);
      Prob.FDVar = dim;
      z          = x(dim);
      x(dim)     = z + h;
      fx_ph      = nlp_f(x,Prob, varargin{:});
      x(dim)     = z - h;
      fx_mh      = nlp_f(x,Prob, varargin{:});
      x(dim)     = z;
   
      XX = [ z-h   z  z+h ];
      YY = [fx_mh fx fx_ph];
   
      if NumDiff == 2       % Use Matlab spline routine
         pp  = spline(XX,YY);
      elseif NumDiff == 3   % Use SPLINE Toolbox routine csaps.m
         pp  = csaps(XX,YY,splineSmooth);
      elseif splineTol >= 0 % Use SPLINE Toolbox routine spaps.m
         pp  = spaps(XX,YY,splineTol);
      else                  % Use SPLINE Toolbox routine csapi.m
         pp  = csapi(XX,YY);
      end
      
      g(dim) = ppval(splineDer(pp, 1), z );
      %g(dim) = fnval(fnder(pp, 1),z);
   end
end

function ppD = splineDer(pp, order)

[breaks,coeffs,L,K,D] = unmkpp(pp);
if K <= order
   ppD = mkpp([breaks(1) breaks(L+1)],zeros(D,1));
else
   Knew = K - order;
   for i=K-1:-1:Knew
       coeffs = coeffs.*repmat([i:-1:i-K+1],D*L,1);
   end
   ppD = mkpp(breaks, coeffs(:,1:Knew),D);
end

% MODIFICATION LOG
%
% 981027  hkh  Return empty g in case of bound error on x
% 981029  hkh  Error in comments. 
% 981110  mbk  Added call to SPLINE Toolbox routines csaps.m and spaps.m
%              if Prob.NumDiff equals 3 or 4.
% 981124  mbk  Use of feval when calling SPLINE Toolbox routines.
% 990306  hkh  Safeguard against x slightly outside bounds
% 990626  hkh  Avoid feval
% 000910  hkh  Written more efficiently, avoid arrays and zero multiplications
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 030127  hkh  Added g as input, estimate NaN elements
% 040407  hkh  NumDiff = 2 now using standard spline, splineTol<0 => csapi
% 040411  hkh  Send Prob.FDVar with the variable index changed


