%
% FDHess
%
% Numerical approximation of the Hessian matrix
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% function H = FDHess(x, Prob, gx, varargin)
%
% Sparsity pattern in Prob.HessPattern is used.
% Preprocessed efficient sparsity input in Prob.HessIx is utilized
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Nov 24, 1998.  Last modified Dec 10, 2004.
%

function H = FDHess(x, Prob, gx, varargin)

if nargin < 3, gx = []; end

x = x(:);

n = Prob.N;

x_L  = Prob.x_L(:);
if isempty(x_L)
   x_L = -inf*ones(n,1);
elseif length(x_L) < n
   x_L(length(x_L)+1:n,1) = -inf;
end
x_U  = Prob.x_U(:);
if isempty(x_U)
   x_U =  inf*ones(n,1);
elseif length(x_U) < n
   x_U(length(x_U)+1:n,1) =  inf;
end

if isempty(gx) 
   gx = nlp_g(x, Prob, varargin{:});
end
Pattern = Prob.HessPattern;
HessIx  = Prob.HessIx;

% HTOL is vector of relative intervals, See Pract. Opt. page 130
HTol    = Prob.GradTolH;
PriLev  = Prob.PriLevOpt;
NumDiff = Prob.NumDiff;

if isempty(HTol)
   HTol(1:n) = Prob.optParam.DiffInt;
end
if HTol(1) <= 0
   global gTol
   if ~isempty(gTol)
      % Use the result from the gradient estimation
      HTol      = 10*gTol;
   else
      h = Prob.optParam.DiffInt;
      if h < 0
         HTol = 1E-6*ones(n,1);
      else
         HTol = h*ones(n,1);
      end
   end
else
   HTol=HTol(:);
end
if length(HTol) < n
   HTol = HTol(1)*ones(n,1);
end

ix=find(x(1:n)+HTol(1:n) > x_U | x(1:n)+HTol(1:n) < x_L);

% If close to a bound, make the step in the opposite direction, changing sign
HTol(ix)=-HTol(ix);

if isempty(Pattern)
   H = zeros(n,n);
else
   [ix,iy,iv]=find(Pattern);
   if isempty(ix)
      % Pattern given, but all zeros. Nothing to calculate.
      H = zeros(n,n);
      return;
   end
   iv = double(iv);
end

if isempty(Pattern) & NumDiff <= 0
   for dim = 1:n
      Prob.FDVar = dim;
      z          = x(dim);
      h          = HTol(dim)*(1+abs(z));
      x(dim)     = z + h;
      gx_ph      = nlp_g(x,Prob, varargin{:});
      x(dim)     = z;
      H(:,dim)   = (gx_ph-gx)/h;
   end   
   H=0.5*(H+H'); % Make perfectly symmetric
elseif isempty(Pattern)
   % Two levels of differentiation, use f(x) values instead
   f             = nlp_f(x,Prob,varargin{:});
   for i = 1:n
      z          = x(i);
      Prob.FDVar = i;
      h          = HTol(i)*(1+abs(z));
      x(i)       = z + h;
      f1         = nlp_f(x,Prob,varargin{:});
      x(i)       = x(i) + h;
      f11        = nlp_f(x,Prob,varargin{:});
      H(i,i)     = ((f-f1)+(f11-f1))/h^2;
      for j = i+1:n
         x(i)       = z + h;
         v          = x(j);
         Prob.FDVar = [i,j];
         h2         = HTol(j)*(1+abs(v));
         x(j)       = v + h2;
         f12        = nlp_f(x,Prob,varargin{:});
         x(i)       = z;
         Prob.FDVar = j;
         f2         = nlp_f(x,Prob,varargin{:});
         x(j)       = v;
         hh         = h*h2;
         H(j,i)     = ((f12-f2)+(f-f1))/(h*h2);
         H(i,j)     = H(j,i);
      end
   end
elseif ~isempty(HessIx)
   ixS=[]; iyS=[]; ivS=[];
   mx = max(HessIx);
   for k = 1:mx
      CI          = find(HessIx==k);
      Prob.FDVar  = CI;
      z           = x(CI);
      h           = HTol(CI).*(1+abs(z));
      if NumDiff > 0 % Also the gradient is estimated numerically
         Prob.g_k = gx; 
         for j = 1:length(CI)
             % only estimate NaN elements in nlp_g
             Prob.g_k(ix(find(CI(j)==iy))) = NaN;  
         end
         Prob.cols = find(isnan(Prob.g_k));
      else
         if length(CI) == 1
            Prob.cols = ix(find(CI==iy));
         else
            Prob.cols = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
         end
      end
      x(CI)       = z + h;
      gx_ph       = nlp_g(x,Prob, varargin{:});
      x(CI)       = z;
      for j = 1:length(CI)
          iz      = find(CI(j)==iy);
          ixz     = ix(iz);
          tmp     = (gx_ph(ixz)-gx(ixz))/h(j);
          iv(iz)  = tmp;
          t       = find(ixz~=CI(j));
          if ~isempty(t)
             ixS  = [ixS;CI(j)*ones(length(t),1)];
             iyS  = [iyS;ixz(t)];
             ivS  = [ivS;tmp(t)];
          end
      end
   end
   H=sparse([ix;ixS],[iy;iyS],[iv;ivS],n,n);
elseif NumDiff > 0
   % Assumption. Only upper rectangle is stored in Pattern
   % Two levels of differentiation, use f(x) values instead
   [ix,iy,iv]=find(Pattern);
   ixS = zeros(length(ix),1); 
   iyS = zeros(length(ix),1); 
   ivS = zeros(length(ix),1);
   nS = 0;
   f                = nlp_f(x,Prob,varargin{:});
   if n == 1
      ic          = 1;
   else
      ic          = find(any(Pattern));
   end
   for ii = 1:length(ic)
      i          = ic(ii);
      iz         = find(i==iy);
      ixz        = ix(iz);
      m          = length(ixz)-1;
      Prob.FDVar = i;
      z          = x(i);
      h          = HTol(i)*(1+abs(z));
      x(i)       = z + h;
      f1         = nlp_f(x,Prob,varargin{:});
      x(i)       = x(i) + h;
      f11        = nlp_f(x,Prob,varargin{:});
      y          = zeros(m+1,1);
      y(end)     = ((f-f1)+(f11-f1))/h^2;
      for jj = 1:m
         j          = ixz(jj);
         x(i)       = z + h;
         v          = x(j);
         Prob.FDVar = [i,j];
         h2         = HTol(j)*(1+abs(v));
         x(j)       = v + h2;
         f12        = nlp_f(x,Prob,varargin{:});
         x(i)       = z;
         Prob.FDVar = j;
         f2         = nlp_f(x,Prob,varargin{:});
         x(j)       = v;
         hh         = h*h2;
         y(jj)      = ((f12-f2)+(f-f1))/(h*h2);
      end
      x(i)          = z;
      iv(iz)        = y;
      if m > 0
         ixS(nS+1:nS+m) = i;
         iyS(nS+1:nS+m) = ixz(1:m);
         ivS(nS+1:nS+m) = y(1:m);
         nS = nS + m;
      end
   end   
   H=sparse([ix;ixS(1:nS)],[iy;iyS(1:nS)],[iv;ivS(1:nS)],n,n);
else
   % Assumption. Only upper rectangle is stored in Pattern
   [ix,iy,iv]=find(Pattern);
   ixS=[]; iyS=[]; ivS=[];
   for dim = 1:n
      iz            = find(dim==iy);
      if ~isempty(iz)
         Prob.FDVar = dim;
         z          = x(dim);
         h          = HTol(dim)*(1+abs(z));
         ixz        = ix(iz);
         Prob.cols  = ixz;
         x(dim)     = z + h;
         if NumDiff > 0 % Also the gradient is estimated numerically
            Prob.g_k      = gx; 
            Prob.g_k(ixz) = NaN;  % only estimate NaN elements in nlp_g
         end
         gx_ph      = nlp_g(x,Prob, varargin{:});
         x(dim)     = z;
         tmp        = (gx_ph(ixz)-gx(ixz))/h;
         iv(iz)     = tmp;
         t          = find(ixz~=dim);
         if ~isempty(t)
            ixS=[ixS;dim*ones(length(t),1)];
            iyS=[iyS;ixz(t)];
            ivS=[ivS;tmp(t)];
         end
      end
   end   
   H=sparse([ix;ixS],[iy;iyS],[iv;ivS],n,n);
end

% MODIFICATION LOG
%
% 990306  hkh  Safeguard against x slightly outside bounds
% 000911  hkh  Major revision
% 001022  hkh  Must use n in sparse command to get right dim of H
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 030206  ango Check for all zeros in HessPattern
% 040125  hkh  Set double(iv) to avoid treatment at logical array in Matlab 6.5
% 040407  hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040407  hkh  If HessIx defined, use a more efficient method with less calls
% 040407  hkh  Set Prob.g_k with NaN, more efficient two-level differentiation
% 040407  hkh  Assume and utilize only upper triagonal of HessPattern
% 040412  hkh  New algorithm using f(x) if doing differentiation twice
% 040413  hkh  Check length of x_L and x_U
% 041210  hkh  Use 10*gTol, not gTol
