%
% FDcHess
%
% Numerical approximation of the nonlinear constraints Hessian
% matrix.
%
% Based on FDHess.m which implements a version of the algorithm FD
% in Practical Optimization, page 343.
%
% function d2c = FDcHess(x,lam,Prob,dc,varargin)
%
% Sparsity pattern in Prob.ConsPattern is used.
% Preprocessed efficient sparsity input in Prob.ConIx is utilized
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written Jan 13, 2003.   Last modified Apr 12, 2004.
%

function d2c = FDcHess(x,lam,Prob,dc,varargin)

if nargin < 4, dc = []; end

x = x(:);
if isempty(dc)
   dc = nlp_dc(x,Prob,varargin{:});
end

lam      = lam(:);
n        = Prob.N;
Pattern  = Prob.ConsPattern;
ConsDiff = Prob.ConsDiff;
ConIx    = Prob.ConIx;
HTol     = Prob.GradTolH;

PriLev = Prob.PriLevOpt;

% Step length, check Prob.GradTolH first

% ... then Prob.optParam.DiffInt, or default 1e-6
if isempty(HTol)
   HTol(1:n) = Prob.optParam.DiffInt;
end
if HTol(1) <= 0
   global gTol
   if ~isempty(gTol)
      % Use the result from the gradient estimation
      HTol      = gTol;
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
% Fix length, using first element if not n elements
if length(HTol) < n
   HTol = HTol(1)*ones(n,1);
end

% DANGEROUS for central differences
% If close to a bound, make step in opposite direction, changing sign
ix = find(x(1:n)+2*HTol(1:n) > Prob.x_U(:) | x(1:n)+2*HTol(1:n) < Prob.x_L(:));

HTol(ix)=-HTol(ix);

d2c = sparse(n,n);

if isempty(Pattern) & ConsDiff <= 0
   ix = find(lam');
   if ~isempty(ix)
      Prob.rows = ix;
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)     = z + h;
        dc_ph      = nlp_dc(x,Prob,varargin{:});
        x(i)     = z;
    
        % Consider only constraints for which lam is nonzero
        for k=ix
          d2c(:,i) = d2c(:,i)+ (lam(k)/h)*( dc_ph(k,:)-dc(k,:) )';
        end
      end
   end
   % Make symmetric matrix
   d2c=0.5*(d2c+d2c');
elseif isempty(Pattern)
   % Two levels of differentiation, use c(x) values instead
   ix = find(lam');
   if ~isempty(ix)
      Prob.rows    = ix;
      c            = nlp_c(x,Prob,varargin{:});
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)       = z + h;
        c1         = nlp_c(x,Prob,varargin{:});
        x(i)       = x(i) + h;
        c11        = nlp_c(x,Prob,varargin{:});
        for k=ix
          %d2c(i,i) = d2c(i,i)+(lam(k)/h^2)*((c1(k)-c(k))+(c(k)-c1b(k)));
          d2c(i,i) = d2c(i,i)+(lam(k)/h^2)*((c(k)-c1(k))+(c11(k)-c1(k)));
        end
        for j = i+1:n
           x(i)       = z + h;
           v          = x(j);
           Prob.FDVar = [i,j];
           h2         = HTol(j)*(1+abs(v));
           x(j)       = v + h2;
           c12        = nlp_c(x,Prob,varargin{:});
           x(i)       = z;
           Prob.FDVar = j;
           c2         = nlp_c(x,Prob,varargin{:});
           x(j)       = v;
           hh         = h*h2;
           % Consider only constraints for which lam is nonzero
           for k=ix
             d2c(j,i) = d2c(j,i)+(lam(k)/hh)*((c12(k)-c2(k))+(c(k)-c1(k)));
           end
           d2c(i,j)   = d2c(j,i);
        end
      end
   end
%elseif 1 & ~isempty(ConIx)
%   % Nonempty Prob.ConsPattern... 
%   iL = find(lam');
%   if ~isempty(iL)
%      ConIx  = findpatt(Pattern(iL,:));
%      [ix,iy] = find(Pattern(iL,:));
%      %[ix,iy] = find(Pattern);
%      mx      = max(ConIx);
%      for j = 1:mx 
%          CI          = find(ConIx==j);
%          Prob.FDVar  = find(ConIx > 0);
%          z           = x(CI);
%          h           = HTol(CI).*(1+abs(z));
%          x(CI)       = z + h;
%          dc_ph       = nlp_dc(x,Prob,varargin{:});
%          x(CI)       = z;
%          for k = 1:length(CI)
%              i       = CI(k);
%    
%              % Consider only constraints for which lam is nonzero
%              % Consider only variables that are dependent on constraint k
%              for r = iL
%                  %iz         = find(r==ix & i==iy);
%                  iz         = find(r==ix);
%                  I          = iy(iz);
%                  d2c(I,i) = d2c(I,i)+ lam(r)*( dc_ph(r,I)-dc(r,I) )'/h(k);
%              end
%          end
%      end
%   end
%   % Make symmetric matrix
%   d2c=0.5*(d2c+d2c');
elseif isempty(Pattern) & ConsDiff > 0
   % Two levels of differentiation, use c(x) values instead
   % Nonempty Prob.ConsPattern... 
   ix = find(lam');
   if ~isempty(ix)
      if length(ix) == 1
         ic          = find(Pattern(ix,:));
      else
         ic          = find(any(Pattern(ix,:)));
      end
      Prob.rows      = ix;
      c              = nlp_c(x,Prob,varargin{:});
      Prob.cols      = ic;
      for ii = 1:length(ic)
         i           = ic(ii);
         z           = x(i);
         Prob.FDVar  = i;
         h           = HTol(i)*(1+abs(z));
         x(i)        = z + h;
         c1          = nlp_c(x,Prob,varargin{:});
         x(i)        = x(i) + h;
         c11         = nlp_c(x,Prob,varargin{:});
         for k = ix
            d2c(i,i) = d2c(i,i)+(lam(k)/h^2)*((c(k)-c1(k))+(c11(k)-c1(k)));
         end
         for jj = ii+1:length(ic)
            j          = ic(jj);
            x(i)       = z + h;
            v          = x(j);
            Prob.FDVar = [i,j];
            h2         = HTol(j)*(1+abs(v));
            x(j)       = v + h2;
            c12        = nlp_c(x,Prob,varargin{:});
            x(i)       = z;
            Prob.FDVar = j;
            c2         = nlp_c(x,Prob,varargin{:});
            x(j)       = v;
            hh         = h*h2;
    
            % Consider only constraints for which lam is nonzero
            for k = ix
               d2c(j,i) = d2c(j,i)+(lam(k)/hh)*((c12(k)-c2(k))+(c(k)-c1(k)));
            end
            d2c(i,j)   = d2c(j,i);
         end
         %x(i)       = z;
      end
   end
else
   % Nonempty Prob.ConsPattern... 
   ix = find(lam');
   if ~isempty(ix)
      if length(ix) == 1
         ic = find(Pattern(ix,:));
      else
         ic = find(any(Pattern(ix,:)));
      end
      Prob.rows     = ix;
      Prob.cols     = ic;
      for i = ic
         z          = x(i);
         Prob.FDVar = i;
         h          = HTol(i)*(1+abs(z));
         x(i)       = z + h;
         dc_ph      = nlp_dc(x,Prob,varargin{:});
         x(i)       = z;
    
         % Consider only constraints for which lam is nonzero
         for k = ix
            %d2c(:,i) = d2c(:,i)+ lam(k)*( dc_ph(k,:)-dc(k,:) )'/h;
            d2c(ic,i) = d2c(ic,i)+ (lam(k)/h)*( dc_ph(k,ic)-dc(k,ic) )';
         end
      end
   end
   % Make symmetric matrix
   d2c=0.5*(d2c+d2c');
end



% MODIFICATION LOG
%
% 030113 ango Wrote file
% 030113 hkh  Use Prob.ConsPattern to find which variables to use
% 030113 hkh  Define d2c sparse
% 030114 ango Fix error in x+HTol < Prob.x_L (was -, not +)
% 030128 hkh  Compute find(Lam') outside loop
% 030220 hkh  If ix empty, avoid all computations
% 040408 hkh  Only consider rows in Pattern for nonzero lam
% 040412 hkh  New algorithm using f(x) if doing differentiation twice
