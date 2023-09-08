%
% FDrHess
%
% Numerical approximation of the residuals Hessian matrix.
%
%   sum(i=1:m) r_i*d2r_i
%
% Based on FDcHess.m which is based on FDHess.m which implements 
% a version of the algorithm FD in Practical Optimization, page 343.
%
% function d2r = FDrHess(x,Prob,r,J,varargin)
%
% Sparsity pattern in Prob.JacPattern is used.
% NOT preprocessed efficient sparsity input in Prob.ConIx is utilized
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written Jan 21, 2005.   Last modified Jan 21, 2005.
%

function d2r = FDrHess(x,Prob,r,J,varargin)

global gurka

gurka = gurka + 1;

if nargin < 5
  J = [];
  if nargin < 4
    r = [];
  end
end

x = x(:);

if isempty(J)
  J = nlp_J(x,Prob,varargin{:});
end
if isempty(r)
  r = nlp_r(x,Prob,varargin{:});
end
  
r        = r(:);
n        = Prob.N;
Pattern  = Prob.JacPattern;
NumDiff  = Prob.NumDiff;
%ConIx    = Prob.ConIx;
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

d2r = sparse(n,n);

if isempty(Pattern) & NumDiff <= 0
   ix = find(r');
   if ~isempty(ix)
      Prob.rows = ix;
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)       = z + h;
        J_ph       = nlp_J(x,Prob,varargin{:});
        x(i)       = z;
    
        % Consider only residuals for which r is nonzero
        for k=ix
          d2r(:,i) = d2r(:,i)+ (r(k)/h)*( J_ph(k,:)-J(k,:) )';
        end
      end
   end
   % Make symmetric matrix
   d2r=0.5*(d2r+d2r');
elseif isempty(Pattern)
   % Two levels of differentiation, use r(x) values instead
   ix = find(r');
   if ~isempty(ix)
      Prob.rows    = ix;
%      r            = nlp_r(x,Prob,varargin{:});
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)       = z + h;
        r1         = nlp_r(x,Prob,varargin{:});
        x(i)       = x(i) + h;
        r11        = nlp_r(x,Prob,varargin{:});
        for k=ix
          %d2c(i,i) = d2c(i,i)+(lam(k)/h^2)*((c1(k)-c(k))+(c(k)-c1b(k)));
          d2r(i,i) = d2r(i,i)+(r(k)/h^2)*((r(k)-r1(k))+(r11(k)-r1(k)));
        end
        for j = i+1:n
           x(i)       = z + h;
           v          = x(j);
           Prob.FDVar = [i,j];
           h2         = HTol(j)*(1+abs(v));
           x(j)       = v + h2;
           r12        = nlp_r(x,Prob,varargin{:});
           x(i)       = z;
           Prob.FDVar = j;
           r2         = nlp_r(x,Prob,varargin{:});
           x(j)       = v;
           hh         = h*h2;
           % Consider only residuals for which r is nonzero
           for k=ix
             d2r(j,i) = d2r(j,i)+(r(k)/hh)*((r12(k)-r2(k))+(r(k)-r1(k)));
           end
           d2r(i,j)   = d2r(j,i);
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
elseif ~isempty(Pattern) & NumDiff > 0
   % Two levels of differentiation, use r(x) values instead
   % Nonempty Prob.ConsPattern... 
   ix = find(r');
   if ~isempty(ix)
      if length(ix) == 1
         ic          = find(Pattern(ix,:));
      else
         ic          = find(any(Pattern(ix,:)));
      end
      Prob.rows      = ix;
%      r              = nlp_r(x,Prob,varargin{:});
      Prob.cols      = ic;
      for ii = 1:length(ic)
         i           = ic(ii);
         z           = x(i);
         Prob.FDVar  = i;
         h           = HTol(i)*(1+abs(z));
         x(i)        = z + h;
         r1          = nlp_r(x,Prob,varargin{:});
         x(i)        = x(i) + h;
         r11         = nlp_r(x,Prob,varargin{:});
         for k = ix
            d2r(i,i) = d2r(i,i)+(r(k)/h^2)*((r(k)-r1(k))+(r11(k)-r1(k)));
         end
         for jj = ii+1:length(ic)
            j          = ic(jj);
            x(i)       = z + h;
            v          = x(j);
            Prob.FDVar = [i,j];
            h2         = HTol(j)*(1+abs(v));
            x(j)       = v + h2;
            r12        = nlp_r(x,Prob,varargin{:});
            x(i)       = z;
            Prob.FDVar = j;
            r2         = nlp_r(x,Prob,varargin{:});
            x(j)       = v;
            hh         = h*h2;
    
            % Consider only residuals for which r is nonzero
            for k = ix
               d2r(j,i) = d2r(j,i)+(r(k)/hh)*((r12(k)-r2(k))+(r(k)-r1(k)));
            end
            d2r(i,j)   = d2r(j,i);
         end
         %x(i)       = z;
      end
   end
else
   % Nonempty Prob.Jacattern... 
   ix = find(r');
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
         J_ph      = nlp_J(x,Prob,varargin{:});
         x(i)       = z;
    
         % Consider only residuals for which r is nonzero
         for k = ix
            %d2c(:,i) = d2c(:,i)+ lam(k)*( dc_ph(k,:)-dc(k,:) )'/h;
            d2r(ic,i) = d2r(ic,i)+ (r(k)/h)*( J_ph(k,ic)-J(k,ic) )';
         end
      end
   end
   % Make symmetric matrix
   d2r=0.5*(d2r+d2r');
end

% MODIFICATION LOG
%
% 050121 frhe File created, based on FDcHess.m