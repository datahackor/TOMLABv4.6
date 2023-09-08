% function HessPattern = estHessPattern(Prob,Trials,varargin)
% 
% estHessPattern estimates Prob.HessPattern, the sparsity pattern of the
% Hessian of the objective
%
% The routine generates Trials random initial points between lower and
% upper bound
%
% It first tries to call the analytic Hessian routine, if available.
% If no analytic Hessian, but analytic gradient routine, 
% it tries TOMLAB /MAD (if installed)
% Otherwise it estimates a numerical finite difference Hessian
%
% INPUT:
% Prob        The Tomlab problem structure
% Trials      Number of trials to increase the probability that no element
%             by chance is 0 in the point tried. Default 2
%
% OUTPUT:
% HessPattern A 0-1 n by n-matrix, sparse or dense, with the pattern of the
%             Hessian of the objective, see the description of Prob.HessPattern
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2004-2004 by Tomlab Optimization Inc., $Release: 4.5.0$
% Written Apr 13, 2004.  Last modified Nov 14, 2004.
%

function HessPattern = estHessPattern(Prob,Trials,varargin)
if nargin < 2
   Trials = [];
end
if isempty(Trials), Trials = 2; end

global NARG
Func  = Prob.USER.g;
FuncH = Prob.USER.H;
if isempty(NARG)
   if isempty(Func)
      p  = 0;
   else
      p  = xnargin(Func);
   end
   if isempty(FuncH)
      pH = 0;
   else
      pH = xnargin(FuncH);
   end
else
   p  = NARG(2);
   pH = NARG(3);
end
N   = Prob.N;
BIG  = 1000;
x_0  = Prob.x_0(:);
if isempty(x_0)
   x_0 = zeros(N,1);
end
x_L  = Prob.x_L(:);
nL   = length(x_L);
if nL < N
   x_L(nL+1:N) = x_0(nL+1:N)-BIG;
   x_L = x_L(:);
end
x_U  = Prob.x_U(:);
nU   = length(x_U);
if nU < N
   x_U(nU+1:N) = x_0(nL+1:N)+BIG;
   x_U = x_U(:);
end
xD              = x_U - x_L;
xD(isinf(xD))   = 2*BIG;
x_L(isinf(x_L)) = x_0(isinf(x_L))-BIG;

Prob.x_L = x_L;
Prob.x_U = x_U;

if Prob.NumDiff == 0 &  pH > 0
   % Try analytic H
   H = sparse(N,N);
   for i = 1:Trials
       x = x_L+rand(N,1).*xD;
       H_x = abs(nlp_H(x,Prob,varargin{:}));
       if isempty(H_x)
          H = [];
          break
       else
          H = H + H_x;
       end
   end
else
   H = [];
end
if isempty(H)
   H = sparse(N,N);
   if checkMAD(0)
      % Try MAD for dc
      try
         for i = 1:Trials
             x    = x_L+rand(N,1).*xD;
             if p > 2
                mad_H=feval(Func, fmad(x(1:N),speye(N)),Prob,varargin{:});
             elseif p > 1
                mad_H=feval(Func, fmad(x(1:N),speye(N)), Prob);
             else
                mad_H=feval(Func, fmad(x(1:N),speye(N)));
             end
             %c=getvalue(mad_c);
             H = H + abs(getinternalderivs(mad_H));
         end
      catch
          % Try FD for H
          if p > 0
             % Utilize analytic gradient routine
             Prob.NumDiff = -1;
          else
             Prob.NumDiff = 1;
          end
          for i = 1:Trials
              x = x_L+rand(N,1).*xD;
              f = nlp_f(x,Prob,varargin{:});
              g = nlp_g(x,Prob,varargin{:});
              H = H + abs(nlp_H(x,Prob,varargin{:}));
          end
      end
   else
      % Try FD for H
      if p > 0
         % Utilize analytic gradient routine
         Prob.NumDiff = -1;
      else
         Prob.NumDiff = 1;
      end
      for i = 1:Trials
          x = x_L+rand(N,1).*xD;
          f = nlp_f(x,Prob,varargin{:});
          g = nlp_g(x,Prob,varargin{:});
          H = H + abs(nlp_H(x,Prob,varargin{:}));
      end
   end
end

HessPattern = spones(H);


% MODIFICATION LOG
%
% 040413  hkh  Algorithm formulated and written
% 040901  med  getvalue lower case
% 041114  hkh  Safe guard for empty nlp_H Hessian
