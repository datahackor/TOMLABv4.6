%
% function f = qp_f(x, Prob)
%
% Compute function value to quadratic problem P
%
% x      Point x where f(x) is evaluated 
% Prob   Problem structure
% f      Function value, f(x).  f(x) = 0.5 * x'*F*x + c'*x;
%
% F*x is stored in global QP_Fx with corresponding x in QP_x, used by qp_g
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written April 28, 1995.  Last modified Aug 29, 1999.
%

function f = qp_f(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic objective function

if isempty(x)
   f = Inf;
else
   QP_x  = x;
   if isempty(Prob.QP.F)
      % LP Problem
      QP_Fx = [];
      if isempty(Prob.QP.c)
         f = 0;
      else
         f = Prob.QP.c(:)'*x;
      end
   else
      % QP Problem
      QP_Fx = Prob.QP.F*x;
      if isempty(Prob.QP.c)
         f = 0.5 * (x'*QP_Fx);
      else
         f = 0.5 * (x'*QP_Fx) + Prob.QP.c(:)'*x;
      end
   end
end

