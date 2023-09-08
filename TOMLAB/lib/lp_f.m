%
% Compute objective function value for linear programming problem
%
% function f = lp_f(x, Prob)
%
% x      Point x where f(x) is evaluated 
% Prob   Problem structure
% f      Function value, f(x).  f(x) =  c'*x;
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Nov 5, 1998.   Last modified Sep 5, 1999.
%

function f = lp_f(x, Prob)

% Evaluate linear objective function

if isempty(x)
   f = Inf;
else
   if isempty(Prob.QP.c)
      f = 0;
   else
      f = Prob.QP.c(:)'*x(:);
   end
end
