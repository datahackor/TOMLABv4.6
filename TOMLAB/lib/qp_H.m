% function H = qp_H(x, Prob)
%
% Compute Hessian to quadratic problem 
%
% x      Point x where H(x) is evaluated 
% Prob   Problem structure
% H      Hessian matrix, H(x) = F, in f(x) = 0.5 * x'*F*x + c'*x;
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1996-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Sept 10, 1996.   Last modified Aug 25, 1999.
%

function H = qp_H(x, Prob)

x=x(:);

if isempty(Prob.QP.F)
   H = zeros(length(x),length(x));
else
   H = Prob.QP.F;
end

