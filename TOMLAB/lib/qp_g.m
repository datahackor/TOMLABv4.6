% function g = qp_g(x, Prob)
%
% Compute gradient to quadratic problem P
%
% x      Point x where g(x) is evaluated 
% Prob   Problem structure
% g      Gradient, g(x).  g(x) = F*x + c;
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1996-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Sept 10, 1996.   Last modified Sep 5, 1999.
%

function g = qp_g(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic gradient

if isempty(x)
   g = [];
   return
end

c=full(Prob.QP.c);

if isempty(c)
   c=zeros(length(x),1);
else
   c=c(:);
end

if ~isempty(QP_Fx) & length(x)==length(QP_x)
   if all(QP_x==x)
         g = QP_Fx + c;
   else
         g = Prob.QP.F*x + c;
   end
else
   if isempty(QP_Fx)
      g = c;
   else
      g = Prob.QP.F*x + c;
   end
end

