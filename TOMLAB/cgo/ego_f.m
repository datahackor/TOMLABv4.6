%
% function f = ego_f(x, Prob)
%
% Expected improvement function
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 2.3.0$
% Written Sep 29, 1998.   Last modified Sept 30, 2004.
%

function f = ego_f(x, Prob)

invR   = Prob.EGO.invR;
yMin   = Prob.EGO.yMin;

k      = Prob.EGO.k;
p      = Prob.EGO.p;

n      = size(invR,1);


if length(p) > 1
   r = exp(-( Prob.EGO.theta'*(abs( ...
       x(1:k)*ones(1,n)-Prob.EGO.X ).^(p*ones(1,n))) ))';
else
   r = exp(-( Prob.EGO.theta'*(abs( ...
       x(1:k)*ones(1,n)-Prob.EGO.X ).^p) ))';
end

%r = zeros(n,1);
%for i = 1:n
%   r(i) = exp(-(  theta'*(abs( x-X(:,i) )).^p ));
%end

Rr = invR * r;

%yHat = my + r'*invR*(Prob.EGO.y-my);

yHat = Prob.EGO.my + Rr'*(Prob.EGO.y-Prob.EGO.my);

tmp = 1-r'*Rr + (1-sum(Rr))^2/sum(sum(invR));

%if tmp <= 0
%   fprintf('Negative term in mean squared error of prediction %12.7e\n',tmp);
%   tmp = 1E-12;    % *********** NOTE ******************
%   % tmp = 1E-30;  % *********** NOTE ******************
%end

s = sqrt(Prob.EGO.sigma2*max(1E-300,tmp));

%  s=1;  % *********** NOTE ******************

% ML function. Minus sign to make it a minimization problem
if s > 0
   f = -( (yMin - yHat)*phi((yMin-yHat)/s) + s*fi((yMin-yHat)/s) );
else
   %s
   %keyboard
   %f = Inf;
   f = 1E10;
end
if ~isfinite(f), f=1E10; end

function y = phi(x) % Normal distribution function
y = 0.5  + 0.5 * erf(x/sqrt(2));

function y = fi(x) % Normal density function
y = 1/sqrt(2*pi)*exp(-x^2/2);


% MODIFICATIONS
% 040930  hkh  Safe guard, set f=1E10 if not finite
