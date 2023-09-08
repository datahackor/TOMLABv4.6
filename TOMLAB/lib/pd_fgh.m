%		pd_fgh.m
%
% function [f, g, H] = pd_fgh(x, Prob, varargin)
%
% pd_fgh calls the TOMLAB gateway functions nlp_f, nlp_g and nlp_H unless
% the problem is an LP problem.
%
% pd_fgh is used when calling pdco and pdsco
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.2.0$
% Written Jan 9, 2002.   Last modified Feb 9, 2003.
%

function [f, g, H] = pd_fgh(x, Prob)

probType = Prob.probType;

n = length(x);
N = Prob.N;

%if checkType('lp',probType) - takes too long to check every time
if probType == 8
   g = Prob.QP.c;
   f = g' * x(1:N);
   if n~=N
      g = [g;zeros(n-N,1)];
   end
   H = zeros(n,1);
else
   f = nlp_f(x(1:N), Prob);
   g = nlp_g(x(1:N), Prob);
   H = diag(nlp_H(x(1:N), Prob));
   if n~=N
      g = [g;zeros(n-N,1)];
      H = [H;zeros(n-N,1)];
   end
end

% MODIFICATION LOG:
%
% 030120 hkh  Only diagonal of Hessian returned
% 030120 hkh  Handle different number of variables because of slacks added
% 030209 hkh  Error computing f for LP, speed up probType test, avoid checkType
% 031204 ango Fix name of this file to pd_fgh