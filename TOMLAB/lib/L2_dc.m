%		L2_dc.m
%
% function dc=L2_dc(x, Prob, varargin)
%
% L2_dc computes the gradient of the constraints c at the point x
% and the Jacobian matrix for the L2 residuals r(x).
% r is the residuals in the original formulation: min 0.5*r'*r
% The extra variables have derivatives -1 exactly 
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Apr 13, 2002. Last modified Nov 30, 2004.

function dc=L2_dc(x, Prob, varargin)

m = Prob.L2.m;
n = length(x) - m;

Prob.x_0 = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L = Prob.x_L(1:n);
Prob.x_U = Prob.x_U(1:n);
Prob.N   = n;

Prob.USER.r = Prob.L2.r;
Prob.USER.J = Prob.L2.J;
args        = Prob.L2.args;

global NARG
NARG(7:8) = args(7:8);

Prob.NumDiff = Prob.L2.NumDiff;
J = nlp_J(x,Prob, varargin{:});

Prob.USER.c  = Prob.L2.c;
Prob.USER.dc = Prob.L2.dc;

if isempty(Prob.USER.c)
   dc = [];
else
   % Must save and reset global variables
   global n_dc NLP_xc NLP_c  
   n_dc1   = n_dc;
   NARG1   = NARG;
   NLP_xc1 = NLP_xc;
   NLP_c1  = NLP_c;
   n_dc    = 0;
   NARG    = args;
   NLP_xc  = [];
   NLP_c   = [];

   Prob.ConsDiff = Prob.L2.ConsDiff;

   dc = nlp_dc(x,Prob, varargin{:});

   n_dc   = n_dc1;
   NARG   = NARG1;
   NLP_xc = NLP_xc1;
   NLP_c  = NLP_c1;
end

if isempty(dc)
   dc=sparse([J,-speye(m,m)]);
else
   dc=sparse([dc,zeros(size(dc,1),m);[J,-speye(m,m)]]);
end

% MODIFICATION LOG
%
% 020413  hkh  Written
% 020416  hkh  Set original Prob.NumDiff and Prob.ConsDiff before call
% 040126  hkh  Wrong global variable, n_dc not n_c should be used
% 040126  hkh  Field Prob.L1.ConsDiff should be Prob.L2.ConsDiff
% 041130  hkh  Wrong check for nonlinear constraints, check Prob.USER.c not dc

