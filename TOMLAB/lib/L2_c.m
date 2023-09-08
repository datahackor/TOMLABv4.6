%		L2_c.m
%
% function cx=L2_c(x, Prob, varargin)
%
% L2_c computes the L2 constraints c and r in the point x 
% r is the residuals in the original formulation: min 0.5 * r'*r
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.2.0$
% Written Apr 13, 1999. Last modified Jan 26, 2004.

function cx=L2_c(x, Prob, varargin)

m = Prob.L2.m;
n = length(x) - m;

Prob.USER.c  = Prob.L2.c;
Prob.USER.dc = Prob.L2.dc;
args         = Prob.L2.args;

Prob.x_0 = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L = Prob.x_L(1:n);
Prob.x_U = Prob.x_U(1:n);
Prob.N   = n;

global NARG
NARG(7:8) = args(7:8);

Prob.USER.r=Prob.L2.r;
Prob.USER.J=Prob.L2.J;

% Send also x(n+1:n+m) to nlp_r. nlp_r calls user routine only with x(1:n),
% where n is Prob.N

r  = nlp_r(x,Prob, varargin{:});
m  = length(r);

if isempty(Prob.USER.c)
   cx = [];
else
   % Must save and reset global variables
   global n_c NLP_xc NLP_c  
   n_c1    = n_c;
   NARG1   = NARG;
   NLP_xc1 = NLP_xc;
   NLP_c1  = NLP_c;
   n_c     = 0;
   NARG    = args;
   NLP_xc  = [];
   NLP_c   = [];

   cx = nlp_c(x,Prob, varargin{:});

   n_c     = n_c1;
   NARG    = NARG1;
   NLP_xc  = NLP_xc1;
   NLP_c   = NLP_c1;
end

cx = [cx;r-x(n+1:n+m)];

% MODIFICATION LOG
%
% 020413  hkh  Written
% 040126  hkh  Missing semi colon

