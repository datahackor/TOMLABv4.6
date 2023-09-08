%		ls_g.m
%
% Computes the gradient to a nonlinear least squares problem,
%
%            J(x)' * r(x)
%
% function g=ls_g(x, Prob, varargin)
%
% ls_g calls the TOMLAB gateway routines nlp_r, that returns the
% residual r(x), and nlp_J, that returns the Jacobian matrix J(x).
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.2.0$
% Written May 14, 1995.  Last modified Dec 1, 2003.
%

function g = ls_g(x, Prob, varargin)

nargin;

if isa(x,'fmad') & Prob.ADObj == -1
    global mad_J mad_r LS_x mad_g
    x = getvalue(x);
    LS_x = []; LS_xJ = [];
    Prob.ADObj = 1;
    nlp_r(x, Prob, varargin{:});
    Prob.ADObj = -1;
    nlp_J(x, Prob, varargin{:});
    g     = mad_J'*mad_r;
    mad_g = g;
else
    r = nlp_r(x, Prob, varargin{:}); 
    J = nlp_J(x, Prob, varargin{:});
    if isempty(r) | isempty(J)
        g = [];
    else
        g = J' * r;
    end
end
% MODIFICATION LOG
%
% 981023  hkh  Simplify this routine, only J'*g. Call nlp_r/nlp_J instead
% 990626  hkh  Avoid feval
% 031201  hkh  Adding AD handling for MAD
% 040901  med  getvalue lower case