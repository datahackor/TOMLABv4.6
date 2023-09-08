%		ls_f.m
%
% Computes the objective function for the least squares problem
%
%          0.5 * r(x)' * r(x)
%
% function f = ls_f(x, Prob, varargin)
%
% ls_f calls the TOMLAB gateway routine nlp_r to evaluate the
% residual r(x).
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.2.0$
% Written May 14, 1995.  Last modified Dec 1, 2003.
%

function f = ls_f(x, Prob, varargin)

%nargin;

if isa(x,'fmad') & Prob.ADObj == 1
    global LS_x mad_f
    x     = getvalue(x);
    LS_x  = [];
    % global mad_r 
    % nlp_r(x, Prob, varargin{:});
    % f = 0.5*mad_r'*mad_r;
    r = nlp_r(x, Prob, varargin{:});
    f = 0.5 * r' * r;
    mad_f = f;
else
    r = nlp_r(x, Prob, varargin{:});

    f = 0.5 * r' * r;
end

%fprintf('f(x) %30.15f\n',f);

% MODIFICATION LOG
%
% 981023  hkh  Changed to use this routine only to compute the square sum
% 990626  hkh  Avoid feval
% 031201  hkh  Add AD handling for MAD
% 040901  med  getvalue lower case