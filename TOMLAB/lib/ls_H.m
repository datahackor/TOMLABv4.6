%		ls_H.m
%
% function H=ls_H(x, Prob, varargin)
%
% ls_H computes the Hessian to a nonlinear least squares problem
%
% First part is computed as J'*J. 
% 2nd order terms are added if available
%
% ls_H calls the TOMLAB gateway routines nlp_r, nlp_J, nlp_d2r
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.6.0$
% Written May 14, 1995.    Last modified Jan 21, 2005.
%

function H=ls_H(x, Prob, varargin)

J = nlp_J(x, Prob);

H = J' * J;

H2=nlp_d2r( x, Prob, [], J, varargin{:});
if isempty(H2) | any(size(H2)~=size(H))
   return
else
   H=H+H2;
end

% MODIFICATION LOG
%
% 981023 hkh  Simplify this routine, only J'*J, but also new call to get
%              2nd part of Hessian, if available.
% 990626 hkh  Avoid feval
% 050121 frhe Removed r evaluation. It is done in nlp_d2r.m.

