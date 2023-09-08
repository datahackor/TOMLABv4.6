%		ls2_rJS.m
%
% function [r_k, J_k]=ls2_rJS(x, Prob, varargin)
%
% ls2_rJS returns both the residual r(x) and the Jacobian J(x) for a
% nonlinear least squares problem.
% 
% ls2_rJS calls the USER function that returns the
% residual, and possibly the user function that returns the Jacobian matrix.
% J_k may be sparse
%
% ls2_rJS is used when implementing the OPTIM TB 2.0 compatibility interface
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.0.3$
% Written July 28, 1999.    Last modified Feb 11, 2003.
%

function [r_k, J_k]=ls2_rJS(x, Prob, varargin)

nargin;
global n_r n_J NARG

r=Prob.USERX.r;
n_r = n_r + 1;

if NARG(7) > 1
   r_k = feval(r, x, Prob, varargin{:});
else
   r_k = feval(r, x);
end

if nargout > 1
   J=Prob.USERX.J;
   if ~isempty(J)
      n_J = n_J + 1;
      if NARG(8) > 1
         J_k = feval(J, x, Prob, varargin{:});
      else
         J_k = feval(J, x, Prob);
      end
   else
      J_k=[];
   end
end

% MODIFICATION LOG:
%
% 030211 hkh Revision for v4.0, using NARG

