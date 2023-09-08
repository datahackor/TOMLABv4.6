%		ls_rJ.m
%
% function [Mode, r_k, J_k]=ls_rJ(x, Prob, Mode, nState)
%
% ls_rJ returns both the residual r(x) and the Jacobian J(x) for a
% nonlinear least squares problem and other vector valued problems.
% 
% ls_rJ calls the Tomlab gateway routines nlp_r, that returns the
% residual, and nlp_J, that returns the Jacobian matrix.
%
% Mode and nState is sent as Prob.Mode, Prob.nState to nlp_c and nlp_dc.
%
% Mode = 0 Assign function values
% Mode = 1 Assign known derivatives, unknown set as -11111 (=missing value)
% Mode = 2 Assign function and known derivatives, unknown derivatives -11111
%
% nState = 1         First call
% nState = 0         Other calls
% nState = 2+Inform  Last call, see Inform parameter for the solver
%
% ls_rJ is used when calling NLSSOL.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Mar 21, 1997.   Last modified Apr 9, 2002.
%

function [Mode, r_k, J_k]=ls_rJ(x, Prob, Mode, nState, dummy)

nargin;

if isfield(Prob.USER,'rJ')
   if Mode==0
      [Mode,r_k]       = feval(Prob.USER.rJ, x, Prob, Mode, nState);
   else
      [Mode, r_k, J_k] = feval(Prob.USER.rJ, x, Prob, Mode, nState);
   end

else
   Prob.Mode   = Mode;
   Prob.nState = nState;

   r_k = nlp_r(x, Prob);

   if Mode > 0
      J_k = full(nlp_J(x, Prob));

      %l = Prob.SOL.optPar(39);
      %if l==1 | l==3 | l==-999
      %   J_k = full(nlp_J(x, Prob));
      %else
      %   J_k = [];
      %end
   end
end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990707  hkh  Make full J_k
% 000918  hkh  New design for f90-MEX interface
% 020409  hkh  Adding Mode to Prob.Mode, nState to Prob.nState
% 020409  hkh  Bug - nstate was used instead of nState in rJ call

