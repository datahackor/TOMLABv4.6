%		nlp_cdc.m
%
% function [Mode, c, dc]=nlp_cdc(x, Prob, Mode, nState, needc)
%
% nlp_cdc calls the TOMLAB gateway function nlp_c, 
% which evaluates the constraints at x for the test problem P (Prob.P).
%
% It also calls the TOMLAB gateway function nlp_dc, 
% which computes the gradient for all constraints at x for the test problem P
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
%
% nlp_cdc is called from NPSOL and NLSSOL
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Mar 5, 1997.    Last modified Apr 9, 2002.
%

function [Mode, c, dc]=nlp_cdc(x, Prob, Mode, nState, needc)

nargin;

if isfield(Prob.USER,'cdc')
   if Mode == 0
      [Mode,c]=feval(Prob.USER.cdc,x,Prob,Mode,nState);
   else
      [Mode,c,dc]=feval(Prob.USER.cdc,x,Prob,Mode,nState);
      dc = full(dc);
   end
else
   Prob.Mode   = Mode;
   Prob.nState = nState;

   c = nlp_c( x, Prob);

   if Mode > 0
      dc=full(nlp_dc(x, Prob));

      %l = Prob.SOL.optPar(39);  % Check on Derivative Level
      %if l==1 | l==3 | l==-999
      %   % Must transpose to have each row correspond to a constraint
      %   dc=dc(:);
      %else
      %   dc=[];
      %end
   end
end
%if any(any(isnan(dc)))
%        keyboard
%end
%if any(any(isinf(dc)))
%        keyboard
%end

% Each row in dc now corresponds to a constraint.

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990628  hkh  Add init
% 990705  hkh  Make dc a full matrix. Add varargin in calls.
% 000915  hkh  Revision for new NPSOL and NLSSOL interfaces
% 020409  hkh  Adding Mode to Prob.Mode, nState to Prob.nState

