%		nlp_fg.m
%
% function [Mode, f, g]=nlp_fg(x, Prob, Mode, nState)
%
% nlp_fg calls the TOMLAB gateway function nlp_f,
% which computes the function value f(x) for the test problem P (Prob.P).
%
% nlp_fg also calls the TOMLAB gateway function nlp_g,
% which computes the gradient g at x for the test problem P.
%
% Mode and nState is sent as Prob.Mode, Prob.nState to nlp_f and nlp_g.
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
% nlp_fg is used when calling MINOS and SNOPT.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.2.1$
% Written Mar 5, 1997.    Last modified June 16, 2002.
%

function [Mode, f, g]=nlp_fg(x, Prob, Mode, nState)

%if Mode > 0
%   save snopt_x x
%end
%fprintf('nlp_fg:   Mode %d nState %d\n',Mode,nState);

if isfield(Prob.USER,'fg')
   if Mode==0
      [Mode,f]=feval(Prob.USER.fg,x,Prob,Mode,nState);
   else
      [Mode,f,g]=feval(Prob.USER.fg,x,Prob,Mode,nState);
   end

else
   % Mode and nState set into the Prob structure
   Prob.Mode = Mode;
   Prob.nState = nState;

   f=nlp_f(x, Prob );
   if Mode > 0
      if Prob.NumDiff < 6
         g = nlp_g(x, Prob );
      else
         g = [];
      end

      %l = Prob.SOL.optPar(39);
      %if any( l == [1 3 -999] )
      %   g=nlp_g(x, Prob );
      %else
      %   g=[];
      %end
   end
end
%if isnan(f) | isinf(f) | isempty(f)
%   f=1E10;
%   if Mode > 0
%      g(isnan(g) | isinf(g)) = 1E3;
%   end
%end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 000820  hkh  Revised for v3.0, changing arguments
% 000828  hkh  Handle numerical differentiation in MINOS, returning empty g
% 020409  hkh  Adding Mode to Prob.Mode, nState to Prob.nState
% 020616  hkh  No call to nlp_g if NumDiff == 6, because no check in nlp_g
