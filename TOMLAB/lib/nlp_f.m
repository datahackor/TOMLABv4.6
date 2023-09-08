%
% function [f,Result] = nlp_f(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of function values f(x)
%
% nlp_f calls the routine Prob.USER.f
% either as f=feval(Prob.USER.f, x) or
%           f=feval(Prob.USER.f, x, Prob or
%           f=feval(Prob.USER.f, x, Prob,varargin)
% depending on the number of inputs
%
% pbuild(x,f) is called to build search directions and search steps
%
% The global counter variable n_f is incremented
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.1$
% Written Oct 10, 1998.   Last modified May 26, 2004.
%

function [f,Result] = nlp_f(x, Prob, varargin)

global n_f BUILDP NARG NARGO PartSep

% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f NLP_pSepIndex

% NOT USED global NLP_xc NLP_c 

x      = x(:);
Result = [];

if Prob.simType > 0
   f = sim_fc(x, Prob, varargin{:} );
   return
end

if ~isempty(NLP_x)
   if length(x)~=length(NLP_x)
      NLP_x=[];
   elseif PartSep
      if NLP_pSepIndex == Prob.PartSep.index
         if all(x == NLP_x) 
            f=NLP_f;
            return
         end
      end
   %elseif all(x == NLP_x) 
   %   f=NLP_f;
   %   return
   end
end

n_f  = n_f+1;
Func = Prob.USER.f;
N    = min(length(x),Prob.N);

if isempty(Func)
   NLP_x=[];
   NLP_f=[];
   f=NaN;
   return
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(1);
   end
   if isempty(NARGO)
      NARGO(1) = nargout(Func);
   end
   
   if Prob.ADObj == 1
      global mad_f
      if NARGO(1) > 1
         if p > 2
            [mad_f,Result]=feval(Func, fmad(x(1:N),speye(N)), Prob,varargin{:});
         elseif p > 1
            [mad_f,Result]=feval(Func, fmad(x(1:N),speye(N)), Prob);
         else
            [mad_f,Result]=feval(Func, fmad(x(1:N),speye(N)));
         end
      else
         if p > 2
            mad_f=feval(Func, fmad(x(1:N),speye(N)), Prob, varargin{:});
         elseif p > 1
            mad_f=feval(Func, fmad(x(1:N),speye(N)), Prob);
         else
            mad_f=feval(Func, fmad(x(1:N),speye(N)));
         end
      end
      f=getvalue(mad_f);
   elseif Prob.ADObj == 2
      global HPI ad_f ad_g ad_H
      if isempty(HPI)
         HPI=getHPI(Func,N,Prob);
         %fprintf('Automatic differentiation using HPI in ADMAT.\n');
      end
      [ad_f ad_g ad_H]=evalH(Func,x(1:N),Prob,HPI);
      f=ad_f;
   else
      BUILDP0=BUILDP;
      BUILDP=[];
      if NARGO(1) > 1
         if p > 2
            [f,Result]=feval(Func, x(1:N), Prob, varargin{:} );
         elseif p > 1
            [f,Result]=feval(Func, x(1:N), Prob);
         else
            [f,Result]=feval(Func, x(1:N));
         end
      else
         if p > 2
            f=feval(Func, x(1:N), Prob, varargin{:} );
         elseif p > 1
            f=feval(Func, x(1:N), Prob);
         else
            f=feval(Func, x(1:N));
         end
      end
      BUILDP=BUILDP0;
   end
end
if isnan(f), f=Inf; end
%if isnan(f)
%   error('nlp_f: The user function returned NaN. Safeguard code!');
%end

NLP_x=x;
NLP_f=f;

if PartSep
   % Number of partially separable functions.
   NLP_pSepIndex = Prob.PartSep.index;    % Function computed
end

if BUILDP > 0
   pbuild(x,f);
end

%xprint(x,'x:')
%f

% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation 
% 981027   hkh   Unnecessary return. No Prob.ctrl used
% 981029   hkh   Communicating NLP_x,f,c,dc to be used in numerical diffs
% 981102   hkh   Must check if partially separable function. Add to global.
% 981119   hkh   Check if x and NLP_x has same length
% 981120   hkh   Check if PartSep values are defined. 
%                Test if Prob.probType empty
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
%                Add printout about the use of automatic differentiation
% 020110   hkh   Set global BUILDP, avoide calls from nlp_r to pbuild
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020528   hkh   More efficient handling of partially separable functions
% 030524   hkh   Add handling of simulation problems
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 031213   hkh   Check f for NaN, set as Inf
% 040318   hkh   Move sim_fc first, before check on NLP_x
% 040409   hkh   More efficient call to MAD
% 040418   hkh   Return 2nd output Result if NARGO(1) > 1 (nargout > 1)
% 040426   hkh   Safeguard if NARGO not defined
% 040526   hkh   Use x(1:N) in function calls
% 040901   med   getvalue lower case