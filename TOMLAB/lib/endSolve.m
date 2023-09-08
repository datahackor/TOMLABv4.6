% Routine to be called after optimization.
%
% function Result=endSolve(Prob,Result)
%
% Catch up search directions and line search steps (optionally).
% Computing CPU time and real time elapsed
% 
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.2.1$
% Written June 1, 1997.   Last modified Apr 7, 2004.
%

function Result=endSolve(Prob,Result)

%global TIME0 TIME1

if isempty(Result)
   return
end
TIME0 = Prob.TIME0;

if ~isempty(TIME0)
   Result.CPUtime = cputime-TIME0;
else
   Result.CPUtime = NaN;
end
%Result.CPUtime=CPUtime;
TIME1 = Prob.TIME1;
if ~isempty(TIME1)
   Result.REALtime = etime(clock,TIME1);
else
   Result.REALtime = NaN;
end
%Result.REALtime=REALtime;

% --------
% endbuild
% --------
% Create search directions, step lengths and variable limits for
% solver routines. Call aVbuild and reset flag BUILDP.
%
%

   global alphaV p_dx F_X BUILDP
   if BUILDP == 1
      alphaV=aVbuild(alphaV,size(p_dx,2));
   end
   BUILDP=0;

% End of endbuild script

if isfield(Prob.LS,'SepAlg')
   if Prob.LS.SepAlg > 0
      global SEP_z SEP_Jz
      Result.SepLS.z = SEP_z;
      Result.SepLS.Jz= SEP_Jz;
   end
end

% Save iterations in result struct
% global p_dx alphaV F_X
global X_min X_max

Result.p_dx=p_dx;
Result.alphaV=alphaV;
Result.x_min=X_min;
Result.x_max=X_max;
Result.F_X=F_X;

% Put Prob into Result struct if this is not already done
if ~isfield(Result,'Prob')
   Result.Prob=Prob;
end

%global solvType
%Result.solvType=solvType;


global n_f n_g n_H
global n_c n_dc n_d2c
global n_r n_J 
%[n_f, n_g, n_H, n_c, n_r, n_J]

if isempty(Result.FuncEv)    | Result.FuncEv==0
   Result.FuncEv=n_f;
end
if isempty(Result.GradEv)    | Result.GradEv==0
   Result.GradEv=n_g;
end
if isempty(Result.HessEv)    | Result.HessEv==0
   Result.HessEv=n_H;
end
if isempty(Result.ConstrEv)  | Result.ConstrEv==0
   Result.ConstrEv=n_c;
end
if isempty(Result.ConJacEv)  | Result.ConJacEv==0
   Result.ConJacEv=n_dc;
end
if isempty(Result.ConHessEv)  | Result.ConHessEv==0
   Result.ConHessEv=n_d2c;
end
if isempty(Result.ResEv)     | Result.ResEv==0
   Result.ResEv=n_r;
end
if isempty(Result.JacEv)     | Result.JacEv==0
   Result.JacEv=n_J;
end

% MODIFICATION LOG
%
% 020702 hkh Return directly if Result is empty
% 040402 hkh Use fields TIME0 and TIME1 in Prob, instead of globals
% 040407 hkh Add code to set ConJacEv and ConHessEv
