% Routine to be called after optimization formulated in AMPL.
%
% function Result=postSolveAMPL(Result)
%
% Includes parts generated by AMPL Presolve in Result.
% 
% Marcus Edvall, Tomlab Optimization Inc, E-mail: medvall@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written July 22, 2003.   Last modified July 22, 2003.
%

function Result = postSolveAMPL(Result)

if isfield(Result,'Prob')
if isfield(Result.Prob, 'AMPL')
        if Result.Prob.AMPL.nzo > 0 % No modification should be made if the objective is a constant.
            Result.f_k = Result.f_k * Result.Prob.AMPL.objType + Result.Prob.AMPL.objConst;
        end
end, end

if isfield(Result,'Prob')
if isfield(Result.Prob, 'AMPL')    
    if Result.Prob.AMPL.solFile == 1 & Result.ExitFlag == 0
        % A .SOL file will be written, so that the solution can be read in
        % AMPL
        if Result.Prob.AMPL.sparseFlag == 1 
            spamfunc(sprintf('Problem: %s, \nSolver: %s \nExit Result: %s, \nObjective: %5.6f',...
		        Result.Name, Result.Solver, Result.ExitText(1,:), Result.f_k), Result.x_k, Result.v_k(length(Result.x_k)+1:end), Result.f_k);
        else
            amplfunc(sprintf('Problem: %s, \nSolver: %s \nExit Result: %s, \nObjective: %5.6f',...
		        Result.Name, Result.Solver, Result.ExitText(1,:), Result.f_k), Result.x_k, Result.v_k(length(Result.x_k)+1:end), Result.f_k);
        end
end, end, end

% MODIFICATION LOG:
% 030801  medvall Written