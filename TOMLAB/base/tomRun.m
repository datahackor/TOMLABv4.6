%
% tomRun.m - General driver routine for TOMLAB
%
% If using the TOMLAB Quick format (TQ), call with:
%
%   function [Result] = tomRun(Solver, Prob, PriLev, ask);
%
% NOTE! The call Result = tomRun(Solver, Prob, [], 2); will be interpreted 
% as (for compatability reasons) PriLev = 2, ask = [] 
% (reverse order of parameters as in Tomlab 1.0 - 4.4)
%   
% The following call will also work (similar 6-input format as below)
%
%   function [Result] = tomRun(Solver, [], [], Prob, PriLev, ask);
%
% If using the TOMLAB Init File format, call with:
%
%   function [Result] = tomRun(Solver, probFile, probNumber, Prob, PriLev, ask);
%
% A third alternative is the call
%
%   function [Result] = tomRun(Solver, probType, probNumber, Prob, PriLev, ask);
%
% Then the default file for problems of type probType is used.
%
% If calling with tomRun; (no arguments), a list of available solvers is given
%
% If calling with tomRun(probType); 
%    a list of available solvers for probType is given
%
% tomRun checks if the second argument is a string, a structure, a number,
% or is empty, to determine which input format is used.
%
% if isempty(Solver),   tomRun is using the TOMLAB default solver for the
%                       problem type as default (calling GetSolver)
% if isempty(probFile), tomRun is using con_prob as default
%
%
% A problem available in the TOMLAB Init File format is defined using a call 
%           Prob=probInit(probFile, probNumber, ask, Prob)
%
% INPUT: (if [] is given or less parameters are given, default values are used)
%
% Solver     The name of the solver that should be used to optimize the problem.
%            If the Solver may run several different optimization algorithms,
%            then the values of Prob.Solver.Alg and Prob.Solver.SubAlg
%            determines which algorithm.
%            The Solver name is put in Prob.Solver.Name 
%
% probFile   User problem initialization file.  The GUI and meny system
%            is using the different probFile's defined by call to nameprob.
%            To make a new probFile permanent, run AddProblemFile.
%            If empty or left out, either probFile=Prob.probFile (if set) or
%            the default probFile returned from nameprob is used.
%          
% probNumber Problem number in probFile. 
%            If empty of left out, either probNumber=Prob.P (if set) or
%            otherwise probNumber=1.
%            When calling the probFile with probNumber=0, probFile must 
%            return a string matrix with the names of the problems defined.
%            (Used by the menu routines and GUI)
%
% Prob       Problem structure. Either define the structure with the
%            call:  Prob=probInit(probFile,probNumber,ask);
%            or set in Prob the parameters with special values.
%            See the manual for a description of the Prob structure
%
%            P        probNumber (YOU MUST SET THIS, IF SETTING PROB AS INPUT!)
%
%            Examples of other fields to set:
%            probFile probFile
%            uP       User problem parameters, which you can use when computing
%                     the functions. See the variable ask .  
%            optParam A substructure with optimization parameters. Default value
%                     from optParamSet(Solver,probType)
%            x_0      Starting point for the problem.
% PriLev     Print level when displaying the result of the optimization in
%            routine PrintResult.
%            =0 No output, =1 Final result, shorter version, 
%            =2 Final result, longer version, =3 All output
%            If isempty(PriLev), Prob.PriLev is used, if nonempty
%         
%            The printing level in the optimization solver is controlled
%            by setting the parameter Prob.PriLevOpt
%
% ask        If ask>=1: ask questions in probFile. If ask = 0: use defaults
%            If ask < 0: use values in user params uP if defined or defaults. 
%            If isempty(ask): If length(uP) > 0, ask=-1, else ask=0
%
% OUTPUT:
% Result Structure with optimization results. See manual for description
%        Result.Prob holds the input Prob structure
%
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Mar 2, 1998.    Last modified Jan 24, 2005.
%

%function [Result] = tomRun(Solver, probFile, probNumber, Prob, PriLev, ask);
%function [Result] = tomRun(Solver, Prob,  PriLev, ask);

function [Result] = tomRun(Solver, Prob, P3, P4, P5, P6)


global n_f n_g n_H    % Count of function, gradient, Hessian evaluations
global n_c n_dc n_d2c % Count of constraint, constr.grad, and 2nd der evals
global n_r n_J n_d2r  % Count of residual, Jacobian and 2nd der evaluations
global xGUI
global GlobalLevel
GlobalLevel = [];     % Initialize to zero depth for recursion
global xGUI % If GUI is used or not
%global solvType
global optType probType
if isempty(xGUI) | xGUI == 0 
   %solvType = [];
   probType = [];
   optType = [];
   xGUI=0; 
end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

if nargin < 2
   if nargin < 1
      for i=1:10
          fprintf('\nSolvers for problem type %d\n',i)
          z=SolverList(i);
          disp(z);
      end
   else
      if ~isstr(Solver)
         fprintf('\nSolvers for problem type %d\n',Solver)
         z=SolverList(max(1,min(10,round(Solver))));
         disp(z);
      else
         for i=1:10
             fprintf('\nSolvers for problem type %d\n',i)
             z=SolverList(i);
             disp(z);
         end
      end
   end
   return
   %fprintf('\n');
   %error('tomRun needs at least two parameters');
end

if nargin < 6
   P6=[];
   if nargin < 5
      P5=[];
      if nargin < 4
         P4=[];
         if nargin < 3
            P3=[];
end, end, end, end

if isstruct(Prob)
   TQ         = 1;
   if nargin == 4 & isempty(P3)
      PriLev     = P4;
      ask        = [];
   else
      PriLev     = P3;
      ask        = P4;
   end
   optType    = Prob.probType;
   probFile   = 0;  % Signal that a Prob struct is used
   probNumber = 1;
elseif isstr(Prob)
   TQ         = 0;
   probFile   = Prob;
   probNumber = P3;
   Prob       = P4;
   PriLev     = P5;
   ask        = P6;
   if isstruct(Prob)
      if isfield(Prob,'probType')
         optType    = Prob.probType;
      else
         optType    = 3;
      end
   else
      Prob = [];
      optType    = 3;
   end
elseif isempty(Prob)
   % Assume 4th argument is the Prob structure
   if nargin < 4
      error('When 2nd argument empty, tomRun needs 4 arguments');
   end
   if ~isstruct(P4)
      error('When 2nd argument empty, 4th argument must be a Prob structure');
   end
   TQ         = 1;
   Prob       = P4;
   PriLev     = P5;
   ask        = P6;
   optType    = Prob.probType;
   probFile   = 0;  % Signal that a Prob struct is used
   probNumber = 1;
else
   TQ         = 0;
   probType   = Prob;

   [DataFile,NameFile,DefFile,probTypeList]=nameprob(probType,0);

   probFile   = DataFile(1,:);
   if isempty(P3)
      probNumber = 1;
   else
      probNumber = P3;
   end
   Prob       = P4;
   PriLev     = P5;
   ask        = P6;
   optType    = probType;
end

if TQ == 0
   %[Prob, ask, PriLev, probFile, probNumber, x_0] = ...
   %      xxxRun(optType, Prob, ask, PriLev, probFile, probNumber);
%

   if isstruct(Prob)
      if ~isempty(probNumber)
         Prob.P=probNumber;
      end
      if isfield(Prob,'x_0')
         x_0=Prob.x_0;  % If input x_0 has correct length, use this as start
      else
         x_0=[];
      end
      if isfield(Prob,'P') & isempty(probNumber)
         probNumber=Prob.P;
      elseif ~isfield(Prob,'P') 
         Prob.P=-1;
      end
      if isfield(Prob,'probFile') & isempty(probFile)
         probFile=Prob.probFile;
      else
         Prob.probFile=[];
      end
      %NOT Prob=ProbCheck(Prob,optType);
      Prob.CHECK=0;
   else
      x_0=[];
   end

   if isempty(probFile)
      [DataFile,NameFile,DefFile,probTypV]=nameprob(optType);
      probFile=DataFile(DefFile,:);
   end

   if isempty(probNumber),  probNumber=1; end
   
   if isempty(PriLev)
      if isfield(Prob,'PriLev')
         if ~isempty(Prob.PriLev)
            PriLev=Prob.PriLev;
         end
      end
   end

   if isempty(PriLev), PriLev=2; end  

   if isempty(ask)
      if isstruct(Prob)
         if isfield(Prob,'uP')
            if isempty(Prob.uP)
               ask=0;
            else
               ask=-1;
            end
         else
            ask=0;
         end
      else
         ask=0;
      end
   end

   if isempty(probFile), probFile='con_prob';end
else
   x_0 = Prob.x_0;

   if isempty(PriLev)
      if isfield(Prob,'PriLev')
         if ~isempty(Prob.PriLev)
            PriLev=Prob.PriLev;
         end
      end
   end

   if isempty(PriLev), PriLev=0; end  
end

%if isempty(Solver) 
%   isType = checkType([],optType);
%   if isfield(Prob,'LargeScale')
%      Solver = GetSolver(isType,Prob.LargeScale);
%   else
%      Solver = GetSolver(isType);
%   end
%else
%   Solver=deblank(Solver);
%end

% Check if Prob structure already defined
if isfield(Prob,'MENU')
   MENU=Prob.MENU;
else
   MENU=0;
end
if isfield(Prob,'P')
   if Prob.P <= 0
      MENU=0;
   end
end

if ~ischar(probFile)
   if probFile==0
      MENU=1; % Accept the Prob structure as it is
      if isfield(Prob,'probType')
         probType = Prob.probType;
      end
   end
end

if MENU==0 
   % Define problem in structure Prob
   probFile=deblank(probFile);
   Prob = probInit (probFile, probNumber(1), ask, Prob);
   probType = Prob.probType;
end

if isempty(probType)
   probType=optType;
end

Prob = mkbound(Prob);

probNumber=Prob.P;

if ~isempty(x_0) & ~isfield(Prob,'N')
   Prob.N=length(x_0);
end

if ~isempty(x_0)
   if isempty(Prob.N), Prob.N=length(x_0); end
   if length(x_0)==Prob.N & Prob.N ~=0
      Prob.x_0 = x_0(:); % Should be column vector. Use input starting values
   end
end

if isempty(Solver) 
   isType = checkType([],probType);
   if isfield(Prob,'LargeScale')
      Solver = GetSolver(isType,Prob.LargeScale);
   else
      Solver = GetSolver(isType);
   end
else
   Solver=deblank(Solver);
end

if 0 % The following part should not be needed
N       = Prob.N;
mLin    = Prob.mLin;
mNonLin = Prob.mNonLin;
if isempty(mLin)
   mLin = size(Prob.A,1);
end
if isempty(mNonLin)
   mNonLin = max(length(Prob.c_L),length(Prob.c_U));
end

M = mLin + mNonLin;

if isfield(Prob,'optParam')
   if isempty(Prob.optParam)
      Prob.optParam=optParamDef(Solver,probType,N,N,M);
   else
      Prob.optParam=optParamSet(Prob.optParam,Solver,probType,N,N,M);
   end
else
   Prob.optParam=optParamDef(Solver,probType,N,N,M);
end
if isfield(Prob,'LineParam')
   Prob.LineParam=LineParamSet(Prob.LineParam);
else
   Prob.LineParam=LineParamDef;
end

if isfield(Prob,'PriLevOpt')
   PriLevOpt=Prob.PriLevOpt;
else
   PriLevOpt      = 0;
   Prob.PriLevOpt = PriLevOpt;
end

% Is this really needed any longer?
if (PriLevOpt > 1) | (PriLevOpt >= 1 & ...
   (strcmpi(Solver,'constr') | strcmpi(Solver,'MINOS')))
   fprintf('\n\n==>==>==>==>==>==>==>==>==>==>==>==>==>==>==>==>==>==>\n\n');
   fprintf('Problem %4.0f: %s',probNumber,Prob.Name);
   fprintf('\n\n');
end

% PRESOLVE
%if Prob.optParam.PreSolve
%   Prob = preSolve(Prob);
%end
end

Result   = [];
NLLSVars = 0;
GlobVars = 0;

[TomV,os,TV]=tomlabVersion;

switch lower(Solver)
 case 'nlpsolve' % SQP. Fletcher-Leyffer

   Result = nlpSolve(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
   NLLSVars = 1;
   GlobVars = 1;

%case 'nlpsolv' % SQP. Fletcher-Leyffer, version II
%
%  Result = nlpSolv(Prob);
%
%  Result.FuncEv=n_f;
%  Result.GradEv=n_g;
%  Result.ConstrEv=n_c;

 case 'consolve' % Schittkowski Augmented Lagrangian SQP
   % alg=Prob.Solver.Alg. Default == 0
   % alg=0 Schittkowski SQP
   % alg=1 Han-Powell SQP. Exact L_1 penalty

   Result = conSolve(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
   NLLSVars = 1;
   GlobVars = 1;

 case 'strustr' % Structured Trust Region
   % sTrustr runs a structured trust region algorithm 
   % (Conn/Gould/Sartenaer/Toint)

   Result = sTrustr(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;

 case 'clssolve' %General TOMLAB constrained LS, clsSolve:
   % alg=Prob.Solver.Alg. Default == 1
   %alg==0 Gauss Newton with Subspace Minimization
   %alg==1 Fletcher-Xu Hybrid Method. Modified Gauss-Newton - BFGS
   %alg==2 Al-Baali - Fletcher Hybrid Method. Modified GN - BFGS
   %alg==3 Huschens TSSM algorithm. SIAM J. Optimization 1994-1 108-129

   Result = clsSolve(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=0;

 case 'glbsolve'
   Result= glbSolve(Prob);

 case 'glbfast'
   Result= glbFast(Prob);

 case 'ego'
   Result= ego(Prob);

 case 'glcsolve'
   Result= glcSolve(Prob);

 case 'glcfast'
   Result= glcFast(Prob);

 case 'glccluster'
   Result= glcCluster(Prob);

 case 'minlpsolve' % General TOMLAB MINLP 

   Result = minlpSolve(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;

   NLLSVars = 1;
   GlobVars = 1;

 case 'rbfsolve'
   Result= rbfSolve(Prob);

 case 'strustr'    % Structured Trust Region

   Result = sTrustr(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;

 case 'lpsolve'         % TOMLAB general LP solver

   Result = lpSolve(Prob);

   if Result.ExitFlag~=0 & PriLev > 0
      if Result.ExitFlag==2
         fprintf('\nlpSolve: Unbounded solution. FLAG = %3.0f\n',...
                  Result.ExitFlag);
      else
         fprintf('\nlpSolve: ERROR FLAG = %3.0f\n',Result.ExitFlag);
      end
   end
 %case {'qpsolve','qld'}      % TOMLAB general QP solver
 case 'qpsolve'      % TOMLAB general QP solver

   Result = qpSolve(Prob);

   %Result.FuncEv=Result.Iter;
   %Result.GradEv=Result.Iter;
   %Result.ConstrEv=Result.Iter;
   if Result.ExitFlag~=0
      fprintf('\nqpSolve: ERROR FLAG = %3.0f\n',Result.ExitFlag);
      if Result.ExitFlag==2
         fprintf('\nThe problem is infeasible\n');
      elseif Result.ExitFlag==4
         fprintf('\nUnbounded solution\n');
      end
   end

 case 'mipsolve'     % TOMLAB general Branch and bound MIP solver

   if PriLev > 1
      fprintf('\n\nCall TOMLAB routine mipSolve\n');
   end

   Result = mipSolve(Prob);

   if Result.ExitFlag~=0 & PriLev > 0
      if Result.ExitFlag==2
         fprintf('\nmipSolve: No integer solution found. FLAG = %3.0f\n',...
                  Result.ExitFlag);
      else
         fprintf('\nmipSolve: ERROR FLAG = %3.0f\n',Result.ExitFlag);
      end
   end

 case 'cutplane'     % TOMLAB cutting plane solver

   if PriLev > 1
      fprintf('\n\nCall TOMLAB routine cutplane\n');
   end

   Result = cutplane(Prob);

   if Result.ExitFlag~=0 & PriLev > 0
      if Result.ExitFlag==2
         fprintf('\ncutplane: No integer solution found. FLAG = %3.0f\n',...
                  Result.ExitFlag);
      else
         fprintf('\ncutplane: ERROR FLAG = %3.0f\n',Result.ExitFlag);
      end
   end

 case 'ucsolve'     % TOMLAB unconstrained optimization solver
   Result = ucSolve(Prob);

   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=0;



% -----------------------
% SOL MEX file interfaces
% -----------------------
 case 'lpopt'      % Run mex interface to LPOPT
   if TV(2)
      Prob   = ProbCheck(Prob,'lpopt',8);
      Result = lpoptTL(Prob);
      if isempty(Result), EmptyResult('lpopt'), end
   else
      fprintf('No valid license for the LPOPT solver\n');
   end

 case 'qpopt'      % Run mex interface to QPOPT
   if TV(2)
      Prob   = ProbCheck(Prob,'qpopt',2);
      Result = qpoptTL(Prob);
      if isempty(Result), EmptyResult('qpopt'), end
   else
      fprintf('No valid license for the QPOPT solver\n');
   end

 case 'sqopt'      % Run mex interface to SQOPT
   if TV(4)
      Prob=ProbCheck(Prob,'sqopt',2);
      Result = sqoptTL(Prob);
      if isempty(Result), EmptyResult('SQOPT'), end
   else
      fprintf('No valid license for the SQOPT solver\n');
   end
   
 case 'sqopt7'      % Run mex interface to SQOPT 7
   if TV(4)
      Prob=ProbCheck(Prob,'sqopt',2);
      Result = sqopt7TL(Prob);
      if isempty(Result), EmptyResult('SQOPT'), end
   else
      fprintf('No valid license for the SQOPT solver\n');
   end
   
 case 'lssol'      % Run mex interface to LSSOL
   if TV(3)
      Prob=ProbCheck(Prob,'lssol',5);
      Result = lssolTL(Prob);
      if isempty(Result), EmptyResult('LSSOL'), end
   else
      fprintf('No valid license for the LSSOL solver\n');
   end

 case 'nlssol'     % Run mex interface to NLSSOL
   if TV(3)
      Prob=ProbCheck(Prob,'nlssol',6);
      Result = nlssolTL(Prob);
      if isempty(Result), EmptyResult('NLSSOL'), end
   else
      fprintf('No valid license for the NLSSOL solver\n');
   end

 case 'minos'
   if TV(2)
      Prob   = ProbCheck(Prob,'minos',3);
      Result = minosTL(Prob);
      if isempty(Result), EmptyResult('minos'), end
   else
      fprintf('No valid license for the MINOS solver\n');
   end
   NLLSVars = 1;
 case 'lp-minos'
   if TV(2)
      Prob   = ProbCheck(Prob,'minos',8);
      Result = minoslpTL(Prob);
      if isempty(Result), EmptyResult('minos'), end
   else
      fprintf('No valid license for the MINOS solver\n');
   end

 case 'qp-minos'
   if TV(2)
      Prob   = ProbCheck(Prob,'minos',2);
      Result = minosqpTL(Prob);
      if isempty(Result), EmptyResult('minos'), end
   else
      fprintf('No valid license for the MINOS solver\n');
   end

 case 'npsol'      % Run mex interface to NPSOL
   if TV(3)
      Prob=ProbCheck(Prob,'npsol',3);
      Result = npsolTL(Prob);
      if isempty(Result), EmptyResult('NPSOL'), end
   else
      fprintf('No valid license for the NPSOL solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case 'snopt'      % Run mex interface to SNOPT
   if TV(4)
      Prob=ProbCheck(Prob,'snopt',3);
      Result = snoptTL(Prob);
      if isempty(Result), EmptyResult('SNOPT'), end
   else
      fprintf('No valid license for the SNOPT solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case 'snopt7'      % Experimental SNOPT 7 interface
    if TV(4)
       Prob=ProbCheck(Prob,'snopt',3);
       Result = snopt7TL(Prob);
       if isempty(Result), EmptyResult('SNOPT7'), end
    else
      fprintf('No valid license for the SNOPT 7 solver\n');
    end
    NLLSVars = 1;
    GlobVars = 1;
% -----------------------
% Other MEX file interfaces
% -----------------------
 case 'tlsqr'      % Run mex interface to LSQR linear least squares
   Prob   = ProbCheck(Prob,'Tlsqr',5);
   Result = TlsqrTL(Prob);
   if isempty(Result), EmptyResult('Tlsqr'), end

 case 'qld'
   Prob   = ProbCheck(Prob,'qld',2);
   Result = qldTL(Prob);
   if isempty(Result), EmptyResult('qld'), end
 case 'lsei'      % Run mex interface to LSEI linear least squares
   Prob   = ProbCheck(Prob,'lsei',5);
   Result = lseiTL(Prob);
   if isempty(Result), EmptyResult('lsei'), end

 case 'bqpd'
   if TV(7)
      Prob   = ProbCheck(Prob,'bqpd',2);
      Result = bqpdTL(Prob);
      if isempty(Result), EmptyResult('bqpd'), end
   else
      fprintf('No valid license for the BQPD solver\n');
   end
 case 'miqpbb'
   if TV(7)
      Prob   = ProbCheck(Prob,'miqpbb',11);
      Result = miqpBBTL(Prob);
      if isempty(Result), EmptyResult('miqpBB'), end
   else
      fprintf('No valid license for the miqpBB solver\n');
   end
 case {'filtersqp', 'filsqp'}
   if TV(7)
      Prob   = ProbCheck(Prob,'filterSQP',3);
      Result = filterSQPTL(Prob);
      if isempty(Result), EmptyResult('filterSQP'), end
   else
      fprintf('No valid license for the filterSQP solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case 'minlpbb'
   if TV(7)
      Prob   = ProbCheck(Prob,'minlpBB',12);
      Result = minlpBBTL(Prob);
      if isempty(Result), EmptyResult('minlpBB'), end
   else
      fprintf('No valid license for the minlpBB solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;

 case 'pensdp'
   if TV(6)
      Prob   = ProbCheck(Prob,'PENSDP',13);
      Result = pensdpTL(Prob);
      if isempty(Result), EmptyResult('PENSDP'), end
   else
      fprintf('No valid license for the PENSDP solver\n');
   end

 case 'penbmi'
   if TV(10)
      Prob   = ProbCheck(Prob,'PENBMI',14);
      Result = penbmiTL(Prob);
      if isempty(Result), EmptyResult('PENBMI'), end
   else
      fprintf('No valid license for the PENBMI solver\n');
   end
 case {'nlpql', 'nlpqlp'}
   if TV(17) 
      Prob   = ProbCheck(Prob,'nlpql',3);
      Result = nlpqlTL(Prob);
      if isempty(Result), EmptyResult('NLPQL'), end
   else
      fprintf('No valid license for the NLPQL solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;

case {'dfnlp', 'dfnlpd'}
   if TV(17) 
      Prob   = ProbCheck(Prob,'dfnlp',3);
      Result = dfnlpTL(Prob);
      if isempty(Result), EmptyResult('DFNLP'), end
   else
      fprintf('No valid license for the DFNLP solver\n');
   end

case {'nlpjob'}
   if TV(17) 
      Prob   = ProbCheck(Prob,'nlpjob',3);
      Result = nlpjobTL(Prob);
      if isempty(Result), EmptyResult('NLPJOB'), end
   else
      fprintf('No valid license for the NLPJOB solver\n');
   end
   
 case 'pdco'
   Prob   = ProbCheck(Prob,'pdco',3);
   Result = pdcoTL(Prob);

   NLLSVars = 1;
   GlobVars = 1;
 case 'pdsco'
   Prob   = ProbCheck(Prob,'pdsco',3);
   Result = pdscoTL(Prob);

   NLLSVars = 1;
   GlobVars = 1;
 
 case 'oqnlp'
   if TV(14)
      Prob    = ProbCheck(Prob,'oqnlp',12);
      Result  = oqnlpTL(Prob);
   else
      fprintf('No valid license for the OQNLP solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
   
 case 'msnlp'
   if TV(22)
      Prob    = ProbCheck(Prob,'msnlp',12);
      Result  = msnlpTL(Prob);
   else
      fprintf('No valid license for the MSNLP solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;   

 case 'lsgrg2'
   if TV(22)
      Prob    = ProbCheck(Prob,'lsgrg2',3);
      Result  = lsgrg2TL(Prob);
   else
      fprintf('No valid license for the LSGRG2 solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;   

case {'knitro','Tknitro'}
   if TV(11) 
      Prob    = ProbCheck(Prob,'knitro',3);
      Result  = knitroTL(Prob);
      if isempty(Result), EmptyResult('KNITRO'), end
   else
      fprintf('No valid license for the KNITRO solver\n');
   end  
   NLLSVars = 1;
   GlobVars = 1;
    
case {'conopt','Tconopt'}
   if TV(12)
      Prob = ProbCheck(Prob,'conopt',3);
      Result = conoptTL(Prob);
      if isempty(Result), EmptyResult('CONOPT'), end
   else
      fprintf('No valid license for the CONOPT solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
    
% -------------------------
% BOEING Solvers
% -------------------------

case {'barnlp'}
   %if TV(12)
      Prob = ProbCheck(Prob,'barnlp',3);
      Result = barnlpTL(Prob);
      if isempty(Result), EmptyResult('BARNLP'), end
   %else
   %   fprintf('No valid license for the BARNLP solver\n');
   %end
   NLLSVars = 1;
   
case {'sprnlp'}
   %if TV(12)
      Prob = ProbCheck(Prob,'sprnlp',3);
      Result = sprnlpTL(Prob);
      if isempty(Result), EmptyResult('SNPNLP'), end
   %else
   %   fprintf('No valid license for the SNPNLP solver\n');
   %end
   NLLSVars = 1;

case {'barqp'}
   %if TV(12)
   Prob = ProbCheck(Prob,'barnlp',2);
   Result = barqpTL(Prob);
   if isempty(Result), EmptyResult('BARQP'), end
   %else
   %   fprintf('No valid license for the BARQP solver\n');
   %end
   NLLSVars = 1;
   
% -----------------------
% Xpress-MP MEX file interface
% -----------------------
 case {'xpress-mp','xpress','xpressmp'}   % Run mex interface to Xpress-MP
   if TV(8)
      Prob=ProbCheck(Prob,'xpress-mp',11);
      Result = xpressTL(Prob);
      if isempty(Result), EmptyResult('Xpress-MP'), end
   else
      fprintf('No valid license for the Xpress-MP solver\n');
   end
% -----------------------
% CPLEX MEX file interface
% -----------------------
 case 'cplex'      % Run mex interface to CPLEX
   if TV(9)
      Prob=ProbCheck(Prob,'CPLEX',11);
      Result = cplexTL(Prob);
      if isempty(Result), EmptyResult('CPLEX'), end
   else
      fprintf('No valid license for the CPLEX solver\n');
   end

   % ------------------------
   % XA MEX file interface
   % ------------------------
case 'xa'
   %if TV(16)
      Prob = ProbCheck(Prob,'xa',2);
      Result = xaTL(Prob);
      if isempty(Result), EmptyResult('XA'), end
   %else
   %   fprintf('No valid license for the XA solver\n');
   %end
  
% -----------------------
% Opt tbx 1.x interfaces
% -----------------------

 case 'constr'
   Result=opt15Run('CONSTR',Prob);
   NLLSVars = 1;
   GlobVars = 1;

 case 'fmins'
   Result=opt15Run('FMINS',Prob);

 case 'leastsq'
   Result=opt15Run('LEASTSQ',Prob);

 case 'lp'
   Result=opt15Run('LP',Prob);

 case 'qp'
   Result=opt15Run('QP',Prob);

 case 'fminu'
   Result=opt15Run('FMINU',Prob);

% -----------------------
% Opt tbx 2.x interfaces
% -----------------------
 case 'fminunc'
   Result=opt20Run('FMINUNC',Prob);

 case 'fmincon'
   Result=opt20Run('FMINCON',Prob);
   NLLSVars = 1;
   GlobVars = 1;

 case 'fminsearch'
   Result=opt20Run('FMINSEARCH',Prob);

 case 'lsqnonlin'
   Result=opt20Run('LSQNONLIN',Prob);

 case 'lsqlin'
   Result=opt20Run('LSQLIN',Prob);

 case 'linprog'
   Result=opt20Run('LINPROG',Prob);

 case 'quadprog'
   Result=opt20Run('QUADPROG',Prob);
   
   % ---------------------
   % LGO interface
   % ---------------------

case 'lgo'
   Prob=ProbCheck(Prob,'lgo',12);
   Result=lgoTL(Prob);
   if isempty(Result), EmptyResult('LGO'); end
   
   % ---------------------
   % PATH interface
   % ---------------------
   
case 'path'
    if Prob.probType == 17
        Prob = ProbCheck(Prob,'PATH',17);
    elseif Prob.probType == 8
        Prob = ProbCheck(Prob,'PATH',8);
    elseif Prob.probType == 2
        Prob = ProbCheck(Prob,'PATH',2);
    else
        Prob = ProbCheck(Prob,'PATH',18);
    end    
    Result=pathTL(Prob);
    if isempty(Result), EmptyResult('PATH'); end
   
 otherwise
   fprintf('Solver %s',Solver);
   fprintf(' NOT found\n');
   error('Illegal solver algorithm!')
end

if any(Prob.probType==[4 5 6])
   if ~isempty(Result)
      Result.ResEv=n_r;
      Result.JacEv=n_J;
   end
end

if isfield(Prob,'AMPL')
    if ~isempty(Prob.AMPL)
        Result = postSolveAMPL(Result);
    end
end

% HKH CHECK THIS CAREFULLY, only used in mex routines before
if NLLSVars
   if any(Prob.probType==[4 6 11])
      %global n_r n_J n_f n_g n_H n_c 
      Result.ResEv=n_r;
      Result.JacEv=n_J;
   end
end
if GlobVars
   if any(Prob.probType==[1 2 3 10])
      %global n_f n_g n_H n_c 
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
      Result.ConstrEv=n_c;
   elseif any(Prob.probType==[1 9])
      %global n_f n_g n_H 
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
   elseif any(Prob.probType==[4 6 11])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
      Result.ConstrEv=n_c;
   end
end

PrintResult(Result,PriLev);


function EmptyResult(Code)
fprintf('\n\n')
fprintf('tomRun: ');
fprintf('ERROR! %s MEX-file interface returned empty Result',Code);
fprintf('\n\n')
fprintf('Is the %s MEX file missing?\n',Code);
which ('-all',Code)
fprintf('\n\n')

% MODIFICATION LOG:
%
% 981013  hkh  Added call to iniSolve and endSolve. Now call MEX with mexRun
%              Deleted 2nd output parameter Prob. Instead in Result.Prob.
%              PrintResult now function, called with Result and Pri.
% 981022  hkh  Delete global p_f, p_g etc.
% 981026  hkh  Redefine input to ???Run files.
% 981102  hkh  Set an initial value on Prob.P. Matlab is giving warnings
%              when setting empty, so set the more dangerous value of -1.
% 981108  hkh  Changed PriLev levels from 0-3
% 981110  hkh  Improve comments, discussing default values
%              Use Prob.probFile, if set in Prob and empty as argument
%              Change to use globals MAX_x, MAX_c, MAX_r
% 981118  hkh  Add call to iniSolve to clear globals.
% 981129  hkh  Add global variable xGUI
%              Write messages to GUI window if xGUI true. 
% 990629  hkh  Made it a function
% 000726  hkh  Add snopt solver
% 000820  hkh  Use mexSOL for SOL solvers, change call to tomRun
% 000820  hkh  Update for v3.0
% 000925  hkh  Revision, remove optParam
% 001106  hkh  General use of tomRun, and no other driver routines
% 010715  hkh  Adding glbFast
% 010815  hkh  Adding glcFast
% 011111  hkh  Adding glcCluster and rbfSolve
% 011204  hkh  Adding lsqlin
% 020512  hkh  Adding Dundee QP solver bqpd
% 020621  hkh  Adding Dundee MIQP solver MIQPbb
% 020630  hkh  Adding Dundee solvers MINLPbb, filterSQP; and PENSDP
% 020702  hkh  Adding CPLEX
% 030110  hkh  Do preSolve if Prob.optParam.PreSolve
% 030116  hkh  Change lsqr to Tlsqr
% 030123  hkh  Add pdco and pdsco
% 030129  hkh  Remove empty varargin in call to opt20Run and opt15Run
% 030129  ango Edit names, Dundee solvers
% 030206  hkh  Change filSQP to filterSQP in mexRun call
% 030317  ango Add OQNLP
% 030427  hkh  Add KS solver nlpql
% 030514  ango Add KNITRO
% 030603  ango Add CONOPT
% 030613  ango Add dfnlp
% 030709  ango Add nlpjob
% 030728  medv Added postSolve for AMPL problems.
% 031023  ango BARNLP, SPRNLP added
% 031029  ango XA added
% 031103  ango Synonyms for xpress-mp added: xpress, xpressmp
% 040101  hkh  NPOPT deleted
% 040102  hkh  Avoid tests, do them in ProbCheck
% 040101  hkh  Make EmptyResult routine, incl. missing dll check
% 040109  ango Add LGO
% 040111  hkh  Removed calls to iniSolve and endSolve, done in solvers and TL
% 040120  ango SNOPT 7 added
% 040419  med  Added PATH
% 041023  hkh  Added minlpSolve
% 041123  hkh  Default PriLev 2 for probFile input, otherwise 0 
% 041123  hkh  Revision to change order of PriLev and ask
% 041216  med  OQNLP and MSNLP added
% 050124  ango SQOPT7 added
