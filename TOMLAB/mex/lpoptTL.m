% TOMLAB LPOPT LP Solver
%
% function Result = lpoptTL(Prob)
%
% INPUT:  
% Prob   Problem structure in TOMLAB format.
%
% Fields used in input structure Prob (call Prob=lpAssign(...) or 
% Prob=ProbDef; to define Prob).
%
% x_L, x_U  Bounds on variables. 
% b_L, b_U  Bounds on linear constraints. 
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% PriLevOpt Print Level.
% WarmStart If true, use warm start, otherwise cold start.
%
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution (and slacks not used) from previous run (Warm start).
% iState    State for the constraints from previous run (Warm start).
% SpecsFile Name of user defined SPECS file, read BEFORE optPar() is used.
% PrintFile Name of SOL Print file. Amount and type of printing determined 
%           by SPECS parameters or optPar parameters.
% SummFile  Name of SOL Summary File.
% optPar    Elements > -999 takes precedence over corresponding TOMLAB
%           params.
%
% -----------------------------------------------
% How optPar is used for setting SPECS parameters:
% -----------------------------------------------
%
% optPar  Structure with optimization parameters. 
% Fields not defined in optParam is defined by a call to SOLSet.m
%
% LPOPT keywords in optPar(#) - and the fields used in optParam
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
% Printing
% 1.  PRINT LEVEL                   0        10                {0,1,5,10,20,30}
% 3.  PRINT FILE                    0        0                 Fortran Unit #
%           SET BY INTERFACE IF PrintFile is given
% 4.  SUMMARY FILE                  0        0                 Fortran Unit #
%           SET BY INTERFACE IF SummFile  is given
% 10. OPTIMALITY TOLERANCE          >0       1.1E-8            sqrt(eps)
%     -optParam.eps_x; 
% 11. FEASIBILITY TOLERANCE         >0       1.1E-8            sqrt(eps)
%     -optParam.bTol; 
% 21. CRASH TOLERANCE               >0       0.01      <1      
% 27. RANK TOLERANCE                >0       1.1E-14           100*eps
%     -optParam.eps_Rank;
% 30. ITERATION LIMIT               >0       max(50,5(n+m))    
%     -optPar(30)=optParam.MaxIter;
% 33. MIN SUM YES (or NO)           0        0         1       1=min infeas
%           IF 1 (MIN SUM YES), minimize the infeasibilities before return
% 36. FEASIBILITY PHASE ITERATIONS  >0       max(50,5(n+m))    
%     -optPar(36)=optParam.MaxMinorIter;
% 45. INFINITE STEP SIZE            >0       1E20
% 51. CHECK FREQUENCY               >0       50
% 52. EXPAND FREQUENCY              >0       5
%
% OUTPUT: 
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Exact gradient (trivially c).
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from lpopt.m (similar to TOMLAB).
% Inform   LPOPT information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (LPOPT).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field
% xs       Solution and slack variables.
% iState   State for variables and constraints in iState.
%
% -----------------------------------------------------------------------
%
% For a problem description, see lpAssign.m
%
% -------------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written July 25, 2000.  Last modified Dec 22, 2004.
%

function Result = lpoptTL(Prob)


if nargin < 1, error('lpoptTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

Prob.solvType = 1; % LP solver

Prob = iniSolve(Prob,1,0,0);

Result=ResultDef(Prob);
Result.Solver='LPOPT';
Result.SolverAlgorithm='LPOPT 1.0-10 LP code';

%
% Define lower and upper bound arrays for LPOPT
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%

BIG=1E20;

[bl, bu, n, m, m2] = defblbu(Prob, BIG, 1);

PriLev=Prob.PriLevOpt;

% Initial checks on the inputs

nb = n+m;

% Check on Warm start

if Prob.WarmStart
   % Warm start for SOL solver
   Warm   = 1;
   x_0    = Prob.SOL.xs(1:n);
   if isempty(Prob.SOL.iState)
      % Use hs field
      iState = Prob.SOL.hs(1:nb);
   else
      iState=Prob.SOL.iState
   end
   % NOT USED cLamda=Prob.SOL.cLamda;
else
   Warm   = 0;
   x_0    = Prob.x_0(:);
   if any(isnan(x_0)),x_0 = []; end
   iState = [];
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
end

[mA,nA] = size(Prob.A);

% Check if any linear part
c = Prob.QP.c(:);

% Determine type of problem
if isempty(c) | all(c==0)
   Result.f_0=0;
else
   Result.f_0=c(1:n)'*x_0(1:n);
end

% Set up the constraint matrix A 

if m==0     % isempty(Prob.A)
   % Construct one linear constraint from some bounded variable when m==0
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);

   bl=[bl;max(bl(i),-BIG)];
   bu=[bu;min(bu(i), BIG)];
   Prob.A = sparse(zeros(1,n));
   Prob.A(1,i)=1;
   if length(iState) == n
      iState=[iState;0];
   end
end

if PriLev > 2
   fprintf('# of non zero entries in linear constraints %d. ',nnz(Prob.A))
   if PriLev > 5
      fprintf('Number of variables          %d.\n',n)
      fprintf('Number of linear constraints %d.\n',m)
   end
   if PriLev >= 1000
      xprinte(bl(1:n),'x_L: ');
      xprinte(bu(1:n),'x_U: ');
      xprinte(bl(n+1:nb),'b_L: ');
      xprinte(bu(n+1:nb),'b_U: ');
      xprinte(c,'c: ');
      if PriLev > 0, pause; end
   end
   if PriLev >= 7
      if PriLev >= 10, pause; end
      printmat(full(Prob.A),'Constraint matrix A:')
      if PriLev >= 10, pause; end
      
      disp('x_L x_0 x_U');
      mPrint([bl(1:n) x_0 bu(1:n)],'xl/x0/xu:')
      mPrint([bl(n+1:n+m) bu(n+1:n+m)],'bl/bu:')
      if PriLev >= 10, pause; end
   end
end

[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('lpopt',8,...
         0, 0, size(Prob.A,1), Prob);

[Inform, Iter, iState, Ax, cLamda, Obj, x_k] = lpopt( ...
         full(Prob.A), bl, bu, c, Warm, x_0, iState, ...
         SpecsFile, PrintFile, SummFile, PriLev, optPar );

switch Inform
   case 4
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 3
     ExitFlag=4;  % Infeasible
   case {1}
     ExitFlag=3;  % ??? May be optimal
   case {5,6,7}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = cLamda;
Result.g_k = c;

% Warm start
% The dummy objective row is last in xs, value -Obj

% Could use Result.x_k instead
Result.SOL.xs=[x_k;zeros(m,1);-Obj];

Result.SOL.iState=iState(1:nb);
% Put iState also in hs field
Result.SOL.hs=[iState;0];

% Could use Result.v_k instead
%Result.SOL.cLamda=cLamda;

Result.SOL.optPar=optPar;

Result.FuncEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

% Compute Result.xState and Result.QP.B only
Result = StateDef(Result, x_k, [], [], Prob.optParam.xTol, [], [], bl, bu, 1); 

% State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%Result.xState=iState(1:n);
% Use already computed iState
Result.bState=iState(n+1:n+m);

switch Inform
   case 0
      Text = 'Optimal solution with unique minimizer found';
   case 1
      Text = 'A dead point was reached';
   case 2
      Text = 'The solution appears to be unbounded (or badly scaled)';
   case 3
      Text = str2mat('The constraints could not be satisfied.' ...
                    ,'The problem has no feasible solution');
   case 4
      Text = 'Too many iterations, in either phase';
   case 5
      Text = str2mat('The Maximum degrees of freedom is too small' ...
             ,'The reduced Hessian must expand if further progress is' ...
             ,'too be made');
   case 6
      Text = 'An input parameter was invalid';
   case 7
      Text = 'The problem type was not recognized';
   otherwise
      Text = 'NOTE: UNKNOWN LPOPT Inform value.';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLPOPT solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LPOPT: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');
   
   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Iterations      %7d. ',Iter);
   fprintf('\n');
   

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector iState for x and constraints = \n');
      xprinti(iState(1:min(length(iState),MAX_x)),'iState: ');
   end

   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (cLamda) = \n');
      xprinte(cLamda(1:min(length(cLamda),MAX_c)),'cLamda:');
   end
end
Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000716 hkh New lpoptMex interface.
% 000916 hkh Return string ExitText with interpretation of Inform flag 
% 001004 hkh Revision for direct call to the mex file
% 010710 hkh Error in comments
% 020417 hkh Send optPar back in Result.SOL
% 020821 hkh Remove code appearing twice
% 040102 hkh Revision for v4.2, call iniSolve and endSolve
% 040102 hkh Return only Iter ~= 0, FuncEv=ConstrEv=0
% 040602 med Help fix PriLev to PriLevOpt
% 041202 hkh Revise calls to defblbu and StateDef
% 041222 med Safeguard added for x_0