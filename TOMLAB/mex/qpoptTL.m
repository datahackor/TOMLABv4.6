% TOMLAB QPOPT QP/LP Solver
%
% function Result = qpoptTL(Prob)
%
% INPUT:  
% Prob   Problem structure in TOMLAB format
%
% Fields used in input structure Prob (call Prob=qpAssign(...) or 
% Prob=ProbDef; to define Prob)
%
% x_L, x_U  Bounds on variables. 
% b_L, b_U  Bounds on linear constraints. 
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% QP.F      Quadratic matrix of size nnObj x nnObj. nnObj < n is OK.
% PriLevOpt Print Level.
% WarmStart If true, use warm start, otherwise cold start.
%
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution (and slacks not used) from previous run (Warm start).
% iState    State for the constraints from previous run (Warm start).
% cLamda    Lagrangian multipliers from previous run (Warm start).
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
% QPOPT keywords in optPar(#) - and the fields used in optParam
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
%     -optParam.MaxIter;
% 33. MIN SUM YES (or NO)           0        0         1       1=min infeas
%           IF 1 (MIN SUM YES), minimize the infeasibilities before return
% 36. FEASIBILITY PHASE ITERATIONS  >0       max(50,5(n+m))    
%     -optPar(36)=optParam.MinorIter;
% 45. INFINITE STEP SIZE            >0       1E20
% 47. HESSIAN ROWS                  0        n         n       0 if FP or LP
%           IMPLICITLY GIVEN BY THE DIMENSIONS OF H IN THE CALL FROM MATLAB
% 48. MAX DEGREES OF FREEDOM        0        n         n       
%           ONLY USED IF HESSIAN ROWS == N
% 51. CHECK FREQUENCY               >0       50
% 52. EXPAND FREQUENCY              >0       5
%
%
% OUTPUT: 
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector)
% ExitFlag Exit status from qpopt.m (similar to TOMLAB).
% Inform   QPOPT information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (QPOPT).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field:
% xs       Solution and slack variables.
% cLamda   Lagrangian multipliers.
% iState   State for variables and constraints in iState.
%
% -----------------------------------------------------------------------
%
% For a problem description, see qpAssign.m
%
% -------------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written July 16, 2000.  Last modified Dec 22, 2004.
%

function Result = qpoptTL(Prob)

if nargin < 1, error('qpoptTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

Prob.solvType = 2; % QP solver

Prob = iniSolve(Prob,2,0,0);

Result=ResultDef(Prob);
Result.Solver='QPOPT';
Result.SolverAlgorithm='QPOPT 1.0 QP/LP code';

PriLev=Prob.PriLevOpt;

%
% Define lower and upper bound arrays for QPOPT
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
   % Initial values
   x_0     = Prob.x_0(:);
   if any(isnan(x_0)),x_0 = []; end
   Warm    = 0;
   iState  = [];
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
end

H     = Prob.QP.F;
nnObj = size(H,1); % number of nonlinear variables

% Check if any linear part
c = Prob.QP.c(:);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m, error('Linear constraints A MUST have m rows!'); end 
end 


% Determine type of problem
if isempty(c) | all(c==0)
   if isempty(H) | all(H == 0)
      Result.f_0=0;
   else
      Result.f_0=0.5*x_0(1:nnObj)'*H*x_0(1:nnObj);
   end
else
   if isempty(H) | all(H == 0)
      Result.f_0=c(1:n)'*x_0(1:n);
   else
      Result.f_0=0.5*x_0(1:nnObj)'*H*x_0(1:nnObj) + c(1:n)'*x_0(1:n);
   end
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

[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('qpopt',2,...
         n, 0, size(Prob.A,1), Prob);

%optPar(10)=1E-6; % Opt
%optPar(11)=1E-5; % Feas
%optPar(52)=10000;
%optPar(1) = 10;
%PrintFile = 'qpopt.txt';

[Inform, Iter, iState, Ax, cLamda, Obj, x_k] = qpopt ( ...
         H, full(Prob.A), bl, bu, c, Warm, x_0, iState, ...
         SpecsFile, PrintFile, SummFile, PriLev, optPar );

switch Inform
   case 4
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 3
     ExitFlag=4;  % Infeasible
%   case {1}
%     ExitFlag=3;  % Rank problem
   case {5,6,7}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = cLamda;
% Change sign when iState == 2
ix             = find(iState==2);
Result.v_k(ix) = -cLamda(ix);

% Warm start
% The dummy objective row is last in xs, value -Obj

% Could use Result.x_k instead
Result.SOL.xs=[x_k;zeros(m,1);-Obj];

Result.SOL.iState=iState(1:nb);
% Put iState also in hs field
Result.SOL.hs=[iState;0];

% Could use Result.v_k instead
Result.SOL.cLamda=cLamda;
%format long
%[cLamda iState]

Result.SOL.optPar=optPar;

if ~isempty(c)
   if ~isempty(H)
      Result.g_k=H*x_k+c;
   else
      Result.g_k=c;
   end
elseif isempty(H)
   Result.g_k=[];
else
   Result.g_k=H*x_k;
end

Result.FuncEv    = 0;
Result.GradEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
Result.MinorIter = 0;
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
      Text = 'NOTE: UNKNOWN QPOPT Inform value.';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nQPOPT solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('QPOPT: Inform = %2d, ',Inform)
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
elseif Prob.optParam.IterPrint 
   fprintf('QPOPT: Inform =%2d, ',Inform)
end
Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000716 hkh New qpoptMex interface.
% 000916 hkh Return string ExitText with interpretation of Inform flag 
% 001004 hkh Revision for direct call to the mex file
% 010710 hkh Error in comments
% 020417 hkh Send optPar back in Result.SOL
% 020821 hkh Remove code appearing twice
% 030406 hkh Print Inform if Prob.optParam.IterPrint = 1, and PriLev == 0
% 040102 hkh Revision for v4.2, call iniSolve and endSolve
% 040102 hkh Return only Iter ~= 0, FuncEv=GradEv=ConstrEv=0
% 040602 med Help fix PriLev to PriLevOpt
% 041126 hkh A dead point changed from ExitFlag = 3 to 0
% 041203 hkh Revise calls to defblbu and StateDef
% 041221 hkh Use n, not nnObj in call to SOLSet, nnObj 0 for LP
% 041222 med Safeguard added for x_0
