% TOMLAB MINOS LP Solver
%
% function Result = minoslpTL(Prob)
%
% INPUT:  
% Prob   Problem structure in TOMLAB format
%
% Fields used in input structure Prob. Call Prob=lpAssign(...) or 
% Prob=ProbDef; to define Prob.
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
% xs        Solution and slacks from previous run.
% hs        State for solution and slacks from previous run.
% nS        Number of superbasics from previous run (always 0 for LP).
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
%
% MINOS keywords in optPar(#) - and the fields used in optParam
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% Printing
% 1.  PRINT LEVEL                   0        0         1     0=brief 1=LU stats
% 3.  PRINT FILE                    0        0         99      Fortran Unit #
% 4.  SUMMARY FILE                  0        0         99      Fortran Unit #
% 5.  PRINT FREQUENCY               0        100
% 6.  SUMMARY FREQUENCY             0        100
% 7.  SOLUTION YES/NO               0        1         1       1 = YES; 0 = NO
% 8.  SUPPRESS PARAMETERS           0        0         1       1 = True
%
% Convergence Tolerances
% 10. OPTIMALITY TOLERANCE          >0       1E-6
%     -optParam.eps_x
% 11. FEASIBILITY TOLERANCE         >0       1E-6  
%     -optParam.bTol
%
% Scaling
% 18. SCALE OPTION                  0        2         2       0,1,2
% 19. SCALE TOLERANCE               >0       0.9       <1
% 20. SCALE PRINT                   0        0         1       1 = True
%
% Other Tolerances
% 21. CRASH TOLERANCE               0        0.1       <1
% 23. LU FACTOR TOLERANCE           1        100.0
% 24. LU UPDATE TOLERANCE           1        10.0
% 25  LU DENSITY TOLERANCE          >0       0.5                            
% 26. LU SINGULARITY TOLERANCE      >0       3.25E-11          eps^(0.67)
% 27. PIVOT TOLERANCE               >0       3.25E-11          eps^(0.67)
%     -optParam.eps_Rank
%
% LP parameters
% 28. CRASH OPTION                  0        3         3       {0,1,2,3}
% 29. WEIGHT ON LINEAR OBJECTIVE    0.0      0.0               during Phase 1 
% 30. ITERATION LIMIT               0        3m
%     -optParam.MaxIter
% 31. PARTIAL PRICE                 1        10
%
% Solve with tight or loose tols (applies to LP but probably not useful)
% 44. COMPLETION                    0        1         1     0=PARTIAL 1=FULL
%
% Frequencies
% 51. CHECK FREQUENCY               >0       60
% 52. EXPAND FREQUENCY              >0       10000
% 53. FACTORIZATION FREQUENCY       >0       100
% 54. SAVE FREQUENCY                >0       100
%
% BASIS files
% 55. OLD BASIS file                0        0
% 56. NEW BASIS file                0        0
% 57. BACKUP BASIS file             0        0
% 58. INSERT file                   0        0
% 59. PUNCH file                    0        0
% 60. LOAD file                     0        0
% 61. DUMP file                     0        0
% 62. SOLUTION file                 0        0
%
% 63. LU COMPLETE PIVOTING          0        0         2    2=rook, 1=complete,0=partial
%     or LU PARTIAL PIVOTING                      
%     or LU ROOK PIVOTING
%
% OUTPUT: 
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Gradient c (linear objective).
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from minos.m (similar to TOMLAB).
% Inform   MINOS information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Not set (default 0)
% GradEv   Number of gradient evaluations. Not set (default 0)
% ConstrEv Number of constraint evaluations. Not set (default 0)
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. Not set (default 0)
% Solver   Name of the solver (minos).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field
% xs       Solution and slack variables.
% hs       State for variables and slacks in xs.
% nS       # of superbasics.
% nInf     # of infeasibilities.
% sInf     Sum of infeasibilities.
%
% -----------------------------------------------------------------------
%
% For a problem description, see lpAssign.m
%
% MINOS creates a print file (print.dat) and a summary file (summary.dat).
%
% -------------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written July 10, 2000.  Last modified Dec 22, 2004.
%

function Result = minoslpTL(Prob)

%#function lp_f lp_g lp_H

if nargin < 1, error('minoslpTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

Prob.solvType = 3;
Prob = iniSolve(Prob,3,0,0);

Result=ResultDef(Prob);
Result.Solver='MINOS';
Result.SolverAlgorithm='MINOS 5.51 LP (NLP) code';

PriLev=Prob.PriLevOpt;

%
% Define lower and upper bound arrays for MINOS
%
% Inf are changed to BIG (=1E10), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%

BIG=1E10;

% m2 = 0 for LPs
[bl, bu, n, m, m2] = defblbu(Prob, BIG, 1);

% Initial checks on the inputs

ObjAdd  = 0;

if isfield(Prob.SOL,'moremem')
   moremem = Prob.SOL.moremem;
else
   moremem = [];
end
if isempty(moremem), moremem = n+m; end

% Check if any linear part
c  = Prob.QP.c;
m3 = length(c) > 0;

nb = n+m+m3;

% Initial checks on the inputs

% Check on Warm start, then set hs,xs

if Prob.WarmStart
   % Warm start for SOL solver
   Warm = 1;
   xs   = Prob.SOL.xs;
   hs   = Prob.SOL.hs;
   x_0  = xs(1:n);
else
   Warm = 0;
   hs   = Prob.SOL.hs;
   % Initial values
   x_0  = Prob.x_0(:);
   if any(isnan(x_0)),x_0 = []; end
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
   xs   = x_0;
end

if isempty(Prob.A)
   % Construct one linear constraint from some bounded variable when m==0
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);

   bl=[bl;-min(bu(i),BIG)];
   bu=[bu;-max(bl(i),-BIG)];
   Prob.A = sparse(zeros(1,n));
   Prob.A(1,i)=1;
   if length(hs) < n+1
      hs=[hs;0];
   end
end

if m3 > 0
   A = [sparse(Prob.A);c'];
   iObj=m+1;  % The constraint row with linear objective term c is put last
   % Add the bounds for the objective in bl and bu
   bl=[bl;-BIG];
   bu=[bu; BIG];
else
   A=sparse(Prob.A);
   iObj=0;
end

if PriLev > 2
   fprintf('# of non zero entries in linear constraints %d. ',nnz(Prob.A))
   if PriLev > 5
      fprintf('Total number of constraints %d.\n',m)
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

% Must change sign, and interchange the constraint bounds, 
% i.e. set correct bounds on slack variables

v=-bl(n+1:n+m);

bl(n+1:n+m)=-bu([n+1:n+m]);

bu(n+1:n+m)=v;

[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('minos',8,...
         0, 0, size(A,1), Prob);

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

if isempty(c)
   Result.f_0=0;
else
   Result.f_0=c(1:n)'*x_0(1:n);
end

% pi = []; H = []; nS = 0; Parameter #14, #15 and #1

[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount] = ...
     minos( [], A, bl, bu, 0, 0, 0, Prob, iObj, optPar, ...
            Warm, hs, xs, [], [], ...
            SpecsFile, PrintFile, SummFile, ...
            PriLev, ObjAdd, moremem, Prob.Name );

Iter=iwCount(2);

switch Inform
   case 3
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 1
     ExitFlag=4;  % Infeasible
   case {10,22}
     ExitFlag=3;  % Rank problem
   case {5,7,8,20,21,30,32,42,43,44}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k   = Obj;
Result.x_k   = xs(1:n);
Result.x_0   = x_0;
Result.v_k   = [rc(1:n);pi(1:m)];
Result.g_k   = c;
Result.y_k   = pi(1:m);


% Saved for a warm start?

% The dummy objective row is last in xs, value -Obj, last in hs is 3 (basic)
Result.SOL.xs    = xs;
Result.SOL.hs    = hs;
Result.SOL.nS    = nS;
Result.SOL.nInf  = nInf;
Result.SOL.sInf  = sInf;


%Result.FuncEv    = 0;
%Result.GradEv    = 0;
%Result.ConstrEv  = 0;
%Result.MinorIter = 0;

Result.Iter      = Iter;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

Result = StateDef(Result, xs(1:n), A(1:m,:)*xs(1:n), [], ...
                  Prob.optParam.xTol, Prob.optParam.bTol, [], ...
                  [bl(1:n);-bu(n+1:n+m)], [bu(1:n);-bl(n+1:n+m)],1); 

switch Inform
   case 0
      Text = 'Optimal solution found';
   case 1
      Text = 'The problem is infeasible';
   case 2
      Text = 'The problem is unbounded (or badly scaled)';
   case 3
      Text = 'Too many iterations';
   case 4
      Text = str2mat('Apparent stall.  The solution has not changed' ...
             ,'for a large number of iterations (e.g. 1000).');
   case 5
      Text = 'The Superbasics limit is too small';
   case 6
      Text = str2mat('User requested termination by returning' ...
             ,'mode <= -2 from funobj or funcon');
   case 7
      Text = 'funobj seems to be giving incorrect gradients';
   case 8
      Text = 'funcon seems to be giving incorrect gradients';
   case 9
      Text = 'The current point cannot be improved';
   case 10
      Text = str2mat('Numerical error in trying to satisfy the linear ' ...
             ,'constraints (or the linearized nonlinear constraints)' ...
             ,'The basis is very ill-conditioned.');
   case 11
      Text = 'Cannot find a superbasic to replace a basic variable';
   case 12
      Text = str2mat('Basis factorization requested twice in a row' ...
             ,'Like case Inform = 9. Possibly convergence?');
   case 13
      Text = str2mat('Near-optimal solution found' ...
             ,'Should probably be treated as inform = 9');
   case 20
      Text = 'Not enough storage for the basis factorization';
   case 21
      Text = 'Error in basis package';
   case 22
      Text = str2mat('The basis is singular after several attempts' ...
            ,'to factorize it (and add slacks where necessary)');
   case 30
      Text = str2mat('An OLD BASIS file had dimensions that did' ...
             ,'not match the current problem');
   case 32
      Text = 'System error. Wrong number of basic variables';
   case 40
      Text = 'Fatal errors in the MPS file';
   case 41
      Text = 'Not enough storage to read the MPS file';
   case 42
      Text = 'Not enough storage to solve the problem';
   otherwise
      Text = 'NOTE: UNKNOWN MINOS Inform value.';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nMINOS LP solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('MINOS: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');
   
   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('iterations      %7d. ',Iter);
   fprintf('\n');
   fprintf('nInf            %7d. ',nInf);
   fprintf('nS              %7d. ',nS);
   fprintf('sInf %14.7e\n',sInf);

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(xs(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector hs for x and constraints = \n');
      xprinti(hs(1:min(length(hs),MAX_x)),'hs: ');
   end
   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (pi) = \n');
      xprinte(pi(1:min(length(pi),MAX_c)),'pi:');
      fprintf('Reduced costs rc: Last %d elements should be -v_k\n',m);
      xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'rc: ',' %14.9f',5);
   end
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000710 hkh  New minosMex interface.
% 000828 hkh  If file name given, and print level in optPar(1) is 0, set to max
% 000916 hkh  Return string ExitText with interpretation of Inform flag 
% 010412 hkh  Clean up. Avoid 3 last output parameters in minos call.
% 010903 hkh  Add complete pivoting option. Allow cold start change of hs
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 031118 hkh  Must change BIG from 1E20 to 1E10, otherwise incorrect solution
% 031120 ango Help for optPar(63) changed, MINOS 5.51
% 040101 hkh  Calls to iniSolve and endSolve
% 040413 hkh  Make full call to iniSolve
% 040602 med  Help fix PriLev to PriLevOpt
% 041202 hkh  Revise call to defblbu and StateDef
% 041222 med  Safeguard added for x_0