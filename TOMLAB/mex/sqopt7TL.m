% TOMLAB SQOPT QP/LP Solver
%
% function Result = sqoptTL(Prob)
%
% INPUT:  
% Prob   Problem structure in TOMLAB format
%
% Fields used in input structure Prob (call Prob=ProbDef; to define Prob)
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
% callback  If 1, use a callback to Matlab to compute QP.F * x for different x.
%           Faster when F is very large and almost dense, avoiding
%           copying of F from Matlab to MEX.
% xs        Solution and slacks from previous run.
% hs        State for solution and slacks from previous run.
% nS        Number of superbasics from previous run.
% hElastic  Defines which variables are elastic in elastic mode. hElastic(j):
%           0 = variable j is non-elastic and cannot be infeasible.
%           1 = variable j can violate its lower bound.
%           2 = variable j can violate its upper bound.
%           3 = variable j can violate either its lower or upper bound.
% moremem   Add more memory if SQOPT stops with not enough storage message.
% SpecsFile Name of user defined SPECS file, read BEFORE optPar() is used.
% PrintFile Name of SOL Print file. Amount and type of printing determined 
%           by SPECS parameters or optPar parameters.
%           Output is written on file sqoptpri.txt, if not given.
% SummFile  Name of SOL Summary File.
%           Output is written on file sqoptsum.txt, if not given.
% To make sqopt to not open and not write anything to file:
%    Set SpecsFile and PrintFile empty.
%    Set optPar(2) = 0 AND optPar(3) = 0.
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
% SQOPT keywords in optPar(#) - and the fields used in optParam
% 2.  PRINT LEVEL                   0        0         10      0, 1 or 10
%     PriLev=1 to PRINT LEVEL 0, 2 to 1, 3 to 10
% 3.  PRINT FILE                    0        0                 Fortran Unit #
%     Set to 9 by default
% 4.  SUMMARY FILE                  0        0                 Fortran Unit #
%     Set to 8 by default
% 5.  PRINT FREQUENCY               0        1                 
% 6.  SUMMARY FREQUENCY             0        100
% 7.  SOLUTION YES/NO               0        1         1       1 = YES; 0 = NO
% 8.  SUPPRESS PARAMETERS           0        0         1       1 = True
%
% Convergence Tolerances
% 11. FEASIBILITY TOLERANCE         >0       1E-6  
%     bTol     Feasibility tol on linear constraints 
% 12. OPTIMALITY TOLERANCE          >0       1E-6  
%     eps_x    Optimality tolerance
%
% Scaling
% 18. SCALE OPTION                  0        1 or 2    2       2 if LP,1 if QP
% 19. SCALE TOLERANCE               >0       0.9       <1
% 20. SCALE PRINT                   0        0         1       1 = True
%
% Other Tolerances
% 21. CRASH TOLERANCE               0        0.1       <1
% 23. LU FACTOR TOLERANCE           1        100 or 5          100 if LP
% 24. LU UPDATE TOLERANCE           1        10  or 5          10  if LP
% 25  LU DENSITY TOLERANCE          >0       0.6                            
% 26. LU SINGULARITY TOLERANCE      >0       3.25E-11          eps^(0.67)
% 27. PIVOT TOLERANCE               >0       3.25E-11          eps^(0.67)
%     eps_Rank Rank tolerance                          
%
% QP subproblems
% 28. CRASH OPTION                  0        0         3       {0,1,2,3}
% 29. ELASTIC WEIGHT                0        1.0
% 30. ITERATION LIMIT               0        3m
%     MaxIter  Iterations (major its)       
% 31. PARTIAL PRICE                 0        10 or 1           10 for LP,1 QP
% 32. MAXIMIZE                      0        0         1       1=maximize
% 39. REDUCED HESSIAN               >0       min(500,1+nnObj)
%     May also set Superbasics Limit instead of Reduced Hessian
% 45. UNBOUNDED STEP SIZE           >0       1E20
% 48. SUPERBASICS LIMIT             >0       min(500,1+nnObj)
% 49. ELASTIC MODE                  0        1                 {0,1,2}
% 50. ELASTIC OBJECTIVE             0        2                 {0,1,2}
%
% Frequencies
% 51. CHECK FREQUENCY               >0       60
% 52. EXPAND FREQUENCY              >0       10000
% 53. FACTORIZATION FREQUENCY       >0       50
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
% 63. LU COMPLETE PIVOTING          0        0         1    1=complete,0=partial
%     or LU PARTIAL PIVOTING                      
%
% OUTPUT: 
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial  solution vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from sqopt.m (similar to TOMLAB).
% Inform   SQOPT information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% Solver   Name of the solver (sqopt).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field:
% xs       Solution and slack variables.
% hs       State for variables and slacks in xs.
% nS       # of superbasics.
% nInf     # of infeasibilities.
% sInf     Sum of infeasibilities.
%
% -----------------------------------------------------------------------
%
% For a problem description, see qpAssign.m
%
% SQOPT creates a print file (print.dat) and a summary file (summary.dat).
%
% -------------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written July 16, 2000.  Last modified Dec 22, 2004.
%

function Result = sqopt7TL(Prob)

if nargin < 1, error('sqoptTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

Prob.solvType = 2; % QP solver

Prob = iniSolve(Prob,2,0,0);

Result=ResultDef(Prob);
Result.Solver='SQOPT';
Result.SolverAlgorithm='SQOPT 7.1-1 QP code';

PriLev=Prob.PriLevOpt;

%
% Define lower and upper bound arrays for SQOPT
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
iObj    = 0;
ObjAdd  = 0;

if isfield(Prob.SOL,'moremem')
   moremem = Prob.SOL.moremem;
else
   moremem = [];
end
if isfield(Prob.SOL,'callback')
   callback = Prob.SOL.callback; 
else
   callback = 0;
end

nb      = n+m;

% Check on Warm start, then set hs,xs,nS

if Prob.WarmStart
   % Warm start for SOL solver
   Warm = 1;
   xs   = Prob.SOL.xs(1:nb);
   hs   = Prob.SOL.hs(1:nb);
   nS   = Prob.SOL.nS;
   x_0  = xs(1:n);
else
   % Initial values
   x_0  = Prob.x_0(:);
   if any(isnan(x_0)),x_0 = []; end
   Warm = 0;
   nS = DefPar(Prob.SOL,'nS',0);
   hs = DefPar(Prob.SOL,'nS',[]);
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0  = max(bl(1:n),min(x_0,bu(1:n)));
   xs   = x_0;
end
% Check on Warm start, then set hs,xs,nS

hElast = DefPar(Prob.SOL,'hElastic',[]);

nnObj  = size(Prob.QP.F,1); % number of nonlinear variables
c      = Prob.QP.c;

if isempty(c), c=zeros(n,1); end

if isempty(Prob.A)
   % Construct one linear constraint from some bounded variable when m==0
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);

   bl=[bl;max(bl(i),-BIG)];
   bu=[bu;min(bu(i),BIG)];
   Prob.A = sparse(zeros(1,n));
   Prob.A(1,i)=1;
   if length(hs) < n+1
      hs=[hs;0];
   end
end
if PriLev > 2
   fprintf('# of non zero entries in linear constraints %d. ',nnz(Prob.A))
   if PriLev > 5
      fprintf('Non linear variables %d. ',n)
      fprintf('Total number of constraints %d.',m)
   end
   fprintf('\n');
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
      %xprint(x_0,'x_0: ',' %14.9f',5);
      %xprinte(bl,'bl: ');
      %xprinte(bu,'bu: ');
      mPrint([bl(1:n) x_0 bu(1:n)],'xl/x0/xu:')
      mPrint([bl(n+1:n+m) bu(n+1:n+m)],'bl/bu:')
      if PriLev >= 10, pause; end
   end
end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

if isempty(c)
   if nnObj==0
      Result.f_0=0;
   else
      Result.f_0=0.5*x_0(1:nnObj)'*Prob.QP.F*x_0(1:nnObj);
   end
elseif nnObj==0
   Result.f_0=c(1:n)'*x_0(1:n);
else
   Result.f_0=0.5*x_0(1:nnObj)'*Prob.QP.F*x_0(1:nnObj) + c(1:n)'*x_0(1:n);
end

%if 0
%   if 1
%      % Put c last in A - WORKS!!!
%      A=[Prob.A;c'];
%      iObj=size(A,1);
%      bl=[bl;-BIG];
%      bu=[bu; BIG];
%   else
%      % Put c first in A - WORKS!!!
%      A=[c';Prob.A];
%      iObj=1;
%      bl=[-BIG;bl];
%      bu=[ BIG;bu];
%   end
%   if 0
%      % Must set c to 0, otherwise 2*c as linear objective
%      c=zeros(n,1);
%   else
%      if 1
%         % Still keep objective in c, not in A - WORKS
%         % Just zero in last row of A
%         A(iObj,:)=0;
%      else
%         %Half in both - WORKS
%         A(iObj,:)=A(iObj,:)/2;
%         c=c/2;
%      end
%   end
%else
%   % Objective in c. iObj == 0 - WORKS
%   A=sparse(Prob.A);
%end

% Objective in c. iObj == 0 
A=sparse(Prob.A);

[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('sqopt',2,...
         nnObj, 0, size(A,1), Prob);

if callback
   [xs, hs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount] = ...
        sqopt7( A, bl, bu, 'HxFunc', c, hElast, iObj, optPar, ...
               Warm, hs, xs,  nS, ...
               SpecsFile, PrintFile, SummFile, ...
               ObjAdd, moremem, Prob.Name, Prob, PriLev );

else
   [xs, hs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount] = ...
        sqopt7( A, bl, bu, Prob.QP.F, c, hElast, iObj, optPar, ...
               Warm, hs, xs,  nS, ...
               SpecsFile, PrintFile, SummFile, ...
               ObjAdd, moremem, Prob.Name, Prob, PriLev );
end

% This is the call if using a callback routine HxFunc to compute
% z= H*x given x, where H is the QP matrix
%
% [xs, hs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount] = ...
%      sqopt( A, bl, bu, 'HxFunc', c, hElast, iObj, optPar, ...
%             Warm, hs, xs,  nS, ...
%             SpecsFile, PrintFile, SummFile, ...
%             ObjAdd, moremem, Prob.Name, Prob );

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

if ~isempty(c)
   if nnObj > 0
      Result.g_k=Prob.QP.F*xs(1:n)+c;
   else
      Result.g_k=c;
   end
elseif nnObj > 0
   Result.g_k=Prob.QP.F*xs(1:n);
else
   Result.g_k=zeros(n,1);
end

% Saved for a warm start
% The dummy objective row is last in xs, value -Obj, last in hs is 3 (basic)
Result.SOL.xs=[xs;-Obj];
Result.SOL.hs=[hs;3];
Result.SOL.nS=nS;
Result.SOL.nInf=nInf;
Result.SOL.sInf=sInf;
Result.SOL.optPar=optPar;

Iter=iwCount(1);

Result.FuncEv    = 0;
Result.GradEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
%Result.MajorIter = iwCount(2);
%Result.MinorIter = iwCount(2);
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

if m > 0
   Result.Ax = A * xs(1:n);
   Result = StateDef(Result, xs(1:n), Result.Ax, [], ...
                     Prob.optParam.xTol, Prob.optParam.bTol, [], bl, bu, 1);
else
   % Compute Result.xState and Result.QP.B only
   Result = StateDef(Result, xs(1:n),[],[], Prob.optParam.xTol,[],[], bl,bu,1);
end



switch Inform
  case 0, Text = 'Finished successfully';
  case 1, Text = 'Optimality conditions satisfied';
  case 2, Text = 'Feasible point found';
  case 4, Text = 'Weak QP minimizer';
   %     
  case 10, Text = 'The problem appears to be infeasible';
  case 11, Text = 'Infeasible linear constraints';
  case 12, Text = 'Infeasible linear equalities';
  case 14, Text = 'Infeasibilities minimized';
   %     
  case 20, Text = 'The problem appears to be unbounded';
  case 21, Text = 'Unbounded objective';
   %     
  case 30, Text = 'Resource limit error';
  case 31, Text = 'Iteration limit reached';
  case 33, Text = 'The superbasics limit is too small';
   %     
  case 40, Text = 'Terminated after numerical difficulties';
  case 42, Text = 'Singular basis';
  case 43, Text = 'Cannot satisfy the general constraints';
  case 44, Text = 'Ill-conditioned null-space basis';
   %     
  case 50, Text = 'Error in the user-supplied functions';
  case 53, Text = 'The QP Hessian is indefinite';
   %     
  case 70, Text = 'User requested termination';
  case 73, Text = 'Terminated during QP objective evaluation';
  case 74, Text = 'Terminated from monitor routine';
   %     
  case 80, Text = 'Insufficient storage allocated';
  case 81, Text = 'Work arrays must have at least 500 elements';
  case 82, Text = 'Not enough character storage';
  case 83, Text = 'Not enough integer storage';
  case 84, Text = 'Not enough real storage';
   %     
  case 90, Text = 'Input arguments out of range';
  case 91, Text = 'Invalid input argument';
  case 92, Text = 'Basis file dimensions do not match this problem';
  case 93, Text = 'The QP Hessian is indefinite';
   %     
  case 140, Text = 'System error';
  case 141, Text = 'Wrong no of basic variables';
  case 142, Text = 'Error in basis package';
  
 otherwise, Text = 'UNKNOWN SQOPT Inform value.';
  
end
  



% switch Inform
%    case 0
%       Text = 'Optimal solution found';
%    case 1
%       Text = 'Problem Infeasible';
%    case 2
%       Text = 'Problem Unbounded (or badly scaled)';
%    case 3
%       Text = 'Too many iterations';
%    case 4
%       Text = str2mat(' The QP Hessian H appears to be indefinite ' ...
%              ,'(the QP is non-convex)');
%    case 5
%       Text = 'The Superbasics limit is too small';
%    case 6
%       Text = str2mat('A weak solution has been found,' ...
%                     ,'i.e. the solution is not unique');
%    case 10
%       Text = str2mat('Numerical error in trying to satisfy' ...
%                     ,'the constraints Ax-s=0.' ...
%                     ,'The basis is very ill-conditioned.');
%    case 20
%       Text = 'Not enough storage for the basis factorization';
%    case 21
%       Text = 'Error in basis package';
%    case 22
%       Text = str2mat('The basis is singular after several attempts' ...
%             ,'to factorize it (and add slacks where necessary)');
%    case 30
%       Text = str2mat('An OLD BASIS file had dimensions that did' ...
%              ,'not match the current problem');
%    case 32
%       Text = 'System error. Wrong number of basic variables';
%    case 42
%       Text = 'Not enough 8-character workspace to solve the problem';
%    case 43
%       Text = 'Not enough integer workspace to solve the problem';
%    case 44
%       Text = 'Not enough real workspace to solve the problem';
%    otherwise
%       Text = 'NOTE: UNKNOWN SQOPT Inform value.';
% end
% 
Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nSQOPT solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('SQOPT: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');
   
   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Sum of iterations%7d. ',iwCount(1));
   %fprintf('Major  iterations%7d. ',iwCount(2));
   %fprintf('Minor iterations%7d. ',iwCount(2));
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

      fprintf('Reduced costs rc: Last %d elements should be v_k\n',m);
      xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'rc: ',' %14.9f',5);
   end
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000718 hkh  New sqoptMex interface.
% 000916 hkh  Return string ExitText with interpretation of Inform flag 
% 010403 hkh  Change superbasics to #48
% 010412 hkh  Cleaning
% 010419 hkh  Fetching moremem from Prob.SOL
% 010710 hkh  Error in comments
% 010903 hkh  Add complete pivoting option. Allow cold start change of hs,ns
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 020227 hkh  Add Prob structure last in output, optional for SQOPT
% 020417 hkh  Send optPar back in Result.SOL
% 020630 hkh  Use Prob.SOL.callback to enable Matlab callback for QP.F*x
% 021217 hkh  Only use Iter, sum of iterations
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040127 ango Descriptive text changed from 6.1-1 to 6.2-2
% 040303 hkh  Safe definition of nS and hs
% 040602 med  Help fix PriLev to PriLevOpt
% 041203 hkh  Revise calls to defblbu and StateDef
% 041203 hkh  Define Result.Ax
% 041222 med  Safeguard added for x_0