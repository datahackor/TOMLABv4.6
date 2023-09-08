% mipSolve.m:
%
% function Result = mipSolve(Prob)
%
% Branch & Bound algorithm for Mixed-Integer Programming (MIP) 
% using LP relaxation (Formulated as min-IP)
%
% Solving MIP on the LP form with generalized lower and upper bounds: 
% 
%        min    c' * x.  x in R^n, b_L <= A x <= b_U, x_L <= x <= x_U
%
% Any of the x could be set as integer valued
%
% Note that the TOMLAB routine cpTransf.m can be used to rewrite
% an LP problem on another form to this standard form
%
% INPUT PARAMETERS
% Fields in Prob:
%   c:      The vector c in c'x. 
%   A:      The linear constraint matrix 
%   b_L:    Lower bounds on linear constraints. 
%           If empty, assumed to be == b_U
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x. If empty assumed to be 0.
%   x_U:    Upper bounds on x. If empty assumed to be Inf.
%
%   x_0:    Starting point x (If EMPTY, the LP solver solves a Phase I LP 
%           to find a feasible point. Some LP solvers chooses own x_0).
% MaxCPU:   Maximal CPU Time (in seconds) to be used by mipSolve, 
%           stops with best point found
%
% PriLevOpt Printing level (in lpSolve and DualSolve, or other interfaces):
%           =0 No output; >0 Convergence results; 
%           >1 Output every iteration  >2 Output each step in the simplex alg
%           For other LP solvers, see the documentation for the solver
%
% QP.B:     Active set B_0 at start. 
%           1  = Include variable x(i) is in basic set.
%           0  = Variable x(i) is set on its lower bound
%           -1 = Variable x(i) is set on its upper bound. NOT USED
%           If EMPTY, lpSolve finds active set (if used for LP solution).
%
% SolverLP  Name of the solver used for initial LP subproblem. If empty,
%           the default solver is used. See GetSolver.m and tomSolve.m
% SolverDLP Name of the solver used for dual LP subproblems. If empty,
%           the default solver is used. See GetSolver.m and tomSolve.m
% SOL.optPar Parameters for the SOL solvers. If this field is set it is
%           sent to the SOL solvers (if they are used as solvers)
%           See help for minosTL.m, lpoptTL.m for how to set these parameters
% SOL.PrintFile Print file name for the SOL solvers. If this field is set it is
%           sent to the SOL solvers (if they are used as solvers)
%
% ----------------------------------------------------------------------
% Special fields for MIP in Prob.MIP
%
% IntVars:  If IntVars is a scalar, then variables 1,...,IntVars are 
%           assumed to be integers. 
%           If empty, all variables are assumed non-integer (LP problem)
%           If length(IntVars) > 1 ==> length(IntVars) == length(c) should hold,
%           Then IntVars(i) ==1 ==> x(i) integer. IntVars(i) ==0 ==> x(i) real.
%           If length(IntVars) < n, IntVars is assumed to be a set of indices.
%
% VarWeight:Weight for each variable in the variable selection phase.
%           A lower value gives higher priority. Setting
%           Prob.MIP.VarWeight = c; for knapsack problems improve convergence.
% DualGap   mipSolve stops if the duality gap is less than DualGap
%           DualGap = 1, stop at first integer solution
%           e.g. DualGap = 0.01, stop if solution < 1% from optimal solution
% fIP       An upper bound on the IP value wanted. Makes it possible to
%           cut branches and avoid node computations.
% xIP       The x-values giving the fIP value.
% KNAPSACK  True if a knapsack problem is to be solved, 
%           then a knapsack heuristic is used.
% ----------------------------------------------------------------------
%
% Prob.Solver.Alg Node selection method
%           = 0 Depth First (default)
%           = 1 Breadth First
%           = 2 Depth First. When integer solution found, use Breadth selection
%
% Prob.Solver.Method  Rule to select new variables in DualSolve/lpSolve:
%           Note - this option is not used for other LP solvers.
%           = 0 Minimum Reduced Cost, sort variables increasing. (DEFAULT)
%           = 1 Bland's Anti-cycling Rule 
%           = 2 Minimum Reduced Cost, Dantzig's rule 
% PriLev    Print level in mipSolve
%
% Fields used in Prob.optParam, in lpSolve and DualSolve (if they are used)
% IterPrint Print one line each iteration
% MaxIter   Maximal number of iterations. max(10*dim(x),100) is DEFAULT.
% wait      Wait flag, pause each iteration if set true
% eps_f     Tolerance used for function value tests
% eps_Rank  Rank tolerance
%
%
% OUTPUT PARAMETERS
% Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag 
%            == 0  => Global optimal solution found, or integer solution with
%                     duality gap less than user tolerance
%            == 1  => Maximal number of iterations reached.
%            == 2  => Empty feasible set, no integer solution found.
%            == 3  => Rank problems (not used)
%            == 4  => No feasible point found running LP relaxation
%            == 5  => Illegal x_0 found in LP relaxation
%            == 99 => Maximal CPU Time used (cputime > Prob.MaxCPU)
%   ExitTest Text string giving ExitFlag and Inform information
%   Inform   If ExitFlag > 0, Inform=ExitFlag, except if
%            duality gap less than given tolerance (Inform = 6)
%   DualGap  Relative duality gap, (fDual - f_k) / f_k, if f_k ~=0
%            (fDual - f_k) if f_k == 0. If f_k ~=0:
%            Scale with 100, 100*DualGap, to get the percentage duality gap.
%            Absolute value duality gap: scale with f_k, f_k * DualGap
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   QP.B     B  Optimal set. B(i)==0, include variable x(i) in basic set.
%            sum(B==1)==length(b)  holds. See QP.B as input.
%   QP.y     Dual parameters y (also part of v_k)
%   p_dx     Search steps in x
%   alphaV   Step lengths for each search step
%   f_k      Function value c'*x
%   g_k      Gradient c
%   Solver   mipSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 1994-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Nov 30, 1994.   Last modified Jan 17, 2005.

% This ExitFlag currently not used, running standard LP instead
%            == 6  => Starting point x_0 and B not dual feasible. 
%                     Some reduced cost is negative.

function ResultMIP = mipSolve(Prob)

if nargin < 1
   error('mipSolve needs input structure Prob');
end

DEBUG=0;

solvType=checkType('mip');

Prob=ProbCheck(Prob,'mipSolve',solvType);

Prob = iniSolve(Prob,solvType,0,0);

Prob.QP.F = [];

A   = Prob.A;
b   = Prob.b_U(:);
c   = Prob.QP.c(:);

MaxCPU = DefPar(Prob,'MaxCPU',Inf);

if isempty(c)
   disp('Empty objective function coefficients in Prob.QP.c');
   error('mipSolve: Not possible to solve MIP without objective function');
end

[m, n] = size(A);

nb=m+n+1;

x = Prob.x_0(:);  
if any(isinf(x) | isnan(x)), x=[]; end

B = Prob.QP.B(:);

x_L = Prob.x_L(:);
if isempty(x_L), x_L=zeros(n,1); end
x_U = Prob.x_U(:);
if isempty(x_U), x_U=Inf*ones(n,1); end

if length(x_L) < n
   if length(x_L)==1 
      x_L=x_L*ones(n,1);
   else 
      x_L=[x_L;zeros(n-length(x_L),1)]; 
   end
end
if length(x_U) < n
   if length(x_U)==1 
      x_U=x_U*ones(n,1);
   else
      x_U=[x_U;Inf*ones(n-length(x_U),1)]; 
   end
end

% Safe-guard starting point
x   = max(x_L(1:n),min(x,x_U(1:n)));
x_L0=x_L;
x_U0=x_U;

%b_L = b;
%Prob.b_L = b;

if m > 0
   if length(Prob.b_U) < m
      Prob.b_U=[Prob.b_U;inf*ones(m-length(Prob.b_U),1)];
   else
      Prob.b_U=Prob.b_U(:);
   end
   if isempty(Prob.b_L)
      Prob.b_L=Prob.b_U(:);
   elseif length(Prob.b_L) < m
      Prob.b_L=[Prob.b_L;-inf*ones(m-length(Prob.b_L),1)];
   else
      Prob.b_L=Prob.b_L(:);
   end
end

b_L = Prob.b_L;
b_U = Prob.b_U;

%if 0
%% Check if any slacks are needed
%
%ixU = find(Prob.b_L ~= Prob.b_U & ~isinf(Prob.b_U));
%ixL = find(Prob.b_L ~= Prob.b_U & ~isinf(Prob.b_L));
%
%if isempty(ixL) & isempty(ixU)
%   N   = n;
%   A   = Prob.A;
%   b   = Prob.b_U(:);
%   c   = Prob.QP.c(:);
%else
%   mL = length(ixL);
%   mU = length(ixU);
%   if isempty(ixL)
%      A   = [Prob.A,zeros(m,mU)];
%      A(ixL,n+1:n+mU)   = eye(mU);
%      b   = Prob.b_U(ixL);
%      c   = [Prob.QP.c(:);zeros(mU,1)];
%   elseif isempty(ixU)
%      A   = [Prob.A,zeros(m,mL)];
%      A(ixL,n+1:n+mU)   = -eye(mL);
%      b   = Prob.b_U(:);
%      c   = [Prob.QP.c(:);zeros(mL,1)];
%   end
%   A   = Prob.A;
%   b   = Prob.b_U(:);
%   c   = Prob.QP.c(:);
%end
%end

Prob.N   = n;
Prob.x_L = x_L;
Prob.x_U = x_U;

optParam=Prob.optParam;

% Pick up SOL parameters
if isfield(Prob.SOL,'optPar')
   optPar=Prob.SOL.optPar;
else
   optPar=[];
end
if isfield(Prob.SOL,'PrintFile')
   PrintFile=Prob.SOL.PrintFile;
else
   PrintFile=[];
end

PriLev    = Prob.PriLev;          % Print level
if isempty(PriLev), PriLev=1; end
PriLevOpt = Prob.PriLevOpt;       % Print level in LP and Dual LP solver
IterPrint = optParam.IterPrint;   % Print short information each iteration

wait     = optParam.wait;         % Pause after printout if 1
%epsRank = optParam.eps_Rank;     % Rank test tolerance
MaxIter  = optParam.MaxIter;      % Maximal number of iterations
eps_f    = optParam.eps_f;    

% NOT USED xTol      Tolerance to judge if x-values are close
%xTol     = optParam.xTol;
%if isempty(xTol), xTol = 100*eps; end

PriLev    = Prob.PriLev;          % Print level
if isempty(PriLev), PriLev=1; end
PriLevOpt = Prob.PriLevOpt;       % Print level in LP and Dual LP solver
IterPrint = optParam.IterPrint;   % Print short information each iteration

if isfield(Prob.MIP,'IntVars')
   IntVars=Prob.MIP.IntVars;
else
   IntVars=[];
end
if isempty(IntVars) 
   Ivar=[];
elseif length(IntVars)==1
   Ivar=[1:min(n,max(0,floor(IntVars)))];
elseif length(IntVars) < n
   Ivar=IntVars;
else
   Ivar=find(IntVars > 0);
end

eps_1 = 1E-8;  % Level when a variable is considered as an integer value

if ~isfield(Prob.MIP,'fIP')
   Prob.MIP.fIP = [];
end
if ~isfield(Prob.MIP,'xIP')
   Prob.MIP.xIP = [];
end
if isempty(Prob.MIP.fIP)
   fIPMax = Inf;     % Current BEST integer solution
   xIPMax = Prob.MIP.xIP(:);
   if ~isempty(xIPMax)
      % User has supplied an x
      if length(xIPMax) < n
         xIPMax=[xIPMax;zeros(n-length(xIPMax))];
      end
      xIPMax(Ivar) = round(xIPMax(Ivar));
      fIPMax= c'*xIPMax;
   end
else
   fIPMax = Prob.MIP.fIP;
   xIPMax = Prob.MIP.xIP;
   if isempty(xIPMax)
      fIPMax = Inf;     % Current BEST integer solution
   else
      xIPMax(Ivar) = round(xIPMax(Ivar));
   end
end
% Check if xIPMax is feasible
if ~isinf(fIPMax) 
   Ax = A*xIPMax;
   if any(xIPMax < x_L | xIPMax > x_U) | ...
      any(b_L-Ax > 1E-8) | any(Ax-b_U > 1E-8)
      if IterPrint | PriLev > 0
         disp('xIPMax is not feasible');
      end
      fIPMax = Inf;     
      xIPMax = [];     
   end
end

if isfield(Prob.MIP, 'DualGap')
   if isempty(Prob.MIP.DualGap)
      DualGap = 0;
   else
      DualGap = Prob.MIP.DualGap;
   end
else
   DualGap = 0;
end
BIPMax=[];
yIPMax=[];
vIPMax=[];

if ~isfield(Prob.MIP,'VarWht')
   Prob.MIP.VarWht = [];
   VarWht=[];
else
   VarWht=Prob.MIP.VarWeight(:);
end
if isfield(Prob.MIP,'KNAPSACK')
   if isempty(Prob.MIP.KNAPSACK)
      KNAPSACK=0;
   else
      KNAPSACK=Prob.MIP.KNAPSACK;
   end
else
   KNAPSACK=2;
end

ResultMIP                 = ResultDef(Prob);
ResultMIP.Solver          = 'mipSolve';
ResultMIP.Prob            = Prob;
ResultMIP.SolverAlgorithm = 'Branch & bound.';

% Node selection method
Alg=Prob.Solver.Alg;

if isempty(Alg), Alg=0; end

Prob.Solver.Alg=Alg;

if Alg==1
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' Breadth First.'];
elseif Alg==2
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' Depth First, then Breadth.'];
else
   Alg=0;
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' Depth First.'];
end

% Index selection rule in DualSolve and lpSolve

Method=Prob.Solver.Method;

if isempty(Method), Method=0; end

Prob.Solver.Method=Method;

if KNAPSACK > 0
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
         ' Knapsack heuristic.'];
end

if ~isempty(VarWht)
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
         ' Priority weights.'];
end

%if Method==1
%   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
%         ' (Bland''s rule).'];
%elseif Method==2
%   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
%         ' (Minimum reduced cost. Dantzig''s rule).'];
%else
%   Method=0;
%   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
%         ' (Minimum reduced cost).'];
%end

if PriLev > 0
   fprintf('%s',ResultMIP.SolverAlgorithm);
   fprintf('\n');
end

IndexRule = Prob.Solver.Method;

ResultMIP.f_0 = [];
ResultMIP.x_k  = [];
ResultMIP.f_k  = [];
ResultMIP.v_k  = [];
ResultMIP.QP.B = [];
ResultMIP.QP.y = [];

if isempty(Prob.SolverLP)
   %SolverLP=GetSolver('lp',Prob.LargeScale);
   % Avoid use of LPOPT by setting LargeScale true
   SolverLP=GetSolver('lp',1);
   Prob.SolverLP = SolverLP;
else
   SolverLP=Prob.SolverLP;
end
if isempty(Prob.SolverDLP)
   %SolverDLP=GetSolver('dlp',Prob.LargeScale);
   % Avoid use of LPOPT by setting LargeScale true
   SolverDLP=GetSolver('dlp',1);
   Prob.SolverDLP = SolverDLP;
else
   SolverDLP=Prob.SolverDLP;
end

% Setup Structure used in LP calls
ProbLP  = CreateProbQP(Prob, 'lp', max(5000,10*n), PriLev-3, optParam.wait);

%HKH
%ProbLP.b_U  = b;
%ProbLP.b_L  = b; 
ProbLP.b_L  = b_L; 
ProbLP.b_U  = b_U;
ProbLP.A    = A;
ProbLP.QP.c = c;
ProbLP.x_L  = x_L; 
ProbLP.x_U  = x_U; 
ProbLP.P    = 0; 

%ProbLP.QP.B = B;
%ProbLP.x_0   =[];
%ProbLP.y    = [];
%ProbLP.QP.DualLimit = fIPMax; % Try stop dual iterations early

if ~isempty(optPar)
   ProbLP.SOL.optPar=optPar;
end
if ~isempty(PrintFile)
   ProbLP.SOL.PrintFile=PrintFile;
end

% Setup Structure used in dual LP calls
ProbDLP = CreateProbQP(Prob, 'dlp', max(5000,10*n), PriLev-3, optParam.wait);
ProbDLP.PriLevOpt   = PriLevOpt;
if ~isempty(optPar)
   ProbDLP.SOL.optPar=optPar;
end
if ~isempty(PrintFile)
   ProbDLP.SOL.PrintFile=PrintFile;
end

% Initial step
if PriLev > 2 & PriLevOpt >= 1
   disp('=+=+=+=+=+=+=+=')
   fprintf('=== mipSolve:    LP relaxation. ')
   fprintf('Call %s solving Phase I and II:\n',SolverLP)
   disp('=+=+=+=+=+=+=+=')
end

if PriLev ==2 & PriLevOpt >=2
   ProbLP.PriLevOpt   = PriLevOpt-1;
else
   ProbLP.PriLevOpt   = max(0,PriLevOpt-1);
end

ResultLP = tomSolve(SolverLP,ProbLP);

% Put back previous value of PriLevOpt
%Prob.PriLevOpt = PriLevOpt;

ExitFlag = ResultLP.ExitFlag;
Inform   = ResultLP.Inform;
x        = min(x_U,max(x_L,ResultLP.x_k));
v_k      = ResultLP.v_k;
if ExitFlag==0
   y = v_k(n+1:n+m);
else
   y   = [];
   SOL = [];
end
if isempty(ResultLP.SOL) 
   SOL=0;
elseif strcmpi(SolverDLP,'minos')  | strcmpi(SolverDLP,'sqopt') | ... 
       strcmpi(SolverDLP,'lpopt')  | strcmpi(SolverDLP,'lssol') | ... 
       strcmpi(SolverDLP,'nlssol') | strcmpi(SolverDLP,'npsol') | ... 
       strcmpi(SolverDLP,'snopt') | strcmpi(SolverDLP,'qpopt') 
   SOL=1;
   ProbDLP.WarmStart=1;
else
   SOL=0;
end

% Variables on upper bound must be set as basic variables for DualSolve
B        = ResultLP.QP.B;
fLP0     = ResultLP.f_k;
Iter     = ResultLP.Iter;
SumIter  = Iter;


if ExitFlag > 0 & ExitFlag ~= 3
   if PriLev >= 2
      fprintf('No solution found to LP relaxation. ')
      fprintf('ExitFlag %d\n', ExitFlag)
   end
   x(Ivar) = round(x(Ivar));
   ResultMIP.x_k  = x;
   ResultMIP.f_k  = fLP0;
   ResultMIP.v_k  = v_k;
   ResultMIP.QP.B = B;
   ResultMIP.QP.y = y;
   ResultMIP.ExitFlag=4;
   ResultMIP.ExitText=ExitText(4);
   ResultMIP.Iter=0;
   ResultMIP.MinorIter=SumIter;
   ResultMIP=endSolve(Prob,ResultMIP);
   return;
end
%HKH NEW
if ~isinf(fIPMax)
    if fIPMax ~= 0
       Gap =  abs((fIPMax-fLP0)/fIPMax);
    else
       Gap =  abs(fIPMax-fLP0);
    end
    if Gap <= DualGap
       ENDTREE=1;
       GapOK = 1;
    else
       GapOK = 0;
    end
else
    GapOK = 0;
    Gap   = Inf;
end


iB = find(~B);
if PriLev > 2
   disp('Phase II solution x:')
   xprint(x,'x:')
end

XX = spalloc(nb,MaxIter+2,min(40000,nb*MaxIter));
BB = spalloc(n, MaxIter+2,min(40000,n *MaxIter));
HS = spalloc(nb,MaxIter+2,min(40000,nb*MaxIter));

BB(1:n,1) = sparse(B(:));
if SOL==0
   XX(1:n,1) = sparse(x);
else
   XX(1:nb,1) = sparse(ResultLP.SOL.xs);
   HS(1:nb,1) = sparse(ResultLP.SOL.hs);
end

% Initialization

L = 1;             % Node List
NODE=1;

% Lowest value possible at node, from LP relaxation
f_min = Inf*ones(MaxIter,1);      
f_min(1)=fLP0;

fDual = -Inf;      % Current best dual value for feasible points
pred  = zeros(MaxIter,1);
Problem.xVar  = zeros(MaxIter,1);
Problem.Bound = zeros(MaxIter,1);
Problem.xBound= zeros(MaxIter,1);

fDP=-Inf;
Iter   = 0;
cpumax = 0;
TIME0  = Prob.TIME0;

%PriLev=2
%wait=1

format compact
while ~(Iter > MaxIter)
   Iter = Iter+1;
   if isempty(L)
      if fIPMax == Inf % Empty feasible set
         if PriLev >= 2
            fprintf('\nTotal number of LP iterations = %d\n',SumIter)
            disp('----------------------------------')
            disp('No feasible solution to IP problem')
            disp('----------------------------------')
         end
         ResultMIP.x_k  = [];
         ResultMIP.f_k  = fLP0;
         ResultMIP.v_k  = v_k;
         ResultMIP.QP.B = [];
         ResultMIP.QP.y = y;
         ResultMIP.ExitFlag=2;
         ResultMIP.ExitText=ExitText(2);
         ResultMIP.Iter=Iter;
         ResultMIP.MinorIter=SumIter;
         ResultMIP.DualGap=Gap;
         ResultMIP=endSolve(Prob,ResultMIP);
         return;
      else
         ResultMIP.x_k  = xIPMax;
         ResultMIP.f_k  = fIPMax;
         ResultMIP.v_k  = vIPMax;
         ResultMIP.QP.B = BIPMax;
         ResultMIP.QP.y = yIPMax;
         if GapOK == 0
            ResultMIP.ExitFlag=0;
            ResultMIP.ExitText=ExitText(0);
         else
            ResultMIP.ExitFlag=0;
            ResultMIP.ExitText=ExitText(6);
         end
         ResultMIP.Iter=Iter;
         ResultMIP.MinorIter=SumIter;
         ResultMIP.DualGap=0;
         if PriLev >= 1
            fprintf('\n--- Branch & Bound converged! ')
            fprintf('Iterations (nodes visited) = %5.0f ',Iter)
            fprintf('Total LP Iterations =%7.0f\n',SumIter)
            fprintf('\n    Optimal Objective function =%28.16f\n',fIPMax);
            if KNAPSACK
               xprinti(xIPMax,'x:',5,15)
            else
               xprinti(xIPMax,'x:',7)
            end
         end
         if PriLev == 1
            xprinti(BIPMax,'B:');
         end
         ResultMIP=endSolve(Prob,ResultMIP);
         return;
      end
   elseif GapOK == 1
      ResultMIP.x_k  = xIPMax;
      ResultMIP.f_k  = fIPMax;
      ResultMIP.v_k  = vIPMax;
      ResultMIP.QP.B = BIPMax;
      ResultMIP.QP.y = yIPMax;
      ResultMIP.ExitFlag=0;
      ResultMIP.ExitText=ExitText(6);
      ResultMIP.Iter=Iter;
      ResultMIP.MinorIter=SumIter;
      ResultMIP.DualGap=Gap;
      if IterPrint
         fprintf('Iter %6.0f ',Iter);
         fprintf('              ');
         fprintf('Best f_IP %8d ',round(fIPMax));
         fprintf('Primal bound on IP %8.1f ',fLP0);
         fprintf('Dual gap %8.2f ',fLP0-fIPMax);
         if fIPMax ~= 0
            fprintf(' %6.3f%% ',100*Gap);
         end
         fprintf('\n');
      end
      if PriLev > 1
         fprintf('\n--- User defined duality gap reached ');
         fprintf('Iterations (nodes visited) = %5.0f ',Iter)
         fprintf('Total LP Iterations =%7.0f\n',SumIter)
         fprintf('\n    Best Integer Value =%28.16f\n',fIPMax);
         if KNAPSACK
            xprinti(xIPMax,'x:',5,15)
         else
            xprinti(xIPMax,'x:',7)
         end
      end
      if PriLev == 1
         xprinti(BIPMax,'B:');
      end
      ResultMIP=endSolve(Prob,ResultMIP);
      return;
   end
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end
   if PriLev > 1
      fprintf('\n========== mipSolve ==========   Iteration %6.0f\n',Iter);
   end
   % Problem selection and relaxation
   if Alg==0 | (Alg==2 & isinf(fIPMax))
      i = L(length(L));
      L = L(1:length(L)-1);
   else
      i = L(1);
      L = L(2:length(L));
   end
   A = A(1:m,1:n);
   c = c(1:n);
   b = b(1:m);
   x_L = x_L0;
   x_U = x_U0;
   j = i;
   ix = [];
   Z  = [];
   while pred(j) > 1     % Generate LP relaxation problem
      j = pred(j);
      ix = [j;ix];
      k = Problem.xVar(j);

      if Problem.Bound(j)==1
         x_U(k)=min(x_U(k),Problem.xBound(j));
      else
         x_L(k)=max(x_L(k),Problem.xBound(j));
      end

   end
   if i == 1     % 1st step
      x = full(XX(1:n,1));
      B = full(BB( : ,1));

      f_i_RP = fLP0;
      f_min(1) = f_i_RP;
   else
      j = pred(i);
      k = Problem.xVar(i);

      if Problem.Bound(i)==1
         x_U(k)=min(x_U(k),Problem.xBound(i));
      else
         x_L(k)=max(x_L(k),Problem.xBound(i));
      end

      x = full(XX(1:n,j));
      B = full(BB( : ,j));


      if PriLev > 2
         xprint(x,'x:')
      end
      if PriLev > 2 & PriLevOpt >= 2
         disp('=+=+=+=+=+=+=+=')
         disp('=== mipSolve: Solve dual feasible LP problem:')
         disp('=+=+=+=+=+=+=+=')
      end
  

      if 1 & (SOL > 0)
         ProbDLP.SOL.xs=full(XX(1:nb,j));
         ProbDLP.SOL.hs=full(HS(1:nb,j));
      else
         ProbDLP.WarmStart=0;
      end

      if Iter > 12000000
         ProbDLP.PriLevOpt      = 2;
         ProbDLP.optParam.wait  = 1;
      end

      % HKH
      %ProbDLP.b_U  = b;
      %ProbDLP.b_L  = b; % DualSolve must know all eq:s are equalities
      ProbDLP.b_L  = b_L; 
      ProbDLP.b_U  = b_U;
      ProbDLP.A    = A;
      ProbDLP.QP.c = c;
      ProbDLP.QP.B = B;
      ProbDLP.x_0   =x;
      ProbDLP.y    = [];
      ProbDLP.x_L  = x_L; 
      ProbDLP.x_U  = x_U; 
      ProbDLP.QP.DualLimit = fIPMax; % Try stop dual iterations early
      ProbDLP.P    = i; 
%      if Iter > 120000
%         % DEBUG part
%         Prob.SolverDLP=2;
%         Result = tomSolve(SolverDLP,ProbDLP);
%ExitFlag=Result.ExitFlag
%         pause
%         Prob.SolverDLP=0;
%      end

%ProbDLP.PriLevOpt   = 1;

      Result = tomSolve(SolverDLP,ProbDLP);

      SumIter=SumIter+Result.Iter;
      ExitFlag=Result.ExitFlag;
      Inform=Result.Inform;
      if (ExitFlag == 6 | ExitFlag == 1) & strcmpi('DualSolve',SolverDLP)
         % MipSolve gave bad starting point (or too many iterations.
         % Try a standard LP solver
         Result = tomSolve('qld',ProbDLP);
         Inform=Result.Inform;
      end
   
      ExitFlag=Result.ExitFlag;
      x   = min(x_U,max(x_L,Result.x_k));
      v_k = Result.v_k;
      %y   = v_k(n+1:length(v_k));
      y    = Result.y_k;
      fDP  = Result.f_k;
      B   = Result.QP.B;

      if ExitFlag == 0 | ExitFlag == 3
         BB(:,i) = sparse(B);
         if SOL == 0 | isempty(Result.SOL)
            XX(1:n,i) = sparse(x);
         else
            XX(:,i)   = sparse(Result.SOL.xs);
            HS(:,i)   = sparse(Result.SOL.hs);
         end

         if PriLev > 2
            fprintf('Dual simplex; Solution x:\n')
            xprint(x,'x:')
         end
         f_i_RP = c'*x;
         f_min(i) = f_i_RP;

%if 0 & abs(f_i_RP-fDP)/max(1,abs(f_i_RP)) > 1E-12
%   disp('BIG DUALITY GAP');
%format long
%(f_i_RP-fDP)/max(1,abs(f_i_RP))
%f_i_RP
%fDP
%b'*y
%Iter
%pause
%end
      else
         if PriLev > 1
            fprintf('   No Dual Solution, EXIT = %4.0f\n',ExitFlag)
         end
         f_i_RP = Inf;
         f_min(i) = Inf;
      end
   end
   if PriLev > 1
      fprintf('   Objective function ');
      fprintf('f_i_RP == f_LP = %30.16f\n',f_i_RP);
      if ~isinf(fIPMax)
         fprintf('   Best integer solution found so far');
         fprintf(' = %20.2f\n',fIPMax);
         fprintf('   Primal bound on integer solution  ');
         fprintf(' = %20.2f\n',fLP0);
         fprintf('   Duality gap                       ');
         fprintf(' = %20.2f\n',fLP0-fIPMax);
      end
      if ~isinf(fDual) & PriLev > 3
         fprintf('   Dual value b''*y           ');
         fprintf('        = %30.16f\n',fDual);
      end
   elseif IterPrint
      fprintf('Iter %6.0f ',Iter);
      fprintf('f_LP %8.1f ',f_i_RP);
      fprintf('Best f_IP %8d ',round(fIPMax));
      fprintf('Primal bound on IP %8.1f ',fLP0);
      fprintf('Dual gap %8.2f ',fLP0-fIPMax);
      if fIPMax ~= 0
         fprintf(' %6.3f%% ',100*abs((fLP0-fIPMax)/fIPMax));
      end
      fprintf('! %s Inform %d ',SolverDLP,Inform);
      if ~isinf(fDual) 
         fprintf('f_Dual %8.2f',fDual);
      end
      fprintf(' ExitFlag %d ',ExitFlag);
      fprintf('\n');
   end
   ENDTREE=1;
   %if f_i_RP <  fIPMax - eps_f 
   if f_i_RP <  fIPMax
      ENDTREE=0;
      % Check if solution is integer
      if isempty(Ivar)
         Iidx = [];
      else
         idx = find(B(Ivar)==0);
         Iidx  = Ivar(idx);
      end
      if isempty(Iidx)
         IntVar=1;
      else
         x_I = floor(x(Iidx)+eps_1*max(1,abs(x(Iidx))));
         x_r = max(0,x(Iidx)-x_I);

%for iii=1:length(Iidx)
%    xz(iii)=floor(x(Iidx(iii))+eps_1);
%fprintf('%30.20f %30.20f\n',x(Iidx(iii)),...
%         x(Iidx(iii))+eps_1*max(1,abs(x(Iidx(iii)))));
%end

         if isempty(VarWht)
            % Variables with frac.part closest to 0.5
            %[best_frac xBest] = min(abs(x_r-0.5));
            r=abs(x_r-0.5);
            [rBest iBest]=sort(r);
            irBest=find(rBest < 0.5-eps_1);
         else
            r=VarWht(Iidx);
            r(x_r < eps_1*max(1,x_r))=Inf;
            [rBest iBest]=sort(r);
            irBest=find(~isinf(rBest));
         end

         if isempty(irBest)
            IntVar=1;
         else
            IntVar=0;
            iCand = iBest(irBest);
            Cand = Iidx(iCand);

            if PriLev > 1
               fprintf('   %d nonintegers. ',length(irBest));
               if isempty(VarWht)
                  fprintf('Best fraction = %8.6f.',x_r(iCand(1)))
               else
                  fprintf('Weight %8.3f. Fraction %8.6f.',...
                           r(iCand(1)),x_r(iCand(1)))
               end
               fprintf(' Base var = %4.0f. Var # = %4.0f',iCand(1),Cand(1))
               fprintf(' |L| %d\n',length(L));
            end
         end
      end
      if IntVar
         if PriLev > 1
            disp('   Integer solution found!!!');
         end
         if fDP > fDual, fDual=fDP; end
         if f_i_RP <  fIPMax 
            fIPMax = f_i_RP;
            x(Ivar) = round(x(Ivar));
            xIPMax = x(1:n);
            BIPMax = B(1:n);
            yIPMax = y;
            vIPMax = v_k;
            DelL = f_min(L) >= fIPMax;
            %xprinti(L(find(DelL)),'d:');
            XX(:,L(find(DelL))) = 0;
            BB(:,L(find(DelL))) = 0;
            HS(:,L(find(DelL))) = 0;

            L = L(find(f_min(L) < fIPMax));

            if PriLev > 1
               fprintf('   Found new BEST IP. ');
               fprintf('Nodes left %6.0f. ',length(L));
               fprintf('Nodes deleted %6.0f.\n',sum(DelL));
            end
            if PriLev > 3
               fprintf('   New List L after fIPMax check:\n');
               xprint(L,'L:',' %4.0f');
            end
            if abs(fIPMax - c(1:n)'*xIPMax) < 1E-10
            %if fIPMax == c(1:n)'*xIPMax
               ENDTREE=1;
            else
               Cand=[];
            end
            if fIPMax ~= 0
               Gap =  abs((fIPMax-fLP0)/fIPMax);
            else
               Gap =  abs(fIPMax-fLP0);
            end
            if Gap <= DualGap
               ENDTREE=1;
               GapOK = 1;
            end
            if PriLev > 2
               xprint(xIPMax,'x:');
               fprintf('   Best IP function value %30.16f\n',fIPMax);
            end
         end % if 
      elseif KNAPSACK > 0
         % Apply rounding knapsack heuristic
         if KNAPSACK ==2
            xKS = floor(x+eps_1*max(1,abs(x)));
            ix = any(A<0);
            xKS(ix) = ceil(x(ix)-eps_1*max(1,abs(x(ix))));
         else
            xKS = floor(x+eps_1*max(1,abs(x)));
         end
         fKS = xKS(1:n-m)'*c(1:n-m);
         if fKS < fIPMax
            xKS(n-m+1:n) = b - A(:,1:n-m)*xKS(1:n-m);
            if all(xKS(n-m+1:n) >=-eps_1) 
               fIPMax = fKS;
               xKS(Ivar) = round(xKS(Ivar));
               xIPMax = xKS;
               BIPMax = B(1:n);
               yIPMax = [];
               vIPMax = [];
               DelL = f_min(L) >= fIPMax;
               %xprinti(L(find(DelL)),'d:');
               XX(:,L(find(DelL))) = 0;
               BB(:,L(find(DelL))) = 0;
               HS(:,L(find(DelL))) = 0;

               L = L(find(f_min(L) < fIPMax));

               if PriLev > 0
                  fprintf('    Found new BEST Knapsack. ');
                  fprintf('Nodes left %6.0f. ',length(L));
                  fprintf('Nodes deleted %6.0f.\n',sum(DelL));
                  fprintf('    Best IP function value %33.16f\n',fIPMax);
               end

               if PriLev > 3
                  fprintf('   New List L after fIPMax check:\n');
                  xprint(L,'L:',' %4.0f');
               end
            end
         end
      end
   end
   if PriLev > 1
      if ENDTREE
         disp('   End of tree');
      end
   end
   if ~ENDTREE & length(Cand) > 0  % Make divisions
       j=i;
       tree=i;
       while pred(j) > 0     % Generate tree
          j = pred(j);
          tree = [j;tree];
       end
       if isempty(tree)
          BRANCH=1;
          k=1;
          xC=Cand(1);
       else
          k=0;
          BRANCH=0;
       end
       while ~BRANCH & k < length(Cand)
          k=k+1;
          xC=Cand(k);
          ixC=iCand(k);
          %ix=find(xC==Problem.xVar(tree));
          ix=[]; % NOT TEST NOW
          BRANCH=1;
       end
       if BRANCH
          NODE = NODE + 1;
          %fprintf('Nodes %d. xC %d. ixC %d\n',NODE,xC,ixC');
          Problem.xVar(NODE:NODE+1)=xC;
          Problem.Bound(NODE:NODE+1)=[-1; 1];
          Problem.xBound(NODE:NODE+1)=[x_I(ixC)+1; x_I(ixC)];
          L = [L NODE NODE+1];
          f_min(NODE:NODE+1) = f_i_RP;
          pred(NODE:NODE+1) = i;
          tree=[tree;NODE;NODE+1];
          NODE=NODE+1;
          if PriLev > 5
             disp('Pred xVar Bound xBound')
             [pred(tree),Problem.xVar(tree),...
              Problem.Bound(tree),Problem.xBound(tree) ]
          end
          if PriLev > 2
              xprinti(Cand,'Cand:');
              disp('  # Tree Pred xVar xBound Bound           Value')
              for ii=1:length(tree)
                  jj=tree(ii);
                  if jj~=1
                     kk=Problem.xVar(jj);
                     fprintf('%3d %4d %4d %4d %6d %5d %15.7f\n',ii,jj,...
                        pred(jj),kk,Problem.xBound(jj),Problem.Bound(jj),x(kk));
                  end
              end
          end
       end

       if PriLev > 2
          disp('New List L after divisions:');
          xprint(L,'L:',' %4.0f');
       end
       if PriLev > 4
          disp([pred(1:NODE) Problem.xVar(1:NODE) Problem.Bound(1:NODE)...
                Problem.xBound(1:NODE)]);
       end
   end
   if wait & PriLev > 1
      if DEBUG
         pause(1)
      else
         disp('Press any key to continue ...')
         pause
      end
   end
end % while
if PriLev >= 1
   fprintf('\n--- TOO MANY Branch and Bound ITERATIONS. ITER = %d\n',Iter)
   fprintf('\n--- Total number of LP iterations = %d\n',SumIter)
end
ResultMIP.x_k  = xIPMax;
ResultMIP.f_k  = fIPMax;
ResultMIP.v_k  = vIPMax;
ResultMIP.QP.B = BIPMax;
ResultMIP.QP.y = yIPMax;
if cpumax
   ResultMIP.ExitFlag=99;
   ResultMIP.ExitText=ExitText(9);
else
   ResultMIP.ExitFlag=1;
   ResultMIP.ExitText=ExitText(1);
end
ResultMIP.Iter=Iter;
ResultMIP.MinorIter=SumIter;
ResultMIP.DualGap=Gap;
ResultMIP=endSolve(Prob,ResultMIP);

% ------------------------------
function Text = ExitText(Inform)
% ------------------------------

switch  Inform
   case 0
     Text = 'Optimal solution found';
   case 1
     Text = 'Maximal number of iterations reached';
   case 2
     Text = 'Empty feasible set, no integer solution found';
   case 4
     Text = 'No solution found to LP relaxation';
   case 5
     Text = 'Illegal x_0 found in LP relaxation';
   case 6
     Text = 'User defined duality gap reached';
   case 9
     Text = 'User defined maximal CPU Time used';
end
% MODIFICATION LOG
%
% 981111  hkh  Change to use lpSolve
% 981114  hkh  Change to Prob/Result format for lpSolve
% 981115  hkh  Setting lower bounds as one, not -Inf.
% 981119  hkh  Change field to Prob.QP.B. Errors in use of B
% 981123  hkh  Change name to mipSolve. 
% 981123  hkh  Change optPar(28) to ExitFlag on 3 places. max_iter to MaxIter
% 990419  hkh  Test on idx = find(B(1:Ivar)) empty (= int vars 0, non basic).
% 990804  hkh  Change to Prob/Result input-output format
% 990810  hkh  Revised for v2.0. Using new DualSolve and lpSolve.
% 990901  hkh  Calling general SolveDLP.
% 990913  hkh  Safeguard against nan and inf in x_0
% 991222  hkh  Chcck if c is empty, stop.
% 000916  hkh  Add text about convergence in ExitText
% 010407  hkh  IntVars and Ivars made correct
% 020204  hkh  Pick up Prob.SOL.optPar/PrintFile and use if SOL solvers
% 020304  hkh  Always use MINOS instead of LPOPT for LP sub problems.
% 020304  hkh  Field error setting Prob.SOL.optPar. Print Inform from solver
% 021223  hkh  Handle lower and upper bounds on linear constraints
% 030107  hkh  Avoid allocating too much memory for XX, HS, BB
% 030107  hkh  Prob.b_L empty now gives Prob.b_U (as comments say)
% 030107  hkh  Correct comments
% 030203  hkh  Correct = b_U to <= b_U
% 030309  hkh  Ensure all x inside bounds: x = min(x_U,max(x_L,Result.x_k));
% 030823  hkh  Fix bugs in handling duality gap, incl. printout last iteration
% 030831  hkh  Duality gap not checked at start. Check feasibility of xIPMax
% 031111  hkh  Add Duality gap as output field Result.DualGap
% 031128  hkh  Gap undefined if fIPMax given. Set DualGap=0 if convergence
% 040111  hkh  Add call to inisolve and endSolve
% 040425  hkh  Add input SmallA, if 1 detect and remove small A elements
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 041123  hkh  Check MIP fields to handle LP problems
% 041222  med  Safeguard added for x_0
% 050117  med  mlint revision