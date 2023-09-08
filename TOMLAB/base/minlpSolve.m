% minlpSolve.m:
%
% function Result = minlpSolve(Prob)
%
% Branch & Bound algorithm for Mixed-Integer Nonlinear Programming (MINLP) 
% using NLP relaxation (Formulated as minlp-IP)
%
% Solving MINLP with generalized lower and upper bounds: 
% 
%        min    f(x).  x in R^n, b_L <= A x <= b_U, x_L <= x <= x_U
%
% Any of the x could be set as integer valued
%
% INPUT PARAMETERS
% Fields in Prob:
%   A:      The linear constraint matrix 
%   b_L:    Lower bounds on linear constraints. 
%           If empty, assumed to be == b_U
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x. If empty assumed to be 0.
%   x_U:    Upper bounds on x. If empty assumed to be Inf.
%
%   x_0:    Starting point x (If EMPTY, the NLP solver solves a Phase I NLP 
%           to find a feasible point. Some NLP solvers chooses own x_0).
% PriLevOpt Printing level:
%           =0 No output; >0 Convergence results; 
%           >1 Output every iteration  >2 Output each step in the NLP alg
%           For other NLP solvers, see the documentation for the solver
%
% QP.B:     Active set B_0 at start. 
%           1  = Include variable x(i) is in basic set.
%           0  = Variable x(i) is set on its lower bound
%           -1 = Variable x(i) is set on its upper bound. NOT USED
%           If EMPTY, the NLP solver finds active set.
%
% SolverNLP Name of the solver used for initial NLP subproblem. If empty,
%           the default solver is used. See GetSolver.m and tomSolve.m
% SOL.optPar Parameters for the SOL solvers. If this field is set it is
%           sent to the SOL solvers (if they are used as solvers)
%           See help for minosTL.m, npsolTL.m and snoptTL.m for how to 
%           set these parameters
% SOL.PrintFile Print file name for the SOL solvers. If this field is set it is
%           sent to the SOL solvers (if they are used as solvers)
%
% ----------------------------------------------------------------------
% Special fields for MIP in Prob.MIP
%
% IntVars:  If IntVars is a scalar, then variables 1,...,IntVars are 
%           assumed to be integers. 
%           If empty, all variables are assumed non-integer (NLP problem)
%           If length(IntVars) > 1 ==> length(IntVars) == length(x_L) should hold,
%           Then IntVars(i) == 1 ==> x(i) integer. IntVars(i) == 0 ==> x(i) real.
%           If length(IntVars) < n, IntVars is assumed to be a set of indices.
%
% VarWeight:Weight for each variable in the variable selection phase.
%           A lower value gives higher priority. Setting
%           Prob.MIP.VarWeight might improve convergence.
% DualGap   minlpSolve stops if the duality gap is less than DualGap
%           DualGap = 1, stop at first integer solution
%           e.g. DualGap = 0.01, stop if solution < 1% from optimal solution
% fIP       An upper bound on the IP value wanted. Makes it possible to
%           cut branches and avoid node computations.
% xIP       The x-values giving the fIP value.
% ----------------------------------------------------------------------
%
% Prob.Solver.Alg Node selection method
%           = 0 Depth First (default)
%           = 1 Breadth First
%           = 2 Depth First. When integer solution found, use Breadth selection
%
% PriLev    Print level in minlpSolve (default 1).
%
% Fields used in Prob.optParam, in sub solvers.
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
%            == 4  => No feasible point found running NLP relaxation
%            == 5  => Illegal x_0 found in NLP relaxation
%   ExitText Text string giving ExitFlag and Inform information
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
%   f_k      Function value at optimum
%   g_k      Gradient
%   Solver   minlpSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 1994-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Mar 31, 2004.   Last modified Jan 17, 2005.

% This ExitFlag currently not used, running standard LP instead
%            == 6  => Starting point x_0 and B not dual feasible. 
%                     Some reduced cost is negative.

function ResultMINLP = minlpSolve(Prob)

if nargin < 1
   error('minlpSolve needs input structure Prob');
end

DEBUG=0;

if ~isfield(Prob,'SolverNLP')
    Prob.SolverNLP = [];
end

if ~isfield(Prob,'SolverDNLP')
    Prob.SolverDNLP = [];
end

solvType=checkType('minlpSolve');

Prob=ProbCheck(Prob,'minlpSolve',solvType);

Prob = iniSolve(Prob,solvType,1,1); % DEPENDS ON SUBSOLVER, SOL ASSUMED
% Avoid call to estimate Lagrange multipliers in PrintResult
Prob.PrintLM = 0;

A     = Prob.A;

if isempty(Prob.USER.f)
   disp('Empty objective.');
   error('minlpSolve: Not possible to solve MINLP without objective function');
end

[m, n] = size(A);
if n == 0
    n = max(length(Prob.x_L), length(Prob.x_U));
end
nc = max(1,max(length(Prob.c_L), length(Prob.c_U)));
nb=m+n+nc; %USED FOR MEMORY ALLOCATION

x = Prob.x_0(:);  
if any(isinf(x) | isnan(x)), x=[]; end

% B = Prob.QP.B(:);

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

x_L0=x_L;
x_U0=x_U;

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

c_L = Prob.c_L;
c_U = Prob.c_U;

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

wait     = optParam.wait;     % Pause after printout if 1
MaxIter  = optParam.MaxIter;  % Maximal number of iterations
eps_f    = optParam.eps_f;    

PriLev    = Prob.PriLev;      % Print level
if isempty(PriLev), PriLev=1; end
PriLevOpt = Prob.PriLevOpt;       % Print level in subsolvers
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

fIPMax = DefPar(Prob.MIP,'fIP',[]);
xIPMax = DefPar(Prob.MIP,'xIP',[]);
if isempty(fIPMax)
   fIPMax = Inf;     % Current BEST integer solution
   if ~isempty(xIPMax)
      % User has supplied an x
      if length(xIPMax) < n
         xIPMax=[xIPMax;zeros(n-length(xIPMax))];
      end
      xIPMax(Ivar) = round(xIPMax(Ivar));
%       fIPMax= c'*xIPMax; REPLACED BY BELOW.
      fIPMax = nlp_f(xIPMax, Prob);
   end
else
   if isempty(xIPMax)
      fIPMax = Inf;     % Current BEST integer solution
   else
      xIPMax(Ivar) = round(xIPMax(Ivar));
   end
end
% Check if xIPMax is feasible
if ~isinf(fIPMax) 
   Ax = A*xIPMax;
   if ~isempty(Prob.USER.c)
       cx = nlp_c(xIPMax,Prob);
   else
       cx = [];
   end
   if any(xIPMax < x_L | xIPMax > x_U) | ...
      any(b_L-Ax > 1E-8) | any(Ax-b_U > 1E-8) | ...
      any(c_L-cx > 1E-6)  | any(c_U-cx > 1E-6)
      % ADDED CHECK FOR NONLINEAR CONSTRAINTS
      if IterPrint | PriLev > 0
         disp('xIPMax is not feasible');
      end
      fIPMax = Inf;     
      xIPMax = [];     
   end
end

if isfield(Prob.MIP,'DualGap')
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

VarWht = DefPar(Prob.MIP,'VarWht',[]);

ResultMINLP                 = ResultDef(Prob);
ResultMINLP.Solver          = 'minlpSolve';
ResultMINLP.Prob            = Prob;
ResultMINLP.SolverAlgorithm = 'Branch & bound.';

% Node selection method
Alg=Prob.Solver.Alg;

if isempty(Alg), Alg=0; end

Prob.Solver.Alg=Alg;

if Alg==1
   ResultMINLP.SolverAlgorithm=[ResultMINLP.SolverAlgorithm ' Breadth First.'];
elseif Alg==2
   ResultMINLP.SolverAlgorithm=[ResultMINLP.SolverAlgorithm ...
         ' Depth First, then Breadth.'];
else
   Alg=0;
   ResultMINLP.SolverAlgorithm=[ResultMINLP.SolverAlgorithm ' Depth First.'];
end

if ~isempty(VarWht)
   ResultMINLP.SolverAlgorithm=[ResultMINLP.SolverAlgorithm ... 
         ' Priority weights.'];
end

if PriLev > 0
   fprintf('%s',ResultMINLP.SolverAlgorithm);
   fprintf('\n');
end

ResultMINLP.f_0 = [];
ResultMINLP.x_k  = [];
ResultMINLP.f_k  = [];
ResultMINLP.v_k  = [];
ResultMINLP.QP.B = [];
ResultMINLP.QP.y = [];

if isempty(Prob.SolverNLP)
   SolverNLP=GetSolver('con',1);
   Prob.SolverNLP = SolverNLP;
else
   SolverNLP=Prob.SolverNLP;
end
if isempty(Prob.SolverDNLP)
   SolverDNLP=GetSolver('con',1);
   Prob.SolverDNLP = SolverDNLP;
else
   SolverDNLP=Prob.SolverDNLP;
end

% Setup Structure used in NLP calls

ProbNLP = Prob; % SAME FOR NOW, IntVars ignored.
ProbNLP.MIP.IntVars = [];

ProbNLP.b_L  = b_L; 
ProbNLP.b_U  = b_U;
ProbNLP.A    = A;
ProbNLP.x_L  = x_L; 
ProbNLP.x_U  = x_U; 

% Setup Structure used in dual NLP calls
ProbDNLP = Prob;
ProbDNLP.MIP.IntVars = [];

ProbDNLP.b_L  = b_L; 
ProbDNLP.b_U  = b_U;
ProbDNLP.A    = A;
ProbDNLP.x_L  = x_L; 
ProbDNLP.x_U  = x_U; 

% Initial step
if PriLev > 2 & PriLevOpt >= 1
   disp('=+=+=+=+=+=+=+=')
   fprintf('=== minlpSolve:    NLP relaxation. ')
   fprintf('Call %s solving Phase I and II:\n',SolverNLP)
   disp('=+=+=+=+=+=+=+=')
end

if PriLev ==2 & PriLevOpt >=2
   ProbNLP.PriLevOpt   = PriLevOpt-1;
else
   ProbNLP.PriLevOpt   = max(0,PriLevOpt-1);
end

ResultNLP = tomSolve(SolverNLP,ProbNLP);

ExitFlag = ResultNLP.ExitFlag;
Inform   = ResultNLP.Inform;
x        = min(x_U,max(x_L,ResultNLP.x_k));
v_k      = ResultNLP.v_k;
if ExitFlag==0
   y = v_k(n+1:n+m);
else
   y   = [];
   SOL = [];
end
if isempty(ResultNLP.SOL) 
   SOL=0;
elseif strcmpi(SolverDNLP,'minos')  | strcmpi(SolverDNLP,'sqopt') | ... 
       strcmpi(SolverDNLP,'lssol') | strcmpi(SolverDNLP,'nlssol') | ...
       strcmpi(SolverDNLP,'npsol') | strcmpi(SolverDNLP,'snopt') | ...
       strcmpi(SolverDNLP,'qpopt') 
   SOL=1;
   ProbDNLP.WarmStart=1;
else
   SOL=0;
end

B        = ResultNLP.QP.B;
fLP0     = ResultNLP.f_k;
Iter     = ResultNLP.Iter;
SumIter  = Iter;

if ExitFlag > 0 & ExitFlag ~= 3
   if PriLev >= 2
      fprintf('No solution found to NLP relaxation. ')
      fprintf('ExitFlag %d\n', ExitFlag)
   end
   x(Ivar) = round(x(Ivar));
   ResultMINLP.x_k  = x;
   ResultMINLP.f_k  = fLP0;
   ResultMINLP.v_k  = v_k;
   ResultMINLP.QP.B = B;
   ResultMINLP.QP.y = y;
   ResultMINLP.ExitFlag=4;
   ResultMINLP.ExitText=ExitText(4);
   ResultMINLP.Iter=0;
   ResultMINLP.MinorIter=SumIter;
   ResultMINLP=endSolve(Prob,ResultMINLP);
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
   XX(1:nb,1) = sparse(ResultNLP.SOL.xs);
   HS(1:nb,1) = sparse(ResultNLP.SOL.hs);
end

% Initialization

L = 1;             % Node List
NODE=1;

% Lowest value possible at node, from NLP relaxation
f_min = Inf*ones(MaxIter,1);      
f_min(1)=fLP0;

fDual = -Inf;      % Current best dual value for feasible points
pred  = zeros(MaxIter,1);
Problem.xVar  = zeros(MaxIter,1);
Problem.Bound = zeros(MaxIter,1);
Problem.xBound= zeros(MaxIter,1);

fDP=-Inf;
Iter = 0;

format compact
while ~(Iter > MaxIter)
   Iter = Iter+1;
   if isempty(L)
      if fIPMax == Inf % Empty feasible set
         if PriLev >= 2
            fprintf('\nTotal number of NLP iterations = %d\n',SumIter)
            disp('----------------------------------')
            disp('No feasible solution to MINLP problem')
            disp('----------------------------------')
         end
         ResultMINLP.x_k  = [];
         ResultMINLP.f_k  = fLP0;
         ResultMINLP.v_k  = v_k;
         ResultMINLP.QP.B = [];
         ResultMINLP.QP.y = y;
         ResultMINLP.ExitFlag=2;
         ResultMINLP.ExitText=ExitText(2);
         ResultMINLP.Iter=Iter;
         ResultMINLP.MinorIter=SumIter;
         ResultMINLP.DualGap=Gap;
         ResultMINLP=endSolve(Prob,ResultMINLP);
         return;
      else
         ResultMINLP.x_k  = xIPMax;
         ResultMINLP.f_k  = fIPMax;
         ResultMINLP.v_k  = vIPMax;
         ResultMINLP.QP.B = BIPMax;
         ResultMINLP.QP.y = yIPMax;
         if GapOK == 0
            ResultMINLP.ExitFlag=0;
            ResultMINLP.ExitText=ExitText(0);
         else
            ResultMINLP.ExitFlag=0;
            ResultMINLP.ExitText=ExitText(6);
         end
         ResultMINLP.Iter=Iter;
         ResultMINLP.MinorIter=SumIter;
         ResultMINLP.DualGap=0;
         if PriLev >= 1
            fprintf('\n--- Branch & Bound converged! ')
            fprintf('Iterations (nodes visited) = %5.0f ',Iter)
            fprintf('Total NLP Iterations =%7.0f\n',SumIter)
            fprintf('\n    Optimal Objective function =%28.16f\n',fIPMax);
         end
         if PriLev == 1
            xprinti(BIPMax,'B:');
         end
         ResultMINLP=endSolve(Prob,ResultMINLP);
         return;
      end
   elseif GapOK == 1
      ResultMINLP.x_k  = xIPMax;
      ResultMINLP.f_k  = fIPMax;
      ResultMINLP.v_k  = vIPMax;
      ResultMINLP.QP.B = BIPMax;
      ResultMINLP.QP.y = yIPMax;
      ResultMINLP.ExitFlag=0;
      ResultMINLP.ExitText=ExitText(6);
      ResultMINLP.Iter=Iter;
      ResultMINLP.MinorIter=SumIter;
      ResultMINLP.DualGap=Gap;
      if IterPrint
         fprintf('Iter %6.0f ',Iter);
         fprintf('              ');
         fprintf('Best f_IP %8d ',round(fIPMax));
         fprintf('Primal bound on MINLP %8.1f ',fLP0);
         fprintf('Dual gap %8.2f ',fLP0-fIPMax);
         if fIPMax ~= 0
            fprintf(' %6.3f%% ',100*Gap);
         end
         fprintf('\n');
      end
      if PriLev > 1
         fprintf('\n--- User defined duality gap reached ');
         fprintf('Iterations (nodes visited) = %5.0f ',Iter)
         fprintf('Total MINLP Iterations =%7.0f\n',SumIter)
         fprintf('\n    Best Integer Value =%28.16f\n',fIPMax);
      end
      if PriLev == 1
         xprinti(BIPMax,'B:');
      end
      ResultMINLP=endSolve(Prob,ResultMINLP);
      return;
   end
   if PriLev > 1
      fprintf('\n========== minlpSolve ==========   Iteration %6.0f\n',Iter);
   end
   % Problem selection and relaxation
   if Alg==0 | (Alg==2 & isinf(fIPMax))
      i = L(length(L));
      L = L(1:length(L)-1);
   else
      i = L(1);
      L = L(2:length(L));
   end
%    A = A(1:m,1:n);
   b_L = b_L(1:m);
   b_U = b_U(1:m);
   x_L = x_L0;
   x_U = x_U0;
   j = i;
   ix = [];
   Z  = [];
   while pred(j) > 1     % Generate NLP relaxation problem
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
         disp('=== minlpSolve: Solve NLP problem:')
         disp('=+=+=+=+=+=+=+=')
      end
  
      if SOL > 0
         ProbDNLP.SOL.xs=full(XX(1:nb,j));
         ProbDNLP.SOL.hs=full(HS(1:nb,j));
      else
         ProbDNLP.WarmStart=0;
      end

      if Iter > 12000000
         ProbDNLP.PriLevOpt      = 2;
         ProbDNLP.optParam.wait  = 1;
      end

      ProbDNLP.b_L  = b_L; 
      ProbDNLP.b_U  = b_U;
      ProbDNLP.A    = A;
      ProbDNLP.QP.B = B;
      ProbDNLP.x_0   =x;
      ProbDNLP.y    = [];
      ProbDNLP.x_L  = x_L; 
      ProbDNLP.x_U  = x_U; 
      ProbDNLP.QP.DualLimit = fIPMax; % Try stop dual iterations early
%       ProbDNLP.P    = i; 

      if ~exist('Result','var') & strcmpi(SolverDNLP, SolverNLP)
        ProbDNLP = WarmDefSOL(SolverDNLP,ProbDNLP,ResultNLP);
      elseif exist('Result','var')
        ProbDNLP = WarmDefSOL(SolverDNLP,ProbDNLP,Result);     
      end
            
      Result = tomSolve(SolverDNLP,ProbDNLP);

      SumIter=SumIter+Result.Iter;
      ExitFlag=Result.ExitFlag;
      Inform=Result.Inform;
      if (ExitFlag == 4)
         % minlpSolve gave poor results, run something else
         Prob.WarmStart = 0;
         Result = tomSolve(SolverDNLP,ProbDNLP);
         Inform=Result.Inform;
      end
   
      ExitFlag=Result.ExitFlag;
      x   = min(x_U,max(x_L,Result.x_k));
      v_k = Result.v_k;
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
            fprintf('Dual problem; Solution x:\n')
            xprint(x,'x:')
         end
         f_i_RP = Result.f_k;
         f_min(i) = f_i_RP;
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
      fprintf('! %s Inform %d ',SolverDNLP,Inform);
      if ~isinf(fDual) 
         fprintf('f_Dual %8.2f',fDual);
      end
      fprintf(' ExitFlag %d ',ExitFlag);
      fprintf('\n');
   end
   ENDTREE=1;
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
            if abs(fIPMax - nlp_f(xIPMax,ProbNLP)) < 1E-10 %CHANGED
            %if abs(fIPMax - c(1:n)'*xIPMax) < 1E-10
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
   fprintf('\n--- Total number of NLP iterations = %d\n',SumIter)
end
ResultMINLP.x_k  = xIPMax;
ResultMINLP.f_k  = fIPMax;
ResultMINLP.v_k  = vIPMax;
ResultMINLP.QP.B = BIPMax;
ResultMINLP.QP.y = yIPMax;
ResultMINLP.ExitFlag=1;
ResultMINLP.ExitText=ExitText(1);
ResultMINLP.Iter=Iter;
ResultMINLP.MinorIter=SumIter;
ResultMINLP.DualGap=Gap;
ResultMINLP=endSolve(Prob,ResultMINLP);

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
     Text = 'No solution found to NLP relaxation';
   case 5
     Text = 'Illegal x_0 found in NLP relaxation';
   case 6
     Text = 'User defined duality gap reached';
end
% MODIFICATION LOG
%
% 040331  hkh  Written, based on generalizing mipSolve
% 040907  med  Tested and refined
% 041023  hkh  KNAPSACK removed, safe guard input for fIP, xIP, VarWht
% 041023  hkh  Set PrintLM=0, no computation of Lagrange multipliers

