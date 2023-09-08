%
% rbfSolve is based on the RBF algorithms presented in 
% 1. Bjorkman, Holmstrom: Global Optimization of Costly Nonconvext Functions 
%    Using Radial Basis Functions, Optimization and Engineering 1,373-397, 2000
% 2. Hans-Martin Gutmann: A radial basis function method for global 
%    optimization, Journal of Global Optimization, 19,201:207, 2001.
% 3. Hans-Martin Gutmann: Radial Basis Function Methods for Global Optimization,
%    Ph.D. Thesis, Cambridge University, Cambridge, UK, September 2001.
%
% rbfSolve solves problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finite
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% Any of the x could be set as integer valued
%
% f(x) are assumed to be a costly function
% c(x) are assumed to be cheaply computed
% If some subset c_I(x) are very costly, create f(x) as a penalty function as
%
%     NEW f(x) = f(x) + beta' * c_I(x), beta positive penalties
%
% Calling syntax:
%
% function Result = rbfSolve(Prob, varargin)
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   USER.f    The routine to compute the function, given as a string, say RBFF
%   USER.c    The routine to compute the nonlinear constraint, say RBFC
%             A call to tomFiles.m or glcAssign.m sets these fields. 
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   WarmStart If true, >0, rbfSolve reads the output from the last run
%             from the mat-file cgoSave.mat, and continues from the last run.
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%   PriLevOpt Print Level 
%             0 = silent. 1 = Summary 2 = Printing each iteration
%             3 = Info about local / global solution 4 = Progress in x
%   user      User field used to send information to low-level functions
% --------------------------------------------
% optParam    Structure in Prob, Prob.optParam 
% ---------------------------------------------
%             Defines optimization parameters. Fields used: 
%  IterPrint  Print one information line each iteration, and the new x tried
%             Default IterPrint = 1.
%             fMinI means the best f(x) is infeasible
%             fMinF means the best f(x) is feasible (also integer feasible) 
%  MaxFunc    Maximal number of costly function evaluations, default 300
%             If WarmStart == 1 and MaxFunc <= nFunc (Number of f(x) used)
%             then MaxFunc = MaxFunc + nFunc
%  MaxIter    Maximal number of iterations used in the local optimization on
%             the response surface in each step. Default 1000, except for
%             pure IP problems, then max(GO.MaxFunc, MaxIter);
%  fGoal      Goal for function value, if empty or Inf not used
%  eps_f      Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%             See the output field maxTri.
%  bTol       Linear constraint tolerance
%  cTol       Nonlinear constraint tolerance
%
% ------------------
% Fields in Prob.CGO
% ------------------
%
% Percent      Type of strategy to get the initial sampled values:
%
%              Percent >= 100: User given initial points x, as a matrix of
%              points in CGO.X. Each column is one sampled point. 
%              If d = length(Prob.x_L), then size(X,1) = d, size(X,2) >= d+1
%              CGO.F should be defined as empty, or contain a vector of 
%              corresponding f(x) values. Any CGO.F value set as NaN will be 
%              computed by rbfSolve.
%
%              Percent > 0 (less than 100): Random strategy, the Percent value 
%              gives the percentage size of an ellipsoid around the so far 
%              sampled points that the new points are not allowed in.
%              Range 1%-50%. Recommended values 10% - 20%.
%
%              Percent == 0: Initial points is the corner points of the box
%              Generates too many points if the dimension is high.
%
%              Percent < 0: Latin hypercube space-filling design. The value
%              of abs(Percent) should in principle be the dimension. The call
%              made is X = daceInit(round(abs(Percent)),Prob.x_L,Prob.x_U);
%              Let k = abs(Percent), then the number of points are:
%              k      : 1  2  3  4  5  6  >6
%              Points : 21 21 33 41 51 65 65
%
%              Percent < -1000: Latin hypercube space-filling design, only keep
%              the points that fulfill the linear and nonlinear constraints. 
%              The algorithm will try up to M = abs(Percent)-1000 points,
%              stopping when it has got length(x_L)+1 feasible points
%
%              Percent == -999: Gutmann strategy: Initial points are x_L and
%              d points x_L  + (x_U(i) - x_L(i))*e_i, i=1,...,d
%
% X            If Percent >= 100, a matrix of initial x values
%              One column for every x value. size(X,2) >= dim(x)+1 needed
% F            If Percent >= 100, a vector of initial f(x) values.
%              If any element is set to NaN, rbfSolve will compute f(x)
% CX           If Percent >= 100, optionally a matrix of nonlinear 
%              constraint c(x) values. If nonempty, then 
%              size(CX,2) == size(X,2). If any element in a column i is set as
%              NaN, the vector c(x) = CX(:,i) will be recomputed
%
% RandState    If >=0, rand('state',RandState) is set to initialize the
%              pseudo-random generator
%              if < 0, rand('state',100*clock) is set to give a new set
%              of random values each run
%              Default RandState = 0
%
% idea         Type of search strategy on the surface, idea 1 (cycle of 6
%              points in target value fnStar) or idea 2 (cycle of 4 points in 
%              alpha, which implicitely gives the target value fnStar).
%              Default is 1.
% rbfType      Type of radial basis function.
%              1-Thin Plate Spline, 2-Cubic (default).
% SCALE        0-Original search space
%              1-Transform search space to unit cube (default).
% PLOT         0-No plotting (default), 1-Plotting sampled points.
% REPLACE      0-No replacement (default for constrained problems)
%              1-Large function values are replaced by the median, (uc default).
% LOCAL        0-No local searches after global search
%              If RBF surface is inaccurate, might be an advantage
%           
%              1-Local search from best points after global search. If equal
%              best function values, up to 20 local searches are done.
%
% N            Cycle length in idea 1 (default N=5 for fStarRule 1 and 2,
%              default N=1 for fStarRule 3) or idea 2 (always N=3)
%
% infStep      If =1, add search step with target value -inf first in cycle
%              Default 0.
% AddMP        If = 1, add the midpoint as extra point in the corner strategy
%              or the Gutmann strategy. Default 0.
% fStarRule    Global-Local search strategy in idea 1. N = cycle length
%              Define min_sn as local minimum on surface. fStar Target value
%              1: fStar = min_sn - ((N-(n-nInit))/N)^2*Delta_n (Default)
%              2: fStar = min_sn -  (N-(n-nInit))/N   *Delta_n
%              Strategy 1 and 2 depends on Delta_n estimate (see DeltaRule)
%              If infStep true, addition of -inf-step first in cycle
%              3: fStar = -inf-step, min_sn-k*0.1*|min_sn| k=N,...,0
%
%              Strategy names in Gutmanns thesis: III, II, I
%
% DeltaRule    1 = Skip large f(x) when computing f(x) interval Delta
%              0 = Use all points. Default 1.
% globalSolver Solver used for global optimization on the RBF surface
%              If the globalSolver is glcCluster, the fields
%              Prob.GO.maxFunc1, Prob.GO.maxFunc2 and Prob.GO.maxFunc3 are used
%              See the help for maxFunc1, maxFunc2, maxFunc3 in glcCluster
% localSolver  Solver used for local optimization on the RBF surface
%
% MaxCycle     Max number of cycles without progress before stopping, default 10
%
% varargin     Additional parameters to rbfSolve are sent to the costly f(x)
%
% ---------------------------------------------------------
% Fields in Prob.GO (Default values are set for all fields)
% ---------------------------------------------------------
%
% MaxFunc      Maximal number of function evaluations in each global search
% MaxIter      Maximal number of iterations in each global search
% maxFunc1     glcCluster parameter, max function evaluations 1st call
%              Only used if globalSolver is glcCluster, see help globalSolver
% maxFunc2     glcCluster parameter, max function evaluations 2nd call
%              Only used if globalSolver is glcCluster, see help globalSolver
% maxFunc3     glcCluster parameter, max sum of function evaluations in 
%              repeated 1st calls trying to get feasible
%              Only used if globalSolver is glcCluster, see help globalSolver
% localSolver  The local solver used by glcCluster
% DIRECT       DIRECT method used in glcCluster: glcSolve or glcFast(default)
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:  If IntVars is a scalar, then variables 1,...,IntVars are 
%             assumed to be integers. 
%             If empty, all variables are assumed non-integer
%             If length(IntVars) >1 ==> length(IntVars) == length(c) should hold
%             Then IntVars(i) ==1 ==> x(i) integer. IntVars(i) ==0 ==> x(i) real
%             If length(IntVars) < n, IntVars is assumed to be a set of indices.
%
% If MIP problems are solved then the only subsolvers working are glcCluster,
% and OQNLP (for both the local and global subproblem)
% e.g.
% Prob.CGO.globalSolver = 'oqnlp';
% Prob.CGO.localSolver  = 'oqnlp';
% will use the OQNLP solver, a license for Tomlab /OQNLP is needed
%
% For pure IP problems, only glcSolve and glcFast (default) may be used. Set
% Prob.CGO.globalSolver = 'glcSolve'; to use glcSolve, otherwise glcFast is used
%
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with the best points as columns, f(x_k) == f_k.
%  f_k      The best function value found so far
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag Always 0
%  Inform   0 = Normal termination
%           1 = Function value f(x) is less than fGoal
%           2 = Error in function value f(x), abs(f-fGoal) <= fTol, fGoal=0
%           3 = Relative Error in function value f(x) is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           4 = Failure in global sub solvers, same point sampled over and over
%               Try lower cTol and or bTol, or other subsolver
%           8 = No progress for MaxCycle*(N+1)+1 function evaluations 
%               (>MaxCycle cycles, input CGO.MaxCycle)
%           9 = Max CPU Time reached
%
% To make a warm start possible, rbfSolve saves the following information in
% the file cgoSave.mat:
% Name      Name of the problem
% O         Matrix with sampled points (in original space)
% X         Matrix with sampled points (in unit space if SCALE == 1)
% F         Vector with function values
% F_m       Vector with function values (replaced)
% nInit     Number of initial points >= d+1 (2^d if center points)
% Fpen      Vector with function values + additional penalty if infeasible
% fMinIdx   Index of the best point found
%
%
% USAGE:
%
% Let the name of the problem be "RBFF Test"
% The function RBFF is best written as
%     function f = RBFF(x, Prob) 
% Then any information, say u and W is easily sent to RBFF (and the constraint
% function RBFC, if present) using the Prob structure. See the example below. 
%
% Assume bounds x_L and x_U are initialized, as well as the linear constraint
% matrix A, lower and upper bounds, b_L and b_U on the linear constraints,
% and lower and upper bounds c_L and c_U on the nonlinear constraints
% (Put [] if all bounds are inf or -inf). Use the TOMLAB Quick format:
%
%    Prob   = glcAssign('RBFF',x_L,x_U,'RBFF Test',A,b_L,b_U,'RBFC',c_L,c_U);
%    Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%    % Default values are now set for PriLevOpt, and structure optParam
%    % To change a value, examples are shown on the two next lines
%    Prob.optParam.MaxFunc = 300;  % Change max number of function evaluations 
%    Prob.optParam.MaxIter = 2000; % Change local iterations to 2000
%    Prob.GO.MaxFunc       = 1000; % 1000 global function values
%
% Direct solver call:
%    Result = rbfSolve(Prob);
%    PrintResult(Result);
%           
% Driver call, including printing with level 2:
%      Result = tomRun('rbfSolve',Prob,2);
%           
% The user function RBFF is written as
%
%    function f = RBFF(x, Prob) 
%    u = Prob.user.u; W = Prob.user.W;
%    f = "some function of x, u and W"
%
% It is also possible to use the function format
%    function f = RBFF(x) 
% but then any additional parameters must be sent as global variables.
%
% The user function RBFC, computing the nonlinear constraints, is written as
%
%    function c = RBFC(x, Prob) 
%    u = Prob.user.u; W = Prob.user.W;
%    c = "some vector function of x, V and W"
%
% Note! If RBFF has the 2nd input argument Prob, also RBFC must have that.
%
% To make a restart, just set the restart flag, and call rbfSolve once again:
%
%    Prob.WarmStart = 1;
%    Result = tomRun('rbfSolve',Prob,2);
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 2.5.0$
% Written Jan 12, 2000. Last modified Nov 30, 2004.
%


function Result = rbfSolve(Prob, varargin)

if nargin < 1
   error('rbfSolve needs one parameter, the structure Prob');
end

time  = fix(clock);

DebugPriLev = 0;  % PriLev in subopt, i.e. gnProb.PriLevOpt, snProb.PriLevOpt

solvType=checkType('glc');

Prob=ProbCheck(Prob,'rbfSolve',solvType);

Prob = iniSolve(Prob,solvType,0,0);

MaxCPU = Prob.MaxCPU;

if isempty(Prob.x_L) | isempty(Prob.x_U)
   disp('rbfSolve requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'rbfSolve requires both lower and upper variable bounds';
   Result=endSolve(Prob,Result);
   return;
end

% Pick up input variables from Prob structure
PriLev    = Prob.PriLevOpt;          % Print level
f_Low     = Prob.f_Low;              % Lower bound on f
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds
MaxIter   = Prob.optParam.MaxIter;   % Iterations in local suboptimization
MaxFunc   = Prob.optParam.MaxFunc;   % Number of costly function evaluations
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration

fTol      = Prob.optParam.eps_f;     % Relative convergence tolerance in f(x)
fGoal     = Prob.optParam.fGoal;     % Goal f(x) for the optimization
epsRank   = Prob.optParam.eps_Rank;  % Rank tolerance in qr-decomposition
bTol      = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol      = Prob.optParam.cTol;      % Constraint feasibility tolerance

%xTol      = Prob.optParam.eps_x;     % Tolerance for rectangle sizes
%                                     % (scaled to (0,1) )

eps_sn    = 1E-6; % Tolerance for min_sn

x_D   = x_U - x_L;
d     = length(x_L);  % Problem dimension

% fnvec = []; % Collect all fnStar values

% Safeguard
if isempty(MaxFunc), MaxFunc = 300; end
if MaxFunc < 0
   MaxFunc = 300;
end
if isempty(MaxIter), MaxIter = 1000; end
if MaxIter < 0
   MaxIter = 1000;
end
if isempty(IterPrint), IterPrint = 1; end

if isempty(Prob.CGO)
   idea      = []; rbfType   = []; SCALE        = []; PLOT        = [];
   REPLACE   = []; Percent   = []; globalSolver = []; localSolver = [];
   LOCAL     = []; N         = []; infStep      = []; AddMP       = []; 
   fStarRule = []; DeltaRule = []; RandState    = []; MaxCycle    = [];

else

   if isfield(Prob.CGO,'idea') 
      idea = Prob.CGO.idea;
   else
      idea = [];
   end
   if isfield(Prob.CGO,'rbfType')
      rbfType = Prob.CGO.rbfType;
   else
      rbfType = [];
   end
   if isfield(Prob.CGO,'SCALE')
      SCALE = Prob.CGO.SCALE;
   else
      SCALE = [];
   end
   if isfield(Prob.CGO,'PLOT')
      PLOT = Prob.CGO.PLOT;
   else
      PLOT = [];
   end
   if isfield(Prob.CGO,'REPLACE')
      REPLACE = Prob.CGO.REPLACE;
   else
      REPLACE = [];
   end
   if isfield(Prob.CGO,'LOCAL')
      LOCAL = Prob.CGO.LOCAL;
   else
      LOCAL = [];
   end
   if isfield(Prob.CGO,'N')
      N = Prob.CGO.N;
   else
      N = [];
   end
   if isfield(Prob.CGO,'infStep')
      infStep = Prob.CGO.infStep;
   else
      infStep = [];
   end
   if isfield(Prob.CGO,'AddMP')
      AddMP = Prob.CGO.AddMP;
   else
      AddMP = [];
   end
   if isfield(Prob.CGO,'fStarRule')
      fStarRule = Prob.CGO.fStarRule;
   else
      fStarRule = [];
   end
   if isfield(Prob.CGO,'DeltaRule')
      DeltaRule = Prob.CGO.DeltaRule;
   else
      DeltaRule = [];
   end
   if isfield(Prob.CGO,'Percent')
      Percent = Prob.CGO.Percent;
   else
      Percent = [];
   end
   if isfield(Prob.CGO,'RandState')
      RandState = Prob.CGO.RandState;
   else
      RandState = [];
   end
   if isfield(Prob.CGO,'globalSolver')
      globalSolver = deblank(Prob.CGO.globalSolver);
   else
      globalSolver = [];
   end
   if isfield(Prob.CGO,'localSolver')
      localSolver = deblank(Prob.CGO.localSolver);
   else
      localSolver = [];
   end
   if isfield(Prob.CGO,'MaxCycle') 
      MaxCycle = Prob.CGO.MaxCycle;
   else
      MaxCycle = [];
   end
end
if isempty(Prob.GO)
   GOMaxFunc = []; GOMaxIter     = []; maxFunc1 = [];  maxFunc2 = [];
   maxFunc3  = []; GOlocalSolver = []; DIRECT   = [];
else
   if isfield(Prob.GO,'MaxFunc')
      GOMaxFunc = Prob.GO.MaxFunc;
   else
      GOMaxFunc = [];
   end
   if isfield(Prob.GO,'MaxIter')
      GOMaxIter = Prob.GO.MaxIter;
   else
      GOMaxIter = [];
   end
   if isfield(Prob.GO,'maxFunc1')
      maxFunc1 = Prob.GO.maxFunc1;
   else
      maxFunc1 = [];
   end
   if isfield(Prob.GO,'maxFunc2')
      maxFunc2 = Prob.GO.maxFunc2;
   else
      maxFunc2 = [];
   end
   if isfield(Prob.GO,'maxFunc3')
      maxFunc3 = Prob.GO.maxFunc3;
   else
      maxFunc3 = [];
   end
   if isfield(Prob.GO,'localSolver')
      GOlocalSolver = Prob.GO.localSolver;
   else
      GOlocalSolver = [];
   end
   if isfield(Prob.GO,'DIRECT')
      DIRECT = Prob.GO.DIRECT;
   else
      DIRECT = [];
   end
end
if ~isempty(Prob.MIP)
   IntVars = Prob.MIP.IntVars;
   if ~isempty(IntVars)
      if length(IntVars)==1
         IntVars=[1:min(d,max(0,floor(IntVars)))];
      elseif length(IntVars) < d
         IntVars=IntVars;
      elseif any(IntVars > 1)
         IntVars=1:d;
      else % logical variables
         IntVars=find(IntVars > 0);
      end
      if length(IntVars) == d
         % Pure IP problem
         if strcmpi(globalSolver,'glcSolve')
            globalSolver = 'glcSolve';
            localSolver  = 'glcSolve';
         elseif strcmpi(globalSolver,'oqnlp')
            globalSolver = 'oqnlp';
            localSolver  = 'oqnlp';
         else
            globalSolver = 'glcFast';
            localSolver  = 'glcFast';
         end
         LOCAL        = 0;
         if isempty(GOMaxFunc)
            GOMaxFunc    = 5000;
         end
         if isempty(GOMaxIter)
            GOMaxIter    = 5000;
         end
         if isempty(MaxIter)
            MaxIter      = max(GOMaxFunc,GOMaxIter);
         else
            MaxIter      = max(MaxIter,GOMaxFunc);
         end
      elseif strcmpi(globalSolver,'oqnlp')
         if ~strcmpi(localSolver,'oqnlp') & ~strcmpi(localSolver,'glcCluster')
            localSolver  = 'oqnlp';
         end
         LOCAL = 0;
      elseif strcmpi(globalSolver,'glcFast') ...
         if isempty(localSolver)
            localSolver  = 'glcFast';
         end
         LOCAL = 0;
      elseif strcmpi(globalSolver,'glcSolve')
         if isempty(localSolver)
            localSolver  = 'glcSolve';
         end
         LOCAL = 0;
      elseif strcmpi(globalSolver,'minlpSolve')
         if isempty(localSolver)
            localSolver  = 'minlpSolve';
         end
         LOCAL = 0;
      elseif strcmpi(globalSolver,'minlpBB')
         if isempty(localSolver)
            localSolver  = 'minlpBB';
         end
         LOCAL = 0;
      else
         globalSolver = 'glcCluster';
         if isempty(GOlocalSolver)
            if isempty(localSolver)
               GOlocalSolver = GetSolver('con',1,0);
            else
               GOlocalSolver  = localSolver;
            end
         end
         localSolver  = 'glcCluster';
      end
      SCALE = 0;
      % Safe guard bounds onto integer values
      x_L(IntVars)   = ceil(x_L(IntVars)); 
      x_U(IntVars)   = floor(x_U(IntVars)); 
      x_D            = x_U - x_L;
   end
else
   IntVars = [];
end

nFunc     = 0;
nCon      = 0;        % Count number of calls to constraint routine, nlp_c
convflag  = 0;
nmax      = [];       % Only used in idea 1
CX        = [];       % Initial set of nonlinear constraint values
DEBUG     = 0;

Result                 = ResultDef(Prob);
Result.Solver          = 'rbfSolve';
Result.SolverAlgorithm = 'Radial Basis Function Interpolation';

% Default strategies for RBF input
if isempty(RandState),    RandState = 0; end
if isempty(idea),         idea = 1; end
if isempty(rbfType),      rbfType = 2; end
if isempty(SCALE),        SCALE = 1; end
if isempty(PLOT),         PLOT = 0; end
if isempty(LOCAL),        LOCAL = 1; end
if isempty(infStep),      infStep = 0; end
if isempty(AddMP),        AddMP = 0; end
if isempty(fStarRule),    fStarRule = 1; end
if isempty(DeltaRule),    DeltaRule = 1; end
if isempty(globalSolver), globalSolver = 'glcFast'; end
if isempty(localSolver),  localSolver = GetSolver('con',1,0); end
if isempty(MaxCycle),     MaxCycle = 10; end

%if isempty(FACTORIZE), FACTORIZE = 1; end


if idea == 1      % Cycle length for direct fnStar strategy
   if fStarRule == 3
      infStep = 1; % Must have infStep in strategy 3, (Gutmann I)
      if isempty(N), N = 1; end
   else
      if isempty(N), N = 5; end
   end
elseif idea == 2  % Cycle length for indirect alpha strategy
   %if isempty(N), N = 3; end
   N = 3; % Always N = 3
end


AmpFac = 1.0;

% Set pseudo random generator
if RandState >=0 
   rand('state',RandState);
else
   rand('state',100*clock);
end

dLin = size(Prob.A,1);
if dLin > 0
   d1 = length(Prob.b_L);
   % Must adjust bounds to have equal length
   if d1 < dLin
      Prob.b_L = [Prob.b_L;-inf*ones(dLin-d1,1)];
   end
   d1 = length(Prob.b_U);
   if d1 < dLin
      Prob.b_U = [Prob.b_U;inf*ones(dLin-d1,1)];
   end
end
dCon = max(length(Prob.c_L),length(Prob.c_U));
if dCon > 0
   d1 = length(Prob.c_L);
   % Must adjust bounds to have equal length, to add extra EGO constraint
   if d1 < dCon
      Prob.c_L = [Prob.c_L;-inf*ones(dCon-d1,1)];
   end
   d1 = length(Prob.c_U);
   if d1 < dCon
      Prob.c_U = [Prob.c_U;inf*ones(dCon-d1,1)];
   end
end
% Default values for GO structure
if dLin + dCon > 0
   if isempty(REPLACE),      REPLACE = 0; end
   if isempty(Percent),      Percent = -5000; end
   if strcmpi(globalSolver,'glcCluster')
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
      if isempty(maxFunc1),     maxFunc1  = 1000; end
      if isempty(maxFunc2),     maxFunc2  = 2000; end
      if isempty(maxFunc3),     maxFunc3  = 3000; end
   elseif strcmpi(globalSolver,'glcFast')
      if isempty(GOMaxFunc),    GOMaxFunc = max(3000,300*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,300*d); end
   elseif strcmpi(globalSolver,'glcSolve')
      if isempty(GOMaxFunc),    GOMaxFunc = max(1500,150*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(1500,150*d); end
   else
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(1000,1000*d); end
      if isempty(maxFunc1),     maxFunc1  = 500; end
      if isempty(maxFunc2),     maxFunc2  = 500; end
      if isempty(maxFunc3),     maxFunc3  = 500; end
   end
   if strcmpi(localSolver,'glcCluster')
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
      if isempty(maxFunc1),     maxFunc1  = 1000; end
      if isempty(maxFunc2),     maxFunc2  = 2000; end
      if isempty(maxFunc3),     maxFunc3  = 3000; end
      MaxIter = max(MaxIter,GOMaxIter);
   end
else
   if isempty(REPLACE),      REPLACE = 1; end
   if isempty(Percent),      Percent = -d; end
   if strcmpi(globalSolver,'glcCluster')
      if isempty(GOMaxFunc),    GOMaxFunc = max(5000,500*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,250*d); end
      if isempty(maxFunc1),     maxFunc1  = 500; end
      if isempty(maxFunc2),     maxFunc2  = 500; end
      if isempty(maxFunc3),     maxFunc3  = 500; end
   elseif strcmpi(globalSolver,'glcFast')
      if isempty(GOMaxFunc),    GOMaxFunc = max(2000,200*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(2000,200*d); end
   elseif strcmpi(globalSolver,'glcSolve')
      if isempty(GOMaxFunc),    GOMaxFunc = max(1000,100*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(1000,100*d); end
   else
      if isempty(GOMaxFunc),    GOMaxFunc = max(5000,500*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,250*d); end
      if isempty(maxFunc1),     maxFunc1  = 500; end
      if isempty(maxFunc2),     maxFunc2  = 500; end
      if isempty(maxFunc3),     maxFunc3  = 500; end
   end
end
% NOTE! Hard coded constants determining current cycle strategy.
% Should be further tested and developed to see which one is best.
% The old rbfSolve
OldStrat = 1; % If false, use a new Idea 1 median strategy
if dLin + dCon > 0
   OldAlpha = 0; % If false, use a new Idea 2 alpha  strategy
else
   OldAlpha = 1; % If false, use a new Idea 2 alpha  strategy
end

% Scaling parameters
if SCALE
   LOWER = 0;
   UPPER = 1;
   x_LL = LOWER*ones(d,1);
   x_UU = UPPER*ones(d,1);
   x_DD = x_UU-x_LL;
else
   x_LL = x_L;
   x_UU = x_U;
   x_DD = x_D;
end

%
%  STEP 1, Initialization
%

if Prob.WarmStart
   % Restart with values from previous run.

   load('cgoSave.mat','Name','O','F','X','F_m','nInit','Fpen','fMinIdx');

   Name1 = deblank(Prob.Name);   % Name for the problem

   if strcmp(Name1,Name)
      [d n] = size(X);
      nFunc = n; % Total count of function evaluations
   
      if PriLev > 0
         fprintf('\n Restarting with %d sampled points from previous run\n',n);
      end
   else
      nFunc = 0; % Total count of function evaluations
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files cgoSave.mat?\n');
      end
   end

   % TOMSOL INIT, send F_m to Fortran, not F
   control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
   if control < 0
       fprintf('Initial interpolation failed');
       tomsol(25) % Deallocates memory
       return     %Something is really wrong
   end
   
   % Add extra function evaluations 
   if MaxFunc <= nFunc
      MaxFunc = MaxFunc + nFunc;
   end
   
end

if ~Prob.WarmStart
   Name = deblank(Prob.Name);  % Problem name
      
   % Three DIFFERENT INITIALIZATION STRATEGIES, plus user defined

   F = [];
   if Percent >= 100
      % User defined points
      X = Prob.CGO.X;
      if size(X,1) ~= d | size(X,2) < d+1
         fprintf('Size of user given X matrix is wrong, size %d %d\n',size(X));
         fprintf('1st dimension must be   %d\n',d);
         fprintf('2nd dimension must be > %d\n',d+1);
         error('Illegal input to rbfSolve')
      end
      F = Prob.CGO.F(:);
      if SCALE
         O = X;
         D = x_U-x_L;
         D(D==0) = 1; 
         for i=1:size(X,2)
             x = X(:,i);
             x = max(x_L,min(x_U,x));
             X(:,i) = (x-x_L)./D;
         end
         %X = tomsol(9, zeros(d,1), O ,ones(d,1)); 
      else
         O = X;
      end
      if isfield(Prob.CGO,'CX')
         CX = Prob.CGO.CX;
         if ~isempty(CX) 
            if size(CX,1) ~= dCon & size(CX,2) ~= size(X,2)
              'size of CX'
               size(CX)
               error('rbfSolve: Illegal size of Prob.CGO.CX')
            end
         end
      end
   elseif Percent > 0
      Percent = max(1,min(50,Percent));
      if PriLev > 0 | IterPrint
         fprintf('Random. Percentage of ellipsoid %f\n',Percent);
      end
      ix = x_LL~=x_UU;
      ix(IntVars) = 0;
      iV = find(ix);
      XX = randomtest(x_LL(iV),x_UU(iV),Percent);
      X  = zeros(d,max(size(XX,2),d+1));
      X(iV,1:size(XX,2)) = XX;
      for i=size(XX,2)+1:size(X,2)
          X(:,i) = x_LL;
          X(iV,i) = 0.7*x_LL(iV)+0.3*x_UU(iV);
      end
      if ~isempty(IntVars)
         X = sampleInts(X,x_LL,x_UU,IntVars);
      end

   elseif Percent == 0
      % Setup initial points (the corners of the box)
      if PriLev > 0 | IterPrint
         fprintf('Corner points %d\n',2^d);
      end
      X = corners(x_LL,x_UU);
      if AddMP
         if isempty(IntVars)
            X = [X,0.5*(x_LL+x_UU)];
         elseif length(IntVars) < d
            X = [X,0.5*(x_LL+x_UU)];
            n = size(X,2);
            X(:,n) = sampleInts(X(:,n),x_LL,x_UU,IntVars);
            % Avoid mid point if all variables are integer
            if any(all(X(:,1:n-1)==X(:,n)*ones(1,n-1)))
               n = n-1;
               X = X(:,1:n);
            end
         end
      end
   elseif Percent == -999
      % Setup initial points, Gutmann coordinate strategy
      if PriLev > 0 | IterPrint
         fprintf('Gutmann coordinate strategy with %d points\n',d+1);
      end
      X = gutmann(x_LL,x_UU);
      if AddMP
         if isempty(IntVars)
            X = [X,0.5*(x_LL+x_UU)];
         elseif length(IntVars) < d
            X = [X,0.5*(x_LL+x_UU)];
            n = size(X,2);
            X(:,n) = sampleInts(X(:,n),x_LL,x_UU,IntVars);
            % Avoid mid point if all variables are integer
            % Avoid mid point if all variables are integer
            if any(all(X(:,1:n-1)==X(:,n)*ones(1,n-1)))
               n = n-1;
               X = X(:,1:n);
            end
         end
      end
   elseif Percent > -999
      dim = round(abs(Percent));
      if PriLev > 0 | IterPrint
         fprintf('Latin Hypercube, dimension %d\n',dim);
      end
      if isempty(IntVars)
         X  = daceInit(dim,x_LL,x_UU);
      elseif length(IntVars) == d
         % Only integer variables
         n  = max(d+1,dim);
         X  = double(zeros(d,n));
         X  = sampleInts(X,x_LL,x_UU,IntVars);
      else
         ix = ones(d,1);
         ix(IntVars) = 0;
         xV = find(ix);
         XX = daceInit(dim,x_LL(xV),x_UU(xV));
         n  = size(XX,2);
         X  = zeros(d,n);
         X(xV,:) = XX;
         X  = sampleInts(X,x_LL,x_UU,IntVars);
      end
   else
      M = max(d+1,abs(Percent) - 1000);
      if PriLev > 0 | IterPrint
         fprintf('Latin Hypercube, ');
         fprintf('try <= %d, to get %d feasible points\n',M,d+1);
      end
      if isempty(IntVars)
         X  = daceInit(d+1,x_LL,x_UU,M);
         n  = size(X,2);
      elseif length(IntVars) == d
         % Only integer variables
         n  = max(d+1,M);
         X  = double(zeros(d,n));
         X  = sampleInts(X,x_LL,x_UU,IntVars);
      else
         ix = ones(d,1);
         ix(IntVars) = 0;
         xV = find(ix);
         XX = daceInit(d+1,x_LL(xV),x_UU(xV),M);
         n  = size(XX,2);
         X  = zeros(d,n);
         X(xV,:) = XX;
         X  = sampleInts(X,x_LL,x_UU,IntVars);
      end
      if SCALE
         O = tomsol(9, x_L, X ,x_D); 
      else
         O = X;
      end
      ixOK = zeros(n,1);
      nOK  = 0;
      PNTS = d+1;
      F    = NaN*ones(PNTS,1);
      cErrBest=Inf;
      ixBest = Inf;
      for i = 1:n
          if nOK >= PNTS, break; end
          cErr = 0;
          if dLin > 0
             L = Prob.A*O(:,i);
             cErr = cErr+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
           end
           if dCon > 0  & cErr == 0
              C = nlp_c(O(:,i), Prob, varargin{:});
              nCon = nCon + 1;
              cErr = cErr+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
           end
           if cErr == 0 
              nOK = nOK + 1; 
              ixOK(i) = 1;
              F(nOK) = nlp_f(O(:,i), Prob, varargin{:});
           elseif cErr < cErrBest
              cErrBest = cErr;
              ixBest   = i;
           end
      end
      X = X(:,find(ixOK));
      
      if size(X,2) < PNTS
         % First add best infeasible point
         if ixOK(ixBest) == 0
            ixOK(ixBest) = 2;
            X = [X,O(:,ixBest)];
         end
         % Then just take the first sampled infeasible points
         if size(X,2) < PNTS
            for i = 1:n
                if ixOK(i) == 0
                   X = [X,O(:,i)];
                end
                if size(X,2) >= PNTS, break; end
            end
         end
      end
      if PriLev > 0 | IterPrint
         iV = find(ixOK == 1);
         fprintf('Found %d out of %d ',length(iV),PNTS);
         fprintf('feasible points at trials ');
         fprintf('%d ',iV);
         fprintf('\n');
         if length(iV) < PNTS
            fprintf('Deviation for best non feasible point %f\n',cErrBest);
         end
      end
      %X(:,1:length(iV))
   end

   if SCALE
      O = tomsol(9, x_L, X ,x_D); 
   else
      O = X;
   end
   % Set dimension d and number of points n
   [d n] = size(X);
       
   % Set initial number of points nInit to n
   nInit = n;

   % Compute objective function value in each undefined point x in X
   % Compute constraint values for all points
   if isempty(F)
      F = NaN*ones(n,1);
   end
   Fpen = F;
   for i = 1:n
       % O(;.i) is vector in original space
       if isnan(F(i))
          F(i) = nlp_f(O(:,i), Prob, varargin{:});
       end
       Fpen(i) = F(i);
       if dLin > 0
          L = Prob.A*O(:,i);
          Fpen(i) = Fpen(i)+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
       end
       if dCon > 0 
          if isempty(CX) 
             C = nlp_c(O(:,i), Prob, varargin{:});
             nCon = nCon + 1;
          else
             C = CX(:,i);
             if any(isnan(C))
                C = nlp_c(O(:,i), Prob, varargin{:});
                nCon = nCon + 1;
             end
          end
          Fpen(i) = Fpen(i)+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
       end
   end
   
   % Replace large function values by median(F)

   if REPLACE
      F_m = min(median(F),F);
   else
      F_m = F;
   end
 
   % TOMSOL INIT, send F_m to Fortran, not F
   control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
   if control < 0
       fprintf('Initial interpolation failed');
       tomsol(25) % Deallocates memory
       return     %Something is really wrong
   end
   nFunc = n; % Number of function evaluations.
end

dc = Prob.USER.dc;

if SCALE > 0 & (dCon > 0 | dLin > 0)
   if dLin > 0
      A   = Prob.A*diag(x_D);
      b_L = Prob.b_L - Prob.A*x_L;
      b_U = Prob.b_U - Prob.A*x_L;
   else
      A   = [];
      b_L = [];
      b_U = [];
   end
   if dCon > 0
      snProb             = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn', ...
                           [],[],[], A,b_L,b_U, 'rbf_c','rbf_dc','rbf_d2c',...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
      snProb.xD          = x_D;
      snProb.xL          = x_L;
      snProb.cNargin     = xnargin(Prob.USER.c);
      if isempty(dc)
         snProb.dcNargin = 0;
      else
         snProb.dcNargin = xnargin(dc);
      end
      snProb.c           = Prob.USER.c;
      snProb.dc          = Prob.USER.dc;
      snProb.d2c         = Prob.USER.d2c;
   else
      snProb             = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn', ...
                           [], [], [], A, b_L, b_U); 
   end
   snProb.dDim           = d;
   snProb.SCALE          = SCALE;

elseif dCon == 0 & dLin == 0
   snProb                = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn'); 
else
   snProb                = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn', ...
                           [], [], [], Prob.A, Prob.b_L, Prob.b_U, ...
                           Prob.USER.c, Prob.USER.dc, Prob.USER.d2c, ...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
end

optParam                  = optParamDef(localSolver,solvType,d,dCon,dCon+dLin);
snProb.optParam           = optParam;
snProb.optParam.MaxIter   = MaxIter;
if strcmpi(localSolver,'glcCluster')
   snProb.optParam.MaxFunc   = GOMaxFunc;
else
   snProb.optParam.MaxFunc   = MaxIter*max(d,10)/10;  
end
snProb.optParam.IterPrint = DebugPriLev > 0;
snProb.optParam.cTol      = cTol;
snProb.optParam.bTol      = bTol;
snProb.GradTolg           = Prob.GradTolg;
snProb.GradTolH           = Prob.GradTolH;
snProb.GradTolJ           = Prob.GradTolJ;
snProb.SOL                = Prob.SOL;
snProb.GO                 = Prob.GO;
if isfield(Prob,'OQNLP')
   snProb.OQNLP           = Prob.OQNLP;
end
if isfield(Prob,'KNITRO')
   snProb.KNITRO          = Prob.KNITRO;
end

snProb.PriLevOpt          = DebugPriLev; 

% Send user info
if isfield(Prob,'user')
   snProb.user            = Prob.user;
end
snProb.MIP                = Prob.MIP;
snProb.uP                 = Prob.uP;
snProb.P                  = Prob.P;
snProb.mLin               = dLin;
snProb.mNonLin            = dCon;


if SCALE > 0 & (dCon > 0 | dLin > 0)
   if dLin > 0
      A   = Prob.A*diag(x_D);
      b_L = Prob.b_L - Prob.A*x_L;
      b_U = Prob.b_U - Prob.A*x_L;
   else
      A   = [];
      b_L = [];
      b_U = [];
   end
   if dCon > 0
      gnProb             = glcAssign('gn_f',x_LL,x_UU,'RBFgn',A,...
                           b_L,b_U,'rbf_c',Prob.c_L,Prob.c_U);
      gnProb.xD          = x_D;
      gnProb.xL          = x_L;
      gnProb.USER.dc     = 'rbf_dc';
      gnProb.c           = Prob.USER.c;
      gnProb.dc          = Prob.USER.dc;
      gnProb.d2c         = Prob.USER.d2c;
      gnProb.cNargin     = xnargin(Prob.USER.c);
      if isempty(dc)
         gnProb.dcNargin = 0;
      else
         gnProb.dcNargin = xnargin(dc);
      end
      gnProb.dDim        = d;
      gnProb.SCALE       = SCALE;
   else
      gnProb             = glcAssign('gn_f',x_LL,x_UU,'RBFgn',A, b_L,b_U);
   end

elseif dCon == 0 & dLin == 0
   gnProb                = glcAssign('gn_f',x_LL,x_UU,'RBFgn');
else
   gnProb                = glcAssign('gn_f',x_LL,x_UU,'RBFgn',Prob.A,...
                           Prob.b_L,Prob.b_U,Prob.USER.c,Prob.c_L,Prob.c_U);
   gnProb.USER.dc        = Prob.USER.dc;
end

% Send user info
if isfield(Prob,'user')
   gnProb.user           = Prob.user;
end
gnProb.MIP               = Prob.MIP;
gnProb.uP                = Prob.uP;
gnProb.P                 = Prob.P;
gnProb.mLin              = dLin;
gnProb.mNonLin           = dCon;

optParam                 = optParamDef(globalSolver,solvType,d,dCon,dCon+dLin);
gnProb.optParam          = optParam;
gnProb.optParam.MaxIter  = GOMaxIter; 
gnProb.optParam.MaxFunc  = GOMaxFunc;  
gnProb.optParam.eps_Rank = epsRank;
gnProb.optParam.IterPrint = DebugPriLev > 0;
gnProb.optParam.cTol     = cTol;
gnProb.optParam.bTol     = bTol;
gnProb.GradTolg          = Prob.GradTolg;
gnProb.GradTolH          = Prob.GradTolH;
gnProb.GradTolJ          = Prob.GradTolJ;
gnProb.GO                = Prob.GO;
gnProb.SOL               = Prob.SOL;
if isfield(Prob,'OQNLP')
   gnProb.OQNLP          = Prob.OQNLP;
end
if isfield(Prob,'KNITRO')
   gnProb.KNITRO         = Prob.KNITRO;
end
gnProb.CGO.idea          = idea;
gnProb.CGO.epsRank       = epsRank;
gnProb.PriLevOpt         = DebugPriLev; 

% Take advantage of any derivatives in local search

if isempty(dc)
   DerLvl_sn             = 1;
   DerLvl_gn             = 0;
   cDiff                 = 6;
else
   DerLvl_sn             = 3;
   DerLvl_gn             = 2;
   cDiff                 = 0;
end

snProb.GO.maxFunc1    = maxFunc1;
snProb.GO.maxFunc2    = maxFunc2;
snProb.GO.maxFunc3    = maxFunc3;
snProb.GO.DIRECT      = DIRECT;

switch lower(localSolver)
 case {'snopt','npsol','nlssol','minos'}
   gnProb.NumDiff        = 6;
   gnProb.SOL.optPar(39) = DerLvl_gn;
   gnProb.ConsDiff       = cDiff;
   snProb.SOL.optPar(39) = DerLvl_sn;
   snProb.ConsDiff       = cDiff;
   %snProb.optPar(10)     = 1E-8;
   snProb.optPar(12)     = 1E-8; % Minor optimality tolerance
 case {'glccluster'}
   if isempty(GOlocalSolver)
      snProb.GO.localSolver = localSolver;
   else
      snProb.GO.localSolver = GOlocalSolver;
   end
 case {'glcFast','glbFast'}
   snProb.NumDiff        = 0;
   snProb.ConsDiff       = 0;
 otherwise
   gnProb.NumDiff        = 1;
   gnProb.ConsDiff       = cDiff > 0;
   snProb.ConsDiff       = cDiff > 0;
end
gnProb.GO.maxFunc1    = maxFunc1;
gnProb.GO.maxFunc2    = maxFunc2;
gnProb.GO.maxFunc3    = maxFunc3;
gnProb.GO.DIRECT      = DIRECT;
switch lower(globalSolver)
 case {'glccluster'}
   LOCAL = 0;
   if isempty(GOlocalSolver)
      gnProb.GO.localSolver = localSolver;
   else
      gnProb.GO.localSolver = GOlocalSolver;
   end
 case {'glcFast','glbFast'}
   gnProb.NumDiff        = 0;
   gnProb.ConsDiff       = 0;
 otherwise
end

Iter       = 0;
minDistEps = 1E-5;
Update     = 1;
fnStar     = Inf;
alphaXXX   = 0;

z = Fpen;
% Set infeasible points as infinity before check on minimum
z(Fpen-F >= 1E-14)=Inf;

% Set integer infeasible points as infinity before check on minimum
if ~isempty(IntVars)
   XX=X(IntVars,:);
   if length(IntVars) == 1
      v = XX~=round(XX);
   else
      v = ~all(XX==round(XX));
   end
   z(v) = Inf;
end

[fMin,fIdx] = min(z);

if isinf(fMin)
   % No feasible point found
   Feasible = 0;
   % Take the best infeasible point
   [fMin,fIdx] = min(Fpen);
else
   Feasible = 1;
end

% Best point found in unit space
x_min = X(:,fIdx(1)); 

% Best point found in original space
if SCALE
   O_min = tomsol(9, x_L, x_min, x_D); 
else
   O_min = x_min;
end

FLOWER    = fIdx(1);
fMinIter  = 0;
fDiff     = Inf;
fDiffOld  = NaN;
NOUPDATE  = 0;
SAME1     = 0;
SAME2     = 0;
O_pre     = inf*ones(d,1);
alpha     = NaN;

if isempty(IntVars)
   snProb = ProbCheck(snProb,localSolver,checkType('con'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('glc'));
else
   snProb = ProbCheck(snProb,localSolver,checkType('minlp'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('minlp'));
end

% *********************************************************
% *************** START MAIN ITERATION LOOP ***************
% *********************************************************

% n     = Number of points used in the interpolation
% nFunc = Total number of costly function evaluations
% nInit = Number of points used for startup, now or before
% Iter  = Number of fresh sampled points 
% modN  = Cycle step

modN = -1-infStep;

if PriLev > 1 | IterPrint
   fprintf('Iter %3d n %3d ', 0, n);
      fprintf('nFunc %3d ', nFunc);
      tt = time([3 2 4 5]);
      if tt(4) < 10
         fprintf('%d/%d %d:0%1d  ', tt);
      else
         fprintf('%d/%d %d:%2d ', tt);
      end
      fprintf('Cycle %2d ', modN);
      fprintf('               ');
      if ~isempty(fGoal) & ~isinf(fGoal)
         fprintf(' fGoal %8.5f', fGoal);
      end
      if Feasible
         fprintf(' fMinF %11.8f', fMin);
      else
         fprintf(' fMinI %11.8f', fMin);
      end
      %if TRANSFORM ~= 0
      %   fprintf('yMin %8.5f yNew %8.5f ', yMin, yNew);
      %end
      %fprintf('fE %8.5f ',fk);
      %fprintf('I%%');
      %fprintf(' %6.4f ',100*abs(EGOResult.f_k/yMin));
      fprintf('\n');
      xprint(O_min,'xMin:',' %12.8f',8)
   end

% -------------- Result saving -------------------
%saveInit;

convflag = 0;
cpumax   = 0;
progress = 1;
TIME0    = Prob.TIME0;
FLOW     = nFunc;

while nFunc < MaxFunc & convflag == 0
   
   if nFunc-FLOW > MaxCycle*(N+1), progress = 0; break; end
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end

   % [c,z1,z2,mode]=lsei(A,[F_m;zeros(d+1,1)],[],[],[],[],0,0,1);
   % lambda = c(1:n);
   % b = c(n+1:n+d);
   % a = c(n+d+1);
   %[c,z1,z2,mode]=lsei(A,[F_m;zeros(d+1,1)],[],[],[],[],0,0,1);
   %lambda = c(1:n);
   %b = c(n+1:n+d);
   %a = c(n+d+1);
   %RBF.lambda = lambda;
   %RBF.a      = a;
   %RBF.b      = b;
   %RBF.rbfType= rbfType;
  
   % Set parameters in global structure CGO
   % gnProb.CGO.X          = X;
   % gnProb.CGO.n          = n;
   % gnProb.CGO.d          = d;

   gnProb.CGO.rbfType    = rbfType;

   % gnProb.CGO.minDistEps = minDistEps;
   
   % Minimize s_n(y)
   %  -  Solve  min[s_n(y)] with a local optimizer by starting from
   %     the interpolation point with the least function value

   % TOMLAB/SOL
      
   Iter        = Iter + 1;
   snProb.snP  = Iter;
   snProb.x_0  = x_min;
   snProb.CGO  = gnProb.CGO;
   %snProb.RBF  = RBF;

   %snProb = glcAssign('sn_f',x_LL,x_UU,'RBFsn');
   %snProb.optParam.fGoal=[];
   %snResult=glbSolve(snProb);
   %min_sn     = snResult.f_k
   %min_sn_y   = snResult.x_k;

   % Local optimum using e.g. NPSOL

%snProb.SOL.PrintFile = 'npsol.txt';
%snProb.SOL.SummFile  = 'npsols.txt';
%snProb.SOL.optPar(1) = 11111;
%snProb.SOL.optPar(13) = 3;

   snResult = tomRun(localSolver,snProb,max(PriLev-4,DebugPriLev*2));

%  Diffx = snResult.x_k(:,1)-x_min;
%  fprintf('Dist-2-1 %f %f \n',norm(Diffx,2),norm(Diffx,1));
%  snResult2 = tomRun('oqnlp',snProb,max(PriLev-4,DebugPriLev*2));
%  Diffx2 = snResult2.x_k(:,1)-x_min;
%  fprintf('Dist-2-1 %f %f \n',norm(Diffx2,2),norm(Diffx2,1));
%  fprintf('Diff %f %f \n ',norm(Diffx,2)-norm(Diffx2,2),...
%                           norm(Diffx,1)-norm(Diffx2,1));
%  if norm(snResult.x_k(:,1)-snResult2.x_k(:,1),1)...
%         /norm(snResult.x_k(:,1),1) > 1E-3
%     'Differerent solutions'
%     fprintf('%f ',Diffx)
%     fprintf('\n')
%     fprintf('%f ',snResult.x_k(:,1)-snResult2.x_k(:,1))
%     fprintf('\n')
%     keyboard
%  elseif abs(norm(Diffx,2)-norm(Diffx2,2))/norm(snResult.x_k(:,1),2) > 1E-3
%     'Differerent differences'
%     fprintf('%f ',Diffx)
%     fprintf('\n')
%     fprintf('%f ',snResult.x_k(:,1)-snResult2.x_k(:,1))
%     fprintf('\n')
%     keyboard
%  elseif abs(norm(Diffx,1)-norm(Diffx2,1)) > 1E-5
%     fprintf('%f ',Diffx)
%     fprintf('\n')
%     fprintf('%f ',snResult.x_k(:,1)-snResult2.x_k(:,1))
%     fprintf('\n')
%  end

   ExitFlag = snResult.ExitFlag;
   %if isempty(snResult.x_k)
   if ExitFlag > 0
      if strcmpi(localSolver,'glcCluster') | strcmpi(localSolver,'glcFast')
         %if GOMaxFunc >= 2500
         %   snProb.optParam.MaxFunc = GOMaxFunc;
         %   snProb.optParam.MaxIter = GOMaxIter;
         %else
         %   snProb.optParam.MaxFunc = 2*GOMaxFunc;
         %   snProb.optParam.MaxIter = 2*GOMaxIter;
         %end
         %if strcmpi(globalSolver,'glcCluster')
         %   snProb.GO.maxFunc1    = 2*maxFunc1;
         %end
         snProb.WarmStart = 1;
         if IterPrint | PriLev > 0
            fprintf('Local min_sn with global solver: No feasible point, ');
            fprintf('Warm Start %s',localSolver);
            fprintf(' with MaxFunc %d',snProb.optParam.MaxFunc);
            fprintf(' and MaxIter %d\n',snProb.optParam.MaxIter);
         end
         snResult    = tomRun(localSolver,snProb);
         PrintResult(snResult,double(PriLev > 0))
         %snProb.optParam.MaxFunc = GOMaxFunc;
         %snProb.optParam.MaxIter = GOMaxIter;
         %snProb.GO.maxFunc1      = maxFunc1;
         snProb.WarmStart = 0;
      elseif 0
         % Try the global solver
         mFunc   = snProb.optParam.MaxFunc;
         snProb.optParam.MaxFunc = GOMaxFunc;
         mIter   = snProb.optParam.MaxIter;
         snProb.optParam.MaxFunc = GOMaxIter;
         snResult    = tomRun(globalSolver,snProb);
         PrintResult(snResult,double(PriLev > 0))
         snProb.optParam.MaxFunc = mFunc;
         snProb.optParam.MaxIter = mIter;
      end
   end
   if ~isempty(IntVars)
      f_k  = snResult.f_k;
      x_k  = snResult.x_k(:,1);
      if any(x_k(IntVars) ~= round(x_k(IntVars)))
         x00 = max(x_L(IntVars),min(x_U(IntVars),round(x_k(IntVars))));
         snProb.x_0 = x_k;
         snProb.x_0(IntVars) = x00;
         snProb.x_L(IntVars) = x00;
         snProb.x_U(IntVars) = x00;
         snResult = tomRun(localSolver,snProb,max(PriLev-4,DebugPriLev*2));
         %snResult = tomRun(localSolver,snProb,3);
         ExitFlag = snResult.ExitFlag;
         snProb.x_L = x_LL;
         snProb.x_U = x_UU;
      end
   end
   if isempty(snResult.x_k)
      % Should not occur
      if snResult.ExitFlag == 0
         fprintf('rbfSolve - local problem: Failure in global solver!')
         fprintf(' Not feasible\n')
      else
         fprintf('rbfSolve - local problem: Failure in global solver!\n')
      end
      fprintf('ExitFlag = %d. ExitText: ', snResult.ExitFlag);
      fprintf('%s.', snResult.ExitText);   
      %pause
      min_sn      = Inf;
      min_sn_y    = x_min;
   else
      min_sn      = snResult.f_k;
      min_sn_y    = snResult.x_k(:,1);
   end
   %if PriLev > 4
   %   PrintResult(snResult,PriLev-4)
   %end
   %PrintResult(snResult,3)
   %snResult    = tomRun('glcCluster',snProb,max(PriLev-4,DebugPriLev*2));
   %PrintResult(snResult,3)
   z           = min_sn_y;

   %if 0
   %   % Check if global solution == local solution on the surface
   %   snProb.optParam.MaxFunc=GOMaxFunc;
   %   snProb.optParam.MaxIter=GOMaxIter;
   %   snResult    = tomRun(globalSolver,snProb);
   %   fGlob       = snResult.f_k;
   %   snProb.optParam.MaxIter  = MaxIter;
   %   snProb.optParam.MaxFunc  = MaxIter*10;  
   %end

   % Transform to original space   
   if SCALE
      min_snc = tomsol(9, x_L, min_sn_y, x_D); 
      c0      = tomsol(9, x_L, x_min, x_D); 
   else
      min_snc = min_sn_y; 
      c0      = x_min;
   end


   % Choose a target value fnStar or compute alpha (implicitly gives fnStar)

   modN = modN + 1;
   if modN > N
      if FLOWER < nFunc - N
         %modN = -1
         modN = 0-infStep;
      else
         modN = 0-infStep;
      end
   end
   ucOK = 1;
   if modN == -1
      alpha  =  inf;
      fnStar = -inf;
   elseif idea == 1 
      if modN == 0
         nmax = n;
      else
         if DeltaRule
            % Current idea used:
            nmax = min(n,max(2,nmax-floor((nFunc-nInit)/N)));
         else
            % Always use ALL points
            nmax = n;
         end
      end
      
      if OldStrat
         F_sort = sort(F_m);
         max_F = F_sort(nmax);
      else
         if REPLACE & nmax > n/2
            max_F = median(F);
         else
            F_sort = sort(F);
            max_F = F(nmax);
         end
      end
      
      fOld = fnStar;
      if f_Low > -1E299
         fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
      else
         fDiff = max_F-min_sn;
      end
      fDiff = max(fDiff,max_F-fMin);

      fDiff = AmpFac * fDiff;

      if fStarRule == 1
         fnStar = min_sn - ( (N-modN)/N )^2*fDiff;
         % Check that the above is exactly the following:
         %fnStar = min_sn - ( (mod(N-(nFunc-nInit),N+1))/N )^2*(max_F-min_sn);
      elseif fStarRule == 2
         fnStar = min_sn - ( (N-modN)/N )*fDiff;
      elseif fStarRule == 3
         if modN == 0
            fnStar = min_sn;
         elseif abs(min_sn) < 1E-12
            fnStar = min_sn - ( (N-modN)/N );
         else
            fnStar = min_sn - ( (N-modN)/N )*abs(min_sn);
         end
      else
         if f_Low > -1E299
            fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
         else
            fDiff = max_F-min_sn;
         end
         switch modN
         case 0
           fnStar = min_sn - fDiff;
         case 1
           fnStar = min_sn - 0.5*fDiff;
         case 2
           fnStar = min_sn - 1E-1*fDiff;
         case 3
           fnStar = min_sn - 1E-2*fDiff;
         case 4
           fnStar = min_sn - 1E-4*fDiff;
         case 5
           fnStar = min_sn;
         end
      end

%     if modN==0
%        if abs(max_F-min_sn)<30
%           fnStar=fnStar-100;
%        end
%     end
      
      % Have made this test looser in old variant
      if abs(fOld - fnStar) < 1E-4  & (fStarRule < 3)
         disp('ALARM. fnStar nearly equal to fnStar in last step');
         nmax = nmax - 1;
         % Create fnStar without REPLACE
         F_sort = sort(F);
         max_F = F_sort(nmax);
         fnStar = min_sn - ( (N-modN)/N )^2*(max_F-min_sn);
         % Check that the above is exactly the following:
         % fnStar = min_sn - ( (mod(N-(nFunc-nInit),N+1))/N )^2*(max_F-min_sn);
      end
       
      if modN == N
         % Unconstrained cycle step. But min_sn must be sufficiently lower
   
         if fMin == 0
            maxF = max(F);
            if min_sn >= -eps_sn*min(1,maxF) 
               if maxF == 0
                  fnStar = -0.01;
               else
                  fnStar = -0.01*min(1,maxF);
               end
               ucOK = 0; 
            end
         elseif min_sn >= fMin-eps_sn*abs(fMin)
            fnStar = min_sn - 1E-2*abs(fMin);
            ucOK = 0;
         end
      end
      % fnvec=[fnvec fnStar];
      gnProb.CGO.fnStar = fnStar;
   elseif idea == 2
      max_F = max(F_m);

      % Generalized the fDiff strategy
      if f_Low > -1E299
         fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
      else
         fDiff = max_F-min_sn;
      end
      fDiff = max(fDiff,max_F-fMin);
      if fDiff <= 0
         %HKH New
         fDiff = max(0,max(F)-fMin);
         if fDiff <= 0
            fDiff = max(1E-4,max(F)-min(F));
         end
      end
      % Try improving search if failure during a whole cycle
      %if fDiffOld == fDiff & FLOWER < nFunc - N
      %   fDiff = fDiff*0.1;
      %end

      alphaPrev = alpha;
      if OldAlpha
       switch modN
          case 0
            alpha = fDiff/2;
            if alpha == 0
               alpha = 1;
            end
          case 1
            alpha = fDiff/2;
            if alpha == 0
               alpha = 0.5;
            end
          case 2
            alpha = min(fDiff/2,1);
            if alpha == 0
               alpha = 0.3;
            end
          case 3
            if fMin == 0
               %if (fMin-min_sn) <= eps_sn, ucOK = 0; end
               if abs(fMin-min_sn) <= eps_sn, ucOK = 0; end
            %elseif (fMin-min_sn) <= eps_sn*abs(fMin)
            elseif abs(fMin-min_sn) <= eps_sn*abs(fMin)
               ucOK = 0;
            end
            if ucOK
               alpha = 0;
            else
               alpha = min(fDiff/2,0.5);
               if alpha == 0
                  alpha = 0.1;
               end
            end
       end
      else
       switch modN
          case 0
            if abs(Update) == 1
               alpha = fDiff/2;
            elseif alphaPrev >0
               alpha = alphaPrev/2;
            else
               alpha = fDiff/2;
            end
            if alpha <= 0
               alpha = 1;
               %alpha = alphaXXX;
               %alpha = 20;
            end
          case 1
            if fMinIter == Iter-1 
               %???alpha = fDiff/2;
               alpha = fDiff/3;
            elseif fDiffOld == fDiff
               alpha = fDiff/4;
            elseif abs(Update) == 1
               %alpha = fDiff/2;
               %alpha = fDiff/3;
               alpha = fDiff/4;
            elseif alphaPrev >0
               alpha = alphaPrev/2;
            else
               alpha =  fDiff/2;
            end
            if alpha <= 0
               alpha = 0.5;
               %alpha = 0.5*alphaXXX;
               %alpha = 10;
            end
          case 2
            if fMinIter == Iter-1 
               %alpha = min(fDiff/2,1);
               % testa
               alpha = alphaPrev/2;
            elseif fDiffOld == fDiff
               if fDiff/2 > 1
                  alpha = 1;
               else
                  alpha = fDiff/4;
               end
            elseif abs(Update) == 1
               %alpha = min(fDiff/2,1);
               alpha = min(fDiff/8,1);
            else
               alpha = alpha/2;
            end
            if alpha <= 0
               alpha = 0.3;
               %alpha = 0.3*alphaXXX;
               %alpha = 5;
            end
          case 3
            if fMin == 0
               %URKif (fMin-min_sn) <= eps_sn, ucOK = 0; end
               if abs(fMin-min_sn) <= eps_sn, ucOK = 0; end
            %URKelseif (fMin-min_sn) <= eps_sn*abs(fMin)
            elseif abs(fMin-min_sn) <= eps_sn*abs(fMin)
               ucOK = 0;
            end
            if ucOK
               alpha = 0;
            else
               if fMinIter == Iter-1 
                  %alpha = min(fDiff/2,1);
                  % testa
                  alpha = alphaPrev/2;
               elseif fDiffOld == fDiff
                   if fDiff/2 > 1
                      alpha = 0.1;
                   else
                      %alpha = min(fDiff/4,0.1);
                      alpha = min(fDiff/8,0.1);
                   end
               elseif abs(Update) == 1
                   %alpha = min(fDiff/4,0.1);
                   alpha = min(fDiff/8,0.1);
               else
                   alpha = alpha/2;
               end
            end
       end
      end
      if ~abs(Update) == 1 & alphaPrev == alpha
         if alpha <= 0
            alpha = 1;
         else
            alpha = 2*alpha;
         end
         if PriLev > -1
            fprintf('Alpha same as previous step, adjust to %10.3f\n',alpha);
         end
      end
      gnProb.CGO.alpha = alpha;
   end

   % No point in not accepting local solution if Pure Integer problem
   if length(IntVars) == d, ucOK = 1; end
   
   % ********** MINIMIZE gn(y) **********
   if (modN == N) & ucOK
      if PriLev > 2
         fprintf('Local search OK, not minimizing gn_f\n')
      end
      xGlob = [];
      fGlob = [];
      xLoc  = min_sn_y;
      fLoc  = min_sn;
   end
   if 0 & (modN == N) & ~ucOK
      if PriLev > 2
         disp('*** Unconstrained solution STILL TOO CLOSE ***');
         disp('*** Try optimizing sn globally ***');
      end
      xQ = 1.00;   % Size of the region around the local solution
      if xQ ~=1
         % Shrink the box
         if SCALE
            snProb.x_L  = max(0,z-xQ);
            snProb.x_U  = min(1,z+xQ);
         else
            snProb.x_L  = max(x_LL,z-xQ*x_D);
            snProb.x_U  = min(x_UU,z+xQ*x_D);
         end
         if ~isempty(IntVars)
            % Must be safe guarded for integers
            snProb.x_L(IntVars) = floor(snProb.x_L(IntVars));
            snProb.x_U(IntVars) = ceil(snProb.x_U(IntVars));
         end
      end
      %fprintf('Number of evaluations = %d\n',snProb.optParam.MaxFunc);
      if strcmpi(localSolver,globalSolver) & strcmpi(globalSolver,'glcCluster')
         % Increase maxFunc1 and total number of function evaluations
         v = snProb.GO.maxFunc1;
         snProb.GO.maxFunc1=v+v;
         snProb.optParam.MaxFunc=snProb.optParam.MaxFunc+v+v;
         %zz=snProb.PriLevOpt;
         %snProb.PriLevOpt = 10;
         %snProb.optParam.IterPrint = 1;
         snResult = tomRun(globalSolver,snProb,max(PriLev-4,DebugPriLev*2));
         %snProb.PriLevOpt=zz;
         %snProb.optParam.IterPrint = 0;
         snProb.GO.maxFunc1=v;
         snProb.optParam.MaxFunc=snProb.optParam.MaxFunc-v-v;
      else
         snResult = tomRun(globalSolver,snProb,max(PriLev-4,DebugPriLev*2));
      end

      ExitFlag = snResult.ExitFlag;
      %if isempty(snResult.x_k)
      if ExitFlag == 7
         snProb.WarmStart = 1;
         if IterPrint | PriLev > 0
            fprintf('Last cycle step: No feasible point, ');
            fprintf('Warm Start %s',globalSolver);
            fprintf(' with MaxFunc %d',snProb.optParam.MaxFunc);
            fprintf(' and MaxIter %d\n',snProb.optParam.MaxIter);
            if strcmpi(globalSolver,'glcCluster')
               fprintf('maxFunc1 %d ',snProb.GO.maxFunc1);
               fprintf('maxFunc2 %d ',snProb.GO.maxFunc2);
               fprintf('maxFunc3 %d ',snProb.GO.maxFunc3);
               fprintf('\n');
            end
         end
         snResult                = tomRun(globalSolver,snProb);
         PrintResult(snResult,double(PriLev > 0))
         snProb.WarmStart = 0;
      end
      if isempty(snResult.x_k)
         % Should never occur
         if snResult.ExitFlag == 0
            fprintf('rbfSolve - global problem: Failure in global solver!')
            fprintf(' Not feasible\n')
         else
            fprintf('rbfSolve - local problem: Failure in global solver!\n')
         end
         fprintf('ExitFlag = %d. ExitText: ', snResult.ExitFlag);
         fprintf('%s.', snResult.ExitText);   
         %pause
      end
      xGlob  = snResult.x_k;
      fGlob  = snResult.f_k;
      if LOCAL
         LocSteps = min(20,size(xGlob,2));
         if LocSteps > 1 & PriLev > 0 
            fprintf('Do %d local search steps ',LocSteps);
            if size(xGlob,2) > LocSteps
               fprintf('out of %d\n',size(xGlob,2));
            end
            fprintf('\n');
         end
         % Do local search from globally best points
         snProb.snP  = Iter;
         % Max 1000 iterations in local solver
         snProb.optParam.MaxFunc  = MaxIter*max(d,10); 
         snProb.optParam.MaxIter  = MaxIter; 

         %snProb.CGO  = snProb.CGO;

         fLoc = inf;
         iBest = 0;
         xBest = [];
         for i = 1:LocSteps
             snProb.x_0  = xGlob(:,i);
             if PriLev == 4
                xprint(snProb.x_0,'x0:   ');
             end
             if PriLev > 4
                fprintf('Try local search #%d\n',i);
             end
             snR = tomRun(localSolver,snProb,max(PriLev-4,DebugPriLev*2));
             OK  = 1;
             %if snR.ExitFlag ~= 0, OK = 0; end
             if dLin > 0
                % Check linear constraints
                Ax = snProb.A*snR.x_k;
                AxMax = max(abs(Ax));
                if any(snProb.b_L - Ax > min(1E-5,1E-5*AxMax)  ...
                     | Ax - snProb.b_U > min(1E-5,1E-5*AxMax))
                   %'DO NOT ACCEPT LOCAL POINT'
                   %'LINEAR CONSTRAINTS NOT FULFILLED'
                   %ExitFlag = snR.ExitFlag
                   %Inform   = snR.Inform
                   OK = 0;
                end
             end
             if snR.f_k < fLoc  & OK
                iBest    = i;
                snResult = snR;
                xBest    = snR.x_k;
                fLoc     = snR.f_k;
                if PriLev == 4
                   xprint(xBest,'xOpt: ');
                end
             end
         end
         xLoc        = xBest;

         if PriLev > 2
            fprintf('Global f(x) %30.20f       with %s\n',fGlob,globalSolver);
            fprintf('Local  f(x) %30.20f (#%2d) with %s\n',...
                     fLoc,iBest,localSolver);
         end
      else
         if isempty(xGlob)
            xLoc       = x_min;
            fLoc       = Inf;
         else
            xLoc       = xGlob(:,1);
            fLoc       = fGlob;
         end
         if PriLev > 3
            xprint(xLoc,'x: ');
         end
         if PriLev > 2
            fprintf('Global f(x) %30.20f with %s\n',fGlob,globalSolver);
         end
      end
      if PriLev > 3
         xprint(xLoc,'x: ');
      end
      if xQ ~=1
         % Set back the full box
         snProb.x_L  = x_LL;
         snProb.x_U  = x_UU;
      end
   end
   if modN ~= N
      
      gnProb.CGO.modN = modN; 

      %fprintf('Number of evaluations = %d\n',gnProb.optParam.MaxFunc);
      % gnProb.optParam.IterPrint = 1;
      gnResult    = tomRun(globalSolver,gnProb,max(PriLev-4,DebugPriLev*2));
      %if isempty(gnResult.x_k)
      ExitFlag = gnResult.ExitFlag;
      if ExitFlag == 7
         % Extra emergency step, no feasible solution found
         gnProb.WarmStart = 1;
         if IterPrint | PriLev > 0
            fprintf('Cycle %d, No feasible point, ',modN);
            fprintf('Warm Start %s',globalSolver);
            fprintf(' with MaxFunc %d',gnProb.optParam.MaxFunc);
            fprintf(' and MaxIter %d\n',gnProb.optParam.MaxIter);
            if strcmpi(globalSolver,'glcCluster')
               fprintf('maxFunc1 %d ',gnProb.GO.maxFunc1);
               fprintf('maxFunc2 %d ',gnProb.GO.maxFunc2);
               fprintf('maxFunc3 %d ',gnProb.GO.maxFunc3);
               fprintf('\n');
            end
         end
         gnResult                = tomRun(globalSolver,gnProb);
         PrintResult(gnResult,double(PriLev > 0))
         gnProb.WarmStart = 0;
      end
      if isempty(gnResult.x_k)
          % Should never occur
         if gnResult.ExitFlag == 0
            fprintf('rbfSolve - global problem: Failure in global solver!')
            fprintf(' Not feasible\n')
         else
            fprintf('rbfSolve - local problem: Failure in global solver!\n')
         end
         fprintf('ExitFlag = %d. ExitText: ', gnResult.ExitFlag);
         fprintf('%s.', gnResult.ExitText);   
      end
      xGlob  = gnResult.x_k;
      fGlob  = gnResult.f_k;
      if PriLev > 3
         if ~isempty(xGlob), xprint(xGlob(:,1),'x: '); end
      end
      if LOCAL
         LocSteps = min(20,size(xGlob,2));
         if LocSteps > 1 & PriLev > -2 
            fprintf('Do %d local search steps ',LocSteps);
            if size(xGlob,2) > LocSteps
               fprintf('out of %d\n',size(xGlob,2));
            end
            fprintf('\n');
         end
         % Do local search from globally best points
         gnProb.snP  = Iter;
         gnProb.optParam.MaxFunc  = MaxIter*max(d,10); 
         gnProb.optParam.MaxIter  = MaxIter; 

         %snProb.CGO  = gnProb.CGO;

         fLoc  = inf;
         iBest = 0;
         xBest = [];
         for i = 1:LocSteps
             gnProb.x_0  = xGlob(:,i);
             if PriLev == 4
                xprint(gnProb.x_0,'x0:   ');
             end
             if PriLev > 4
                fprintf('Try local search #%d\n',i);
             end
             gnR = tomRun(localSolver,gnProb,max(PriLev-4,DebugPriLev*2));
             OK  = 1;
             %if gnR.ExitFlag ~= 0 & dLin > 0
             if dLin > 0
                % Check linear constraints
                Ax = gnProb.A*gnR.x_k;
                AxMax = max(abs(Ax));
                if any(gnProb.b_L - Ax > max(1E-5,1E-5*AxMax)  ...
                     | Ax - gnProb.b_U > max(1E-5,1E-5*AxMax))
                   %'DO NOT ACCEPT LOCAL POINT'
                   %'LINEAR CONSTRAINTS NOT FULFILLED'
                   %ExitFlag = gnR.ExitFlag
                   %Inform   = gnR.Inform
                   OK = 0;
                end
             end
             if gnR.f_k < fLoc  & OK
                iBest    = i;
                gnResult = gnR;
                xBest    = gnR.x_k;
                fLoc     = gnR.f_k;
                if PriLev == 4
                   xprint(xBest,'xOpt: ');
                end
             end
         end
         gnProb.optParam.MaxIter  = GOMaxIter; 
         gnProb.optParam.MaxFunc  = GOMaxFunc; 
         if iBest == 0
            xLoc                     = xGlob(:,1);
         else
            xLoc                     = xBest;
         end
         %if LocSteps > 1 & iBest > 1 & PriLev > -2
         %   fprintf('Best local point found in search step %d\n',iBest);
         %end

         if PriLev > 2
            fprintf('Global f(x) %30.20f with %s\n',fGlob,globalSolver);
            fprintf('Local  f(x) %30.20f (#%2d) with %s\n',...
                     fLoc,iBest,localSolver);
         end
      else
         if isempty(xGlob)
            xLoc       = [];
            fLoc       = Inf;
         else
            xLoc       = xGlob(:,1);
            fLoc       = fGlob;
         end
         if PriLev > 3
            xprint(xLoc,'x: ');
         end
         if PriLev > 2
            fprintf('Global f(x) %30.20f with %s\n',fGlob,globalSolver);
         end
      end
      if PriLev > 3
         xprint(xLoc,'x: ');
      end
   end
   
   % Best point found on surface is xLoc - set as xNew
   xNew = xLoc;
   % New point in original space
   if SCALE
      O_new   = tomsol(9, x_L, xNew, x_D); 
   else
      O_new   = xNew;
   end
   
   % ********** PLOTTING **********
   if PLOT & d > 1 
      
      if 10 & ~(modN == N)
         switch lower(globalSolver)
         case 'glcfast'
              glb=load('glcFastSave.mat','C');
         case 'glbfast'
              glb=load('glbFastSave.mat','C');
         case 'glbsolve'
              glb=load('glbSave.mat','C');
         case 'glcsolve'
              glb=load('glcSave.mat','C');
         otherwise
              break;
         end
         % Transform to original space
         if SCALE
            CCC   = tomsol(9, x_L, glb.C, x_D); 
         else
            CCC   = glb.C;  
         end
         CCC = zeros(size(glb.C));
         
         plot(CCC(1,:),CCC(2,:),'.r');
         hold on
         plot(O(1,:),O(2,:),'*');
      else
         plot(O(1,:),O(2,:),'*');
         hold on
         TT = 0:0.1:1;
         YY(1,:) = c0(1)+TT*(O_new(1)-c0(1));
         YY(2,:) = c0(2)+TT*(O_new(2)-c0(2));
         plot(YY(1,:),YY(2,:),'-g')
      end   
      plot(O_new(1),O_new(2),'*g');
      x_opt = Prob.x_opt;
      if ~isempty(x_opt)
         if min(size(x_opt))==1
            x_opt = x_opt(:);
         end
         plot(x_opt(1,:),x_opt(2,:),'*k');
      end
      hold off
      grid on
      zoom on
   end % if plot
   if PLOT & d == 1 
         switch lower(globalSolver)
         case 'glcfast'
              glb=load('glcFastSave.mat','C');
         case 'glbfast'
              glb=load('glbFastSave.mat','C');
         case 'glbsolve'
              glb=load('glbSave.mat','C');
         case 'glcsolve'
              glb=load('glcSave.mat','C');
         otherwise
              break;
         end
         % Transform to original space
         if SCALE
            CCC   = tomsol(9, x_L, glb.C, x_D); 
         else
            CCC   = glb.C;  
         end
         CCC = zeros(size(glb.C));
      plot(CCC(1,:),zeros(size(CCC,2)),'.r');
      plot(O_new(1),fLoc,'*g');
      x_opt = Prob.x_opt;
      if ~isempty(x_opt)
         if min(size(x_opt))==1
            x_opt = x_opt(:);
         end
         plot(x_opt(1,:),x_opt(2,:),'*k');
      end
      hold off
      grid on
      zoom on
   end
   
   % *************** UPDATE ***************
   % Remove information for response surface problem
   global NLP_x NLP_f NARG
   NLP_x=[]; NLP_f=[]; NARG = [];

   if isempty(O_new) % Infeasibility problem
      fNew   = Inf;
      fPen   = Inf;
      SAME1  = 0;
      SAME2  = SAME2+1;
      Update = -5;
   else
      if all(O_new == O_min), SAME1 = SAME1 + 1; else SAME1 = 0; end
      if all(O_new == O_pre), SAME2 = SAME2 + 1; else SAME2 = 0; end
      ix = find(all(X==xNew*ones(1,size(X,2))));
      if isempty(ix)
         Update = 1;
         fNew  = nlp_f(O_new, Prob, varargin{:});
         nFunc = nFunc + 1; % Total number of sampled points
         fPen = fNew;
         if dLin > 0
            L = Prob.A*O_new;
            fPen = fPen+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
         end
         if dCon > 0 
            C = nlp_c(O_new, Prob, varargin{:});
            nCon = nCon + 1;
            fPen = fPen+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
         end
      else
         Update = -4;
         fNew   = NaN;
      end
   end
   time  = fix(clock);

   if REPLACE & Update == 1
      fnewV=min(median([F;fNew(1)]),[F;fNew(1)]);
      Update = tomsol(24, xNew, fnewV);    
   elseif  Update == 1
      Update = tomsol(24, xNew, fNew);    
   end
     
   %if Update < 0
   %   Update=0;
   %else
   %   Update=1;
   %end
   
   % if minDist < minDistEps
   %    if fNew >= F(ixC)
   %       Update =  0;
   %    else
   %       fprintf('interchange %d with new %d\n',ixC,n+1);
           % Update =  -1;
           % pause
   %    end
   % else
   %    Update = 1;
   % end

   if Update == 1 | Update == -1
      % All is OK, Update =-1, refactorization made the trick
      F   = [F;fNew];
      Fpen= [Fpen;fPen];
      X   = [X xNew];
      O   = [O O_new];
      O_pre = O_new;
      n   = n + 1;
      
      if REPLACE 
         F_m = min(median(F),F);
      else
         F_m = F; % No replacement
      end
      % Is the function value in the new point less than fMin?
      % If feasible, compare fNew,fMin, otherwise fPen and fMin
      % Or now feasible, but previously not?
      if (Feasible & fNew < fMin & fPen==fNew) | ...
         (~Feasible & fPen < fMin) | (~Feasible & fPen == fNew)
         Feasible = fNew == fPen;
         fMin     = fPen;
         fMinIter = Iter;
         fIdx     = length(Fpen);
         x_min    = xNew;
         O_min    = O_new;
         FLOWER   = nFunc;
         FLOW     = nFunc;
         alphaXXX = alpha;
      end
      NOUPDATE = 0;
   elseif Update == -2
      % Update == -2 New point bad, even refactorization failed
      % Infeasible problem, ill-conditioning?

      Dist = tomsol(30,X(:,fIdx),X);
      Dist(fIdx) = Inf;

      control = -1;
      VALUE = 1;
      while control < 0
         tomsol(25) % Deallocates memory
         fprintf('Minimal distance to X set');
         [minDist ixDist] = min(Dist);
         fprintf(' %s',minDist);
         fprintf('\n');
         ix = ones(n,1);
         ix(ixDist)=VALUE;
         % Remove most infeasible point
         [maxPen ixF] = max(Fpen-F);
         if maxPen == 0
            % If all points feasible, remove point with largest f(x)
            [maxPen ixF] = max(Fpen);
         end
         ix(ixF)=VALUE;
         ix   = find(ix);
         F    = F(ix);
         Fpen = Fpen(ix);
         X    = X(:,ix);
         O    = O(:,ix);
         n    = size(X,2);
         Dist = Dist(ix);
         fIdx = find(isinf(Dist));
         if VALUE == 1
            REPLACE=1-REPLACE
         end
         if REPLACE 
            F_m = min(median(F),F);
         else
            F_m = F; % No replacement
         end
         % TOMSOL INIT, send F_m to Fortran, not F
         'make init again'
         control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
         control
         keyboard
         VALUE = 0;
      end
      if control < 0
          fprintf('New initial interpolation failed');
          tomsol(25) % Deallocates memory
          return     %Something is really wrong
      end
      NOUPDATE = 0;
   elseif Update == -3 | Update == -4 | Update == -5
      % Update == -3 Point too close to old point
      % Update == -4 Point identical to old point
      % Update == -5 No feasible point found on surface
      NOUPDATE = NOUPDATE+1;
      if PriLev > 1
         fprintf('!!!!! Update impossible\n')
      end
   end
   %if FLOWER < nFunc - N
   %   % Trying to improve a nonworking strategy, by changing replacement
   %   REPLACE = 1 - REPLACE;
   %end
   
   if PriLev > 1 | IterPrint
      fprintf('Iter %3d n %3d ', Iter, n);
      fprintf('nFunc %3d ', nFunc);
      tt = time([3 2 4 5]);
      if tt(4) < 10
         fprintf('%d/%d %d:0%1d  ', tt);
      else
         fprintf('%d/%d %d:%2d ', tt);
      end
      fprintf('Cycle %2d ', modN);
      if idea == 1
         fprintf(' fnStar%7.3f',fnStar);
      else
         fprintf(' alpha %7.3f',alpha);
      end
      if ~isempty(fGoal) & ~isinf(fGoal)
         fprintf(' fGoal %8.5f', fGoal);
      end
      if Feasible
         fprintf(' fMinF %11.8f ', fMin);
      else
         fprintf(' fMinI %11.8f ', fMin);
      end
      fprintf('at %3d/It %d fNew %11.8f ', FLOWER, fMinIter,fNew);
      fprintf(' fLoc %11.8f', min_sn);
      %NEWHKHfprintf(' %11.8f', fGlob);
      fprintf('\n');
      xprint(O_new,'xNew:',' %12.8f',8)
   end
   
   % ********** CONVERGENCE TEST **********
   if convflag == 0
      convflag = isClose(fGoal,fMin,fTol,nFunc,Iter,PriLev);
   end
   if NOUPDATE > N
      convflag = 4;
   elseif SAME1 > N
      convflag = 5;
   elseif SAME2 > N
      convflag = 6;
   end
   
   if PLOT
      pause(1)
   end
   if abs(Update) == 1
      fDiffOld = fDiff;
   end

   % -------------- Result saving -------------------
   %saveIter;
end

% *******************************************************
% *************** END MAIN ITERATION LOOP ***************
% *******************************************************

% SAVE RESULTS

fMinIdx = fIdx(1);
save('cgoSave.mat','Name','O','F','X','F_m','nInit','Fpen','fMinIdx');

% All points i with F(i)=f_min
if SCALE
   Result.x_k      = tomsol(9,x_L,X(:,find(F==fMin)),x_D);    
else
   Result.x_k      = X(:,find(F==fMin));
end
if isempty(Result.x_k)
   if SCALE
      Result.x_k      = tomsol(9,x_L,X(:,find(Fpen==fMin)),x_D);    
   else
      Result.x_k      = X(:,find(Fpen==fMin));
   end
   fprintf('Warning: Optimal point is found to be infeasible: ');
   fprintf('constraints are violated\n')
end

Result.f_k      = fMin;     % Best function value
cMin = [];
if dCon > 0
   for i=1:size(Result.x_k,2)
       cMin = [cMin,nlp_c(Result.x_k(:,i), Prob, varargin{:})];
       nCon = nCon + 1;
   end
end
if dLin > 0 & size(Result.x_k,2)==1
   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at best x_k
end
Result.c_k      = cMin;     % Constraint value at best x_k
Result.Iter     = Iter;     % Number of iterations
Result.FuncEv   = nFunc;
Result.ConstrEv = nCon;
Result.SolverAlgorithm = [Result.SolverAlgorithm ...
                          '. Global solver ' globalSolver  ...
                          '. Local solver ' localSolver];
Result.ExitFlag = 0;
Result.ExitText = ['Tried ' num2str(nFunc) ' f(x), using ' ...
                    num2str(n) ', startup ' num2str(nInit)];
if cpumax & convflag == 0
   Result.Inform   = 9;
   Result.ExitText = [Result.ExitText '. Max CPU reached. '];
elseif ~progress & convflag == 0
   Result.Inform   = 8;
   Result.ExitText = [Result.ExitText '. No progress for ' ...
                      num2str(nFunc-FLOWER) ' function evaluations'];
else
   Result.Inform   = convflag;
end
%TIME0=TIME0Save;
%TIME1=TIME1Save;
Result          = endSolve(Prob,Result);

% -------------- Result saving -------------------
%saveResult;

tomsol(25) % Deallocates memory

function convflag = isClose(fGoal,f,fTol,nFunc,Iter,PriLev)

convflag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convflag = 1;
elseif fGoal == 0
   %if abs(f-fGoal) < fTol
   if abs(f) < fTol
      convflag = 2;
   end
elseif abs(f-fGoal) <= abs(fGoal) * fTol
   convflag = 3;
end

if convflag > 0 & PriLev > 0 
   if convflag == 1
      fprintf('\n\nFunction value %f is less than fGoal %f \n',f,fGoal);
   elseif convflag == 2
      fprintf('\n\nError in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal),fTol);
   elseif convflag == 3
      fprintf('\n\nRelative error in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal)/abs(fGoal),fTol);
   end
   fprintf('Number of function evaluations:  %d\n',nFunc);
   fprintf('Number of iterations:            %d\n',Iter);
end

% ====================================================================
function X = corners(x_L,x_U)
% ====================================================================

d  = length(x_L);

ix = x_L~=x_U;
iV = find(ix);
iF = find(~ix);
m  = length(iV);

n = 2^m;

X = zeros(d,n);
for i = 1:length(iF)
    j = iF(i);
    X(j,:)=x_L(j);
end

for j=1:m
    var = iV(j);
    for i=1:2^j
        l=n/2^j;
        if mod(i,2)==1
           X(var,l*(i-1)+1:l*i)=x_L(var);
        else
           X(var,l*(i-1)+1:l*i)=x_U(var);
        end
     end
end
% ====================================================================
function X = randomtest(x_L,x_U,proc)
% ====================================================================

d = length(x_L);

n = 2*d;

dist=abs(x_U-x_L);

minDist=proc/100*min(dist);

%X=rand(d,1).*dist+x_L;

X=x_L+dist/2; %center point
x_n=rand(d,1).*dist+x_L;

for i=1:n-1
    [t,m]=size(X);
    while min(sqrt(sum((x_n*ones(1,m)-X).^2)))<minDist
        x_n=rand(d,1).*dist+x_L;
    end
    X=[X x_n];
end

% ====================================================================
function X = gutmann(x_L,x_U)
% ====================================================================

d  = length(x_L);
iV = find(x_L~=x_U);
m  = length(iV);

n = length(iV)+1;

X = x_L*ones(1,n);

for i=1:m
    vars = iV(i);
    X(i,i+1) = x_U(i);
end

% ====================================================================
function X = sampleInts(X,x_L,x_U,IntVars);
% ====================================================================
d = length(x_L);
n = size(X,2);
for i = 1:length(IntVars)
    j = IntVars(i);
    X(j,:) = x_L(j)+floor(rand(1,n)*(x_U(j)-x_L(j)+1));
end

% NOT USED:
% FACTORIZE    0-No factorization (default), 1-Factorisation.
%              NOTE! Not used in this version, only in old rbfSolve.
% rbflib uses tiny = 1E-8, as tolerance for the distance being too small
%?? eps_x     Convergence tolerance in x. All possible rectangles are 
%             less than this tolerance (scaled to (0,1) )

% MODIFICATION LOG
%
% 000112  mbk  First version written
% 001101  hkh  Revision for new tests, speedups
% 010431  hkh  Total revision, making it similar to glbSolve
% 010723  hkh  Send MaxFunc to Mex
% 011110  hkh  Fixes for GUI. Revision.
% 011206  hkh  Alternative cycle strategies are tried.
% 020103  hkh  Use general field CGO and cgoSave.mat for warm start
% 020614  hkh  Correcting errors in print statements
% 020820  hkh  Print sub problem results dependent on PriLev
% 020820  hkh  Try local search for all possible global solutions
% 020823  hkh  Error in transforming linear constraints if SCALE = 1
% 020824  hkh  Accept local solution only if linear constraints fulfilled
% 020828  hkh  Use GetSolver to select default local solver
% 020831  hkh  Add parameters in Prob.GO as input parameters
% 020831  hkh  Add CGO.LOCAL, cycle length CGO.N as input parameters
% 020831  hkh  Add CGO.infStep, add inf target value step in cycle
% 020831  hkh  Add CGO.AddMP, add midpoint to corner strategy if true
% 020831  hkh  Add CGO.fStarRule and CGO.DeltaRule
% 020831  hkh  Added gutmann d+1 initial value strategy
% 020901  hkh  Create two more fStar search strategies, now all Gutmann
% 020904  hkh  Set GOMaxFunc and GOMaxIter empty if not initialized.
% 021020  hkh  MaxFunc and MaxIter wrongly set empty if Prob.GO empty
% 021020  hkh  MaxFunc and MaxIter wrongly set empty if Prob.GO empty
% 021020  hkh  Use Fpen=F+sum(max(0,dev)) to handle infeasible initial points
% 030221  hkh  isempty(GOLocalSolver) should be isempty(GO.LocalSolver)
% 030222  hkh  Add check for empty Result.x_k, use best infeasible solution
% 030223  hkh  randomtest used size as variable name, changed to m
% 030223  hkh  Changed variable names, (xNew,fNew) for f(x) point
%              (xGlob,fGlob) for all global opt solutions on surface, 
%              (xLoc,fLoc) for best solution found on surface
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define fields mLin and mNonLin in snProb and gnProb
% 040226  hkh  Set field Prob.MIP in gnProb.MIP and snProb.MIP
% 040226  hkh  Set globalSolver='glcCluster' if Prob.MIP.IntVars is non-empty
% 040303  hkh  Making rbfSolve handle mixed-integer problems
% 040306  hkh  Add CGO.RandState, to initialize the random generator
% 040306  hkh  Check for cycle of same points
% 040306  hkh  Set Result.Inform to convflag
% 040307  hkh  Revised functions gutmann and corners, added new sampleInts
% 040307  hkh  Revise feasiblity handling, print fMinI, fMinF, use Feasible
% 040308  hkh  Special treatment of pure IP, with glcSolve and glcFast
% 040318  hkh  Major tuning of algorithm, Update codes etc.
% 040318  hkh  Set CGO.F as column always. Revised handling of initial F,X
% 040319  hkh  New optional input CGO.CX for initial constraint values
% 040321  hkh  Set initial FLOWER to fIdx(1), not nFunc
% 040323  hkh  Set cTol and bTol in suboptimizations, avoiding infeasibility
% 040323  hkh  Default for Percent if dLin+dCon > 0, use -5000 instead of -d 
% 040404  hkh  Major revision, global solvers do not return empty
% 040406  hkh  Save fMinIdx (=fIdx(1)) on cgoSave.mat
% 040412  hkh  Add calls to ProbCheck
% 040414  hkh  Do not set MENU or LargeScale in gnProb and snProb
% 040414  hkh  Return Result.c_k and Result.Ax
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 040928  hkh  sampleInts, x_L(j), offset from 0, must be added
% 040928  hkh  If localSolver is glcCluster, GOMaxFunc must be used
% 041005  hkh  Use abs(fMin-min_sn) <= eps_sn*abs(fMin), not fMin-min_sn
% 041005  hkh  Revise STILL TOO CLOSE handling if glcCluster
% 041005  hkh  Adding ExitFlag=8 for no progress for MaxCycle*(N+1) func evals
% 041005  hkh  CGO.MaxCycle new input, default 10
% 041005  hkh  Add localSolver and globalSolver to SolverAlgorithm output
% 041006  hkh  If Unc Sol too close, use global solver on snProb, not gnProb
% 041018  hkh  Set snProb.GO = Prob.GO, because global solver sometimes used
% 041019  hkh  Avoid global step if "TOO CLOSE", set eps_sn = 1E-6
% 041114  hkh  Add parameters for solver minlpSolve
% 041114  hkh  Rerun with rounded IntVars if local solution not int feasible
% 041114  hkh  Change default strategy to idea 1
% 041123  hkh  Change call to tomRun
% 041130  hkh  Remove unused code (previously used when ~ucOK and modN==N)
