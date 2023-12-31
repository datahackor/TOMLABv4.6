% ego.m
% 
% ego implements the algorithm EGO by D. R. Jones,
% Matthias Schonlau and William J. Welch presented in the paper
% "Efficient Global Optimization of Expensive Black-Box Functions".
% All page references and notations are taken from this paper.
%
% ego solves problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finite
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% f(x) are assumed to be a costly function
% c(x) are assumed to be cheaply computed
% Some or all x may be integer valued as specified by other input variables
%
% If some subset c_I(x) are very costly, create f(x) as a penalty function as
%
%     NEW f(x) = f(x) + beta' * c_I(x), beta positive penalties
%
% Calling syntax:
%
% function Result = ego(Prob, varargin)
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   USER.f    The routine to compute the function, given as a string, say EGOF
%   USER.c    The routine to compute the nonlinear constraint, say EGOC
%             A call to tomFiles.m or glcAssign.m sets these fields. 
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   WarmStart If true, >0, ego reads the output from the last run
%             from the mat-file cgoSave.mat, and continues from the last run.
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%   PriLevOpt Print Level 
%             0 = silent. 1 = Summary 2 = Printing each iteration
%             3 = Info about local / global solution 4 = Progress in x
%   user      User field used to send information to low-level functions
% ---------------------------------------------------------
% optParam    Structure in Prob, Prob.optParam 
% ---------------------------------------------------------
%             Defines optimization parameters. Fields used: 
%  IterPrint  Print one information line each iteration, and the new x tried
%             Default IterPrint = 1.
%  MaxIter    Maximal number of iterations used in the global optimization on
%             the response surface in each step. Default 10000.
%  MaxFunc    Maximal number of function evaluations in ego, default 200
%             If WarmStart == 1 and MaxFunc <= nFunc (Number of f(x) used)
%             then MaxFunc = MaxFunc + nFunc
%  fGoal      Goal for function value, if empty not used
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
% SCALE        0-Original search space
%              1-Transform search space to unit cube (default).
% PLOT         0-No plotting (default), 1-Plotting sampled points.
% REPLACE      0-No replacement
%              1-Large function values are replaced by the median, (default).
% LOCAL        0-No local searches after global search
%              If EGO surface is inaccurate, might be an advantage
%           
%              1-Local search from best points after global search. If equal
%              best function values, up to 20 local searches are done.
%
% TolExpI      Convergence tolerance for expected improvement, default 1E-6
%
% pEst         1-Estimate d-vector of || ||_p parameters (default), 0-fix p=2
%
% TRANSFORM    Function value transformation
%              TRANSFORM = 0  % No transformation made.
%              TRANSFORM = 1  % log(y) transformation made.
%              TRANSFORM = 2  % -log(-y) transformation made.
%              TRANSFORM = 3  % -1/y transformation made.
%              Default EGO is computing the best possible transformation
%              from the initial set of data.
%              Note! No check is made on illegal y if user gives TRANSFORM
%
% globalSolver Solver used for global optimization on the response surface
%              If the globalSolver is glcCluster, the fields
%              Prob.GO.maxFunc1, Prob.GO.maxFunc2 and Prob.GO.maxFunc3 are used
%              See the help for maxFunc1, maxFunc2, maxFunc3 in glcCluster
% localSolver  Solver used for local optimization on the response surface
%
% Percent      Type of strategy to get the initial sampled values:
%
%              Percent >= 100: User given initial points x, as a matrix of
%              points in CGO.X. Each column is one sampled point. 
%              If d = length(Prob.x_L), then size(X,1) = d, size(X,2) >= d+1
%              CGO.F should be defined as empty, or contain a vector of 
%              corresponding f(x) values. Any CGO.F value set as NaN will be 
%              computed by ego.
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
%              The algorithm will try up to M = abs(Percent-1000) points,
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
% AddMP        If = 1, add the midpoint as extra point in the corner strategy
%              or the Gutmann strategy. Default 0.
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
% will use the OQNLP solver, a license for TOMLAB /OQNLP is needed
%
% For pure IP problems, only glcSolve and glcFast (default) may be used. Set
% Prob.CGO.globalSolver = 'glcSolve'; to use glcSolve, otherwise glcFast is used
%
% varargin     Additional parameters to ego are sent to the costly f(x)
%
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
%           9 = Max CPU Time reached
%
% To make a warm start possible, ego saves the following information in
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
% USAGE:
%
% Let the name of the problem be "EGOF Test"
% The function EGOF is best written as
%     function f = EGOF(x, Prob) 
% Then any information, say u and W is easily sent to EGOF (and the constraint
% function EGOC, if present) using the Prob structure. See the example below. 
%
% Assume bounds x_L and x_U are initialized, as well as the linear constraint
% matrix A, lower and upper bounds, b_L and b_U on the linear constraints,
% and lower and upper bounds c_L and c_U on the nonlinear constraints
% (Put [] if all bounds are inf or -inf). Use the TOMLAB Quick format:
%
%    Prob   = glcAssign('EGOF',x_L,x_U,'EGOF Test',A,b_L,b_U,'EGOC',c_L,c_U);
%    Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%    % Default values are now set for PriLevOpt, and structure optParam
%    % To change a value, examples are shown on the two next lines
%    Prob.optParam.MaxFunc = 300; % Change max number of function evaluations 
%    Prob.optParam.MaxIter = 200; % Change the number of iterations to 200
%
% Direct solver call:
%    Result = ego(Prob);
%    PrintResult(Result);
%           
% Driver call, including printing with level 2:
%      Result = tomRun('ego',Prob,2);
%           
% The user function EGOF is written as
%
%    function f = EGOF(x, Prob) 
%    u = Prob.user.u; W = Prob.user.W;
%    f = "some function of x, u and W"
%
% It is also possible to use the function format
%    function f = EGOF(x) 
% but then any additional parameters must be sent as global variables.
%
% The user function EGOC, computing the nonlinear constraints, is written as
%
%    function c = EGOC(x, Prob) 
%    u = Prob.user.u; W = Prob.user.W;
%    c = "some vector function of x, V and W"
%
% Note! If EGOF has the 2nd input argument Prob, also EGOC must have that.
%
% To make a restart, just set the restart flag, and call ego once again:
%
%    Prob.WarmStart = 1;
%    Result = tomRun('ego',Prob,2);
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 2.5.0$
% Written Sep 27, 1998.   Last modified Nov 23, 2004.
%

% idea         Type of search strategy on the surface, idea 1 (cycle of 6
%              points in target value fnStar) or idea 2 (cycle of 4 points in 
%              alpha, which implicitely gives the target value fnStar).
% rbfType      Type of radial basis function.
%              1-Thin Plate Spline, 2-Cubic Spline (default).

function Result = ego(Prob, varargin)

if nargin < 1
   error('ego needs one parameter structure Prob');
end

time  = fix(clock);

DebugPriLev = 0;  % PriLev in subopt, i.e. gnProb.PriLevOpt, snProb.PriLevOpt

solvType=checkType('glc');

Prob=ProbCheck(Prob,'ego',solvType);

Prob = iniSolve(Prob,solvType,0,0);

MaxCPU = Prob.MaxCPU;

global GOMaxFunc GOMaxIter globalSolver maxFunc1 maxFunc2 maxFunc3
global NARG 
NARGSAVE = NARG;

if isempty(Prob.x_L) | isempty(Prob.x_U)
   disp('ego requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'ego requires both lower and upper variable bounds';
   Result=endSolve(Prob,Result);
   return;
end

% Pick up input variables from Prob structure
PriLev    = Prob.PriLevOpt;          % Print level
f_Low     = Prob.f_Low;              % Lower bound on f
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds
MaxIter   = Prob.optParam.MaxIter;   % Iterations used in global subopt
MaxFunc   = Prob.optParam.MaxFunc;   % Number of function evaluations
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration

fTol      = Prob.optParam.eps_f;     % Relative convergence tolerance in f(x)
fGoal     = Prob.optParam.fGoal;     % Goal f(x) for the optimization
epsRank   = Prob.optParam.eps_Rank;  % Rank tolerance in qr-decomposition
bTol      = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol      = Prob.optParam.cTol;      % Constraint feasibility tolerance
% HARD CODED LOW RIGHT NOW
epsRank   = 1E-14;

%xTol      = Prob.optParam.eps_x;     % Tolerance for rectangle sizes
%                                     % (scaled to (0,1) )

%EpsF      = Prob.optParam.eps_f;     % Relative convergence tolerance in f(x)
%fOpt      = Prob.f_opt;
%if isempty(fOpt), fOpt = -Inf; end
%if fOpt == 0
%   StopTol = 0.01*EpsF;
%else
%   StopTol = 0.01*fOpt;
%end

x_D   = x_U - x_L;
d     = length(x_L);  % Problem dimension
nCon  = 0;            % Count number of calls to constraint routine, nlp_c

% Safeguard
if isempty(MaxIter), MaxIter = 1000; end
if MaxIter < 0
   MaxIter = 1000;
end
if isempty(MaxFunc), MaxFunc = 200; end
if MaxFunc < 0
   MaxFunc = 200;
end
if isempty(IterPrint), IterPrint = 1; end

if isempty(Prob.CGO)
   idea      = []; rbfType   = []; SCALE        = []; PLOT        = [];
   REPLACE   = []; Percent   = []; globalSolver = []; localSolver = [];
   LOCAL     = []; TolExpI   = []; pEst         = []; TRANSFORM   = []; 
   AddMP     = []; RandState = [];

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
   if isfield(Prob.CGO,'AddMP')
      AddMP = Prob.CGO.AddMP;
   else
      AddMP = [];
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
   if isfield(Prob.CGO,'TolExpI')
      TolExpI = Prob.CGO.TolExpI;
   else
      TolExpI = [];
   end
   if isfield(Prob.CGO,'pEst')
      pEst = Prob.CGO.pEst;
   else
      pEst = [];
   end
   if isfield(Prob.CGO,'TRANSFORM')
      TRANSFORM = Prob.CGO.TRANSFORM;
   else
      TRANSFORM = [];
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

Result                 = ResultDef(Prob);
Result.Solver          = 'ego';
Result.SolverAlgorithm = 'Efficient Global Optimization';

if any(isinf(x_L)) | any(isinf(x_U))
   disp('ego solves box-bounded problems.');
   disp('Found some bound to be Inf');
   Result.ExitFlag = 2;
   Result.ExitText =str2mat('ego solves box-bounded problems' ...
                           ,'Found some bound to be Inf');
   Result=endSolve(Prob,Result);
   return
end

DEBUG  = 0;     % Debugging flag
DEBUG  = 0;
DEBUG2 = 0;     % Debugging flag for extreme debugging.

if DEBUG | DEBUG2
   PriLev = 1;
end

x_D   = x_U - x_L;
d     = length(x_L);  % Problem dimension

CX    = [];

% Default strategies for CGO input
if isempty(RandState),    RandState = 0; end
if isempty(SCALE),        SCALE = 1; end
if isempty(PLOT),         PLOT = 0; end
if isempty(REPLACE),      REPLACE = 0; end
if isempty(LOCAL),        LOCAL = 1; end
if isempty(TolExpI),      TolExpI = 1E-6; end
if isempty(pEst),         pEst = 0; end
if isempty(AddMP),        AddMP = 0; end
if isempty(globalSolver), globalSolver = 'glcFast'; end
if isempty(localSolver),  localSolver = GetSolver('con',1,0); end
%NOT, keep empty,if isempty(TRANSFORM),    TRANSFORM = -1; end

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
   Prob.c_L = [Prob.c_L;-inf*ones(dCon-length(Prob.c_L),1)];
   Prob.c_U = [Prob.c_U; inf*ones(dCon-length(Prob.c_U),1)];
end
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
   end
else
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
   % control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
   % if control < 0
   %    fprintf('Initial interpolation failed');
   %    tomsol(25) % Deallocates memory
   %    return     %Something is really wrong
   % end
   
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
         error('Illegal input to EGO')
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
% -----------------------------------------------------------------------
%   STEP 1    Fit a DACE model to a set of initial points specified by
%             a "space-filling" experimental design.
% -----------------------------------------------------------------------
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
           end
      end
      X = X(:,find(ixOK));
      
      if size(X,2) < PNTS
         for i = 1:n
             if ixOK(i) == 0
                X = [X,O(:,i)];
             end
             if size(X,2) >= PNTS, break; end
         end
      end
      if PriLev > 0 | IterPrint
         iV = find(ixOK);
         fprintf('Found %d out of %d ',length(iV),PNTS);
         fprintf('feasible points at trials ');
         fprintf('%d ',iV);
         fprintf('\n');
      end
      % X(:,1:length(iV))
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

   nFunc = n; % Number of function evaluations.
end


dc  = Prob.USER.dc;
d2c = Prob.USER.d2c;

if SCALE > 0 & (dCon > 0 | dLin > 0)
   if dLin > 0
      b_L = Prob.b_L - Prob.A*x_L;
      b_U = Prob.b_U - Prob.A*x_L;
      if pEst
         A  = [ Prob.A*diag(x_D), zeros(dLin,d)];
         XL = [Prob.x_L;zeros(d,1)];
         XU = [Prob.x_U;3*ones(d,1)];
      else
         A  = Prob.A*diag(x_D);
         XL = Prob.x_L;
         XU = Prob.x_U;
      end
   else
      A   = [];
      b_L = [];
      b_U = [];
      XL = Prob.x_L;
      XU = Prob.x_U;
   end
   if dCon > 0
      if pEst
         DACEProb = conAssign('dace_f',[],[],[],XL,XU,'DACE',...
                               [Prob.x_0,2*ones(d,1)], [],[], ...
                               A,b_L,b_U, 'rbf_c','rbf_dc','rbf_d2c',...
                               Prob.ConsPattern, Prob.c_L, Prob.c_U);
      else
         DACEProb = conAssign('dace_f',[],[],[],XL,XU,'DACE', Prob.x_0,...
                               [],[],A,b_L,b_U, 'rbf_c','rbf_dc','rbf_d2c',...
                               Prob.ConsPattern, Prob.c_L, Prob.c_U);
      end

      DACEProb.xD          = x_D;
      DACEProb.xL          = x_L;
      DACEProb.cNargin     = xnargin(Prob.USER.c);
      if isempty(dc)
         DACEProb.dcNargin = 0;
      else
         DACEProb.dcNargin = xnargin(dc);
      end
      if isempty(d2c)
         DACEProb.d2cNargin = 0;
      else
         DACEProb.d2cNargin = xnargin(d2c);
      end
      DACEProb.c           = Prob.USER.c;
      DACEProb.dc          = dc;
      DACEProb.d2c         = d2c;
      DACEProb.SCALE       = 1;
   else
      if pEst
         DACEProb = conAssign('dace_f',[],[],[],XL,XU,'DACE',...
                               [Prob.x_0,2*ones(d,1)], [],[], A,b_L,b_U);
      else
         DACEProb = conAssign('dace_f',[],[],[],XL,XU,'DACE', Prob.x_0,...
                               [],[],A,b_L,b_U);
      end
   end

elseif dCon == 0 & dLin == 0
   if pEst
      DACEProb = conAssign('dace_f',[],[],[],[],[],'DACE',...
                           [Prob.x_0,2*ones(d,1)]);
   else
      DACEProb = conAssign('dace_f',[],[],[],[],[],'DACE',Prob.x_0);
   end
else
   if pEst
      if dLin > 0
         A = [Prob.A, zeros(dLin,d)];
      else
         A = Prob.A;
      end
      XL = [Prob.x_L;zeros(d,1)];
      XU = [Prob.x_U;3*ones(d,1)];
      DACEProb = conAssign('dace_f',[],[],[],XL,XU,'DACE',...
                           [Prob.x_0,2*ones(d,1)], ...
                           [], [], A, Prob.b_L, Prob.b_U, ...
                          'rbf_c','rbf_dc','rbf_d2c',...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
                           %Prob.USER.c, Prob.USER.dc, Prob.USER.d2c, ...
      DACEProb.cNargin     = xnargin(Prob.USER.c);
      if isempty(dc)
         DACEProb.dcNargin = 0;
      else
         DACEProb.dcNargin = xnargin(dc);
      end
      if isempty(d2c)
         DACEProb.d2cNargin = 0;
      else
         DACEProb.d2cNargin = xnargin(d2c);
      end
      DACEProb.SCALE       = 0;
      DACEProb.c           = Prob.USER.c;
      DACEProb.dc          = dc;
      DACEProb.d2c         = d2c;
   else
      DACEProb = conAssign('dace_f',[],[],[],[],[],'DACE',Prob.x_0, ...
                           [], [], Prob.A, Prob.b_L, Prob.b_U, ...
                           Prob.USER.c, Prob.USER.dc, Prob.USER.d2c, ...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
   end
end

DACEProb.dDim    = d;

optParam                  = optParamDef(localSolver,solvType,d,dCon,dCon+dLin);

DACEProb.optParam           = optParam;
DACEProb.optParam.IterPrint = DebugPriLev > 0;
DACEProb.optParam.MaxIter   = MaxIter;
DACEProb.optParam.MaxFunc   = MaxIter*max(d,10)/10;  
DACEProb.optParam.bTol      = bTol;
DACEProb.optParam.cTol      = cTol;
%DACEProb.optParam.wait = 0;
DACEProb.GradTolg           = Prob.GradTolg;
DACEProb.GradTolH           = Prob.GradTolH;
DACEProb.GradTolJ           = Prob.GradTolJ;

DACEProb.SOL                = Prob.SOL;
if isfield(Prob,'OQNLP')
   DACEProb.OQNLP              = Prob.OQNLP;
end
if isfield(Prob,'KNITRO')
   DACEProb.KNITRO             = Prob.KNITRO;
end

%DACEProb.CHECK              = 1; % Avoid test of structure
%DACEProb.LargeScale        = 1; % Avoid saving of steps
DACEProb.PriLevOpt          = DebugPriLev;

% Send user info
if isfield(Prob,'user')
   DACEProb.user         = Prob.user;
end
DACEProb.MIP             = Prob.MIP;
DACEProb.uP              = Prob.uP;
DACEProb.P               = Prob.P;
DACEProb.mLin    = dLin;
DACEProb.mNonLin = dCon;


theta=[]; my=[]; p=[]; invR=[];

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
      EGOProb = glcAssign('ego_f',x_LL,x_UU,'EGO', A, b_L, b_U,...
                          'ego_cc', [Prob.c_L;1],[Prob.c_U;1E12]);
      EGOProb.xD          = x_D;
      EGOProb.xL          = x_L;
      %??????EGOProb.USER.dc     = 'rbf_dc';
      EGOProb.c           = Prob.USER.c;
      EGOProb.dc          = Prob.USER.dc;
      EGOProb.d2c         = Prob.USER.d2c;
      EGOProb.cNargin     = xnargin(Prob.USER.c);
      EGOProb.SCALE       = 1;
      if isempty(dc)
         EGOProb.dcNargin = 0;
      else
         EGOProb.dcNargin = xnargin(dc);
      end
   else
      EGOProb = glcAssign('ego_f',x_LL,x_UU,'EGO',A,b_L,b_U,'ego_c',1,1E12);
   end

elseif dCon == 0 & dLin == 0
   EGOProb = glcAssign('ego_f',x_LL,x_UU,'EGO',[],[],[],'ego_c',1,1E12);
else
   EGOProb = glcAssign('ego_f',x_LL,x_UU,'EGO', Prob.A, Prob.b_L, ...
                        Prob.b_U,'ego_cc', [Prob.c_L;1],[Prob.c_U;1E12]);
   EGOProb.SCALE         = 0;

   EGOProb.c           = Prob.USER.c;
   EGOProb.dc          = Prob.USER.dc;
   EGOProb.d2c         = Prob.USER.d2c;
   if dCon > 0
      EGOProb.cNargin     = xnargin(Prob.USER.c);
      if isempty(dc)
         EGOProb.dcNargin = 0;
      else
         EGOProb.dcNargin = xnargin(dc);
      end
   else
      EGOProb.cNargin     = 0;
      EGOProb.dcNargin    = 0;
   end
   %??????EGOProb.USER.dc     = 'rbf_dc';
end

% Send user info
if isfield(Prob,'user')
   EGOProb.user          = Prob.user;
end
EGOProb.MIP              = Prob.MIP;
EGOProb.uP               = Prob.uP;
EGOProb.P                = Prob.P;
EGOProb.mLin             = dLin;
EGOProb.mNonLin          = dCon+1;
EGOProb.EGO.k              = d;

optParam                  = optParamDef(globalSolver,solvType,d,dCon,dCon+dLin);
EGOProb.optParam          = optParam;
EGOProb.optParam.MaxIter  = GOMaxIter; 
EGOProb.optParam.MaxFunc  = GOMaxFunc;
EGOProb.optParam.epsRank  = epsRank;
EGOProb.optParam.IterPrint= DebugPriLev > 0;
EGOProb.optParam.bTol     = bTol;
EGOProb.optParam.cTol     = cTol;
EGOProb.GradTolg          = Prob.GradTolg;
EGOProb.GradTolH          = Prob.GradTolH;
EGOProb.GradTolJ          = Prob.GradTolJ;
EGOProb.GO                = Prob.GO;
EGOProb.SOL               = Prob.SOL;
if isfield(Prob,'OQNLP')
   EGOProb.OQNLP          = Prob.OQNLP;
end
if isfield(Prob,'KNITRO')
   EGOProb.KNITRO         = Prob.KNITRO;
end
EGOProb.CGO.idea          = idea;
EGOProb.CGO.epsRank       = epsRank;
EGOProb.PriLevOpt         = DebugPriLev;

% Take advantage of any derivatives in local search
dc                        = Prob.USER.dc;

EGOProb.USER.dc           = dc;

if isempty(dc)
   DerLvl_sn             = 0;
   DerLvl_gn             = 0;
   cDiff                 = 6;
else
   DerLvl_sn             = 3;
   DerLvl_gn             = 2;
   cDiff                 = 0;
end

DACEProb.GO.maxFunc1    = maxFunc1;
DACEProb.GO.maxFunc2    = maxFunc2;
DACEProb.GO.maxFunc3    = maxFunc3;
DACEProb.GO.DIRECT      = DIRECT;

switch lower(localSolver)
 case {'snopt','npsol','nlssol','minos'}
   EGOProb.NumDiff         = 6;
   EGOProb.SOL.optPar(39)  = DerLvl_gn;
   EGOProb.ConsDiff        = cDiff;

   DACEProb.NumDiff        = 6;
   DACEProb.SOL.optPar(39)= DerLvl_sn;
   DACEProb.ConsDiff       = cDiff;
   %DACEProb.optPar(10)     = 1E-8;
   DACEProb.optPar(12)     = 1E-8; %Minor optimality tolerance
 case {'glccluster'}
   if isempty(GOlocalSolver)
      DACEProb.GO.localSolver = localSolver;
   else
      DACEProb.GO.localSolver = GOlocalSolver;
   end
 case {'glcFast','glbFast'}
   snProb.NumDiff        = 0;
   snProb.ConsDiff       = 0;
 otherwise
   EGOProb.NumDiff         = 1;
   EGOProb.ConsDiff        = cDiff > 0;

   DACEProb.NumDiff        = 1;
   DACEProb.ConsDiff       = cDiff > 0;
end
EGOProb.GO.maxFunc1    = maxFunc1;
EGOProb.GO.maxFunc2    = maxFunc2;
EGOProb.GO.maxFunc3    = maxFunc3;
EGOProb.GO.DIRECT      = DIRECT;
switch lower(globalSolver)
 case {'glccluster'}
   LOCAL = 0;
   if isempty(GOlocalSolver)
      EGOProb.GO.localSolver  = localSolver;
   else
      EGOProb.GO.localSolver  = GOlocalSolver;
   end

 case {'glcFast','glbFast'}
   gnProb.NumDiff        = 0;
   gnProb.ConsDiff       = 0;
 otherwise
end

z = Fpen;
% Set infeasible points as infinity before check on minium
z(Fpen-F >= 1E-14)=Inf;

% Set integer infeasible points as infinity before check on minium
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

% Keep point and value in transformed space
% fIdx was previously yIdx in EGO code and documentation
yIdx = fIdx;

% EGO is using transformed values in variables
% y
% yNew
% yMin

y               = F;      % Used in fitmodel

if PriLev > 1 | IterPrint
   fprintf('Iter %3d n %3d ', 0, n);
   fprintf('nFunc %3d ', nFunc);
   fprintf('%d/%d %d:%d ', time([3 2 4 5]));
   %fprintf('Cycle %2d ', modN);
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
   %   fprintf('fMin %8.5f yNew %8.5f ', fMin, yNew);
   %end
   %fprintf('fE %8.5f ',fk);
   %fprintf('I%%');
   %fprintf(' %6.4f ',100*abs(EGOResult.f_k/fMin));
   fprintf('\n');
   xprint(O_min,'xMin:',' %12.8f',8)
end

LowExp = 0;
ExitFlag = 1;
ExitText = 'Maximal number of iterations reached';

%needed?
FLOWER    = fIdx(1);
fDiff     = Inf;
fDiffOld  = NaN;
NOUPDATE  = 0;
SAME1     = 0;
SAME2     = 0;
O_pre     = inf*ones(d,1);

% This can maybe be removed later on:
if isempty(IntVars)
   EGOProb  = ProbCheck(EGOProb,globalSolver,checkType('glc'));
   DACEProb = ProbCheck(DACEProb,localSolver,checkType('con'));
else
   EGOProb  = ProbCheck(EGOProb,globalSolver,checkType('minlp'));
   DACEProb = ProbCheck(DACEProb,localSolver,checkType('minlp'));
end


% *********************************************************
% *************** START MAIN ITERATION LOOP ***************
% *********************************************************

% n     = Number of points used in the interpolation
% nFunc = Total number of function evaluations
% nInit = Number of points used for startup, now or before
% Iter  = Number of fresh sampled points 

MaxIter = MaxFunc-nFunc;

convflag  = 0;
cpumax    = 0;
TIME0     = Prob.TIME0;

%while nFunc < MaxFunc & convflag == 0

for Iter = 1:MaxIter  % Main iteration loop
   
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end

   DACEProb.DACE.X = X;      % Used in dace_f 
   DACEProb.DACE.y = y;      % Used in dace_f   
   
%if DEBUG2 & d==2   
%   plot(X(1,:),X(2,:),'.r');
%   grid on
%   if Iter > 1
%      r_DB = zeros(n-1,1);
%      for i = 1:n-1
%        % r(i) = exp(-( sum( theta.*(abs( x-X(:,i) )).^p) ));
%        x_DB = [0;-1];
%        r_DB(i) = exp(-(  theta'*(abs( x_DB-X(:,i) )).^p ));
%      end
%      y_hat_DB = my + r_DB'*invR*(y(1:n-1)-my)
%   end
%end  

%if DEBUG2 & d==2  % contour plot of the Likelihood function
%   [xLL,xUU] = dace_bounds(X,y) % Compute variable bounds
%   u = linspace(xLL(1),xUU(1),50);
%   v = linspace(xLL(1),xUU(2),50);
%   for j=1:length(u)
%      for i=1:length(v)
%         x=[u(j);v(i)];
%         ff(i,j)=dace_f(x,Prob);
%      end
%   end
%   contour(u,v,ff,75);
%   save 'vars' u v ff 
%   Result=endSolve(Prob,Result);
%   return
%end

   % Maximize the Likelihood function. 
   if PriLev > 1
      disp('Maximize the Likelihood function');
   end
   [theta p my sigma2 invR fk] = fitmodel(...
          X, y, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
   
   % -----------------------------------------------------------------------
   %   STEP 2    % Cross-validation of the DACE model (first iteration only)
   % -----------------------------------------------------------------------
   
   if Iter == 1 & isempty(TRANSFORM)
      if PriLev > 0
         fprintf('Cross-validation of the DACE model\n');
      end
      
      [valid1,maxv1,sumv1] = crossval(X,y,my,sigma2,theta,p,epsRank);
      if all(F > 0)
         z = log(F);
         DACEProb.DACE.y = z;
         [theta p my sigma2 invR fk] = fitmodel( ...
                X, z, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
         [valid2,maxv2,sumv2] = crossval(X,z,my,sigma2,theta,p,epsRank);
         TRANSFORM = 1;
      elseif all(F < 0)
         z = -log(-F);
         DACEProb.DACE.y = z;
         [theta p my sigma2 invR fk] = fitmodel(...
                X, z, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
         [valid2,maxv2,sumv2] = crossval(X,z,my,sigma2,theta,p,epsRank);
         TRANSFORM = 2;
      else
         TRANSFORM = 0;
         DACEProb.DACE.y = F;
         valid2 = 0; maxv2 = inf; sumv2 = inf;
      end
      if all(F ~= 0)
         z = -1./F;
         DACEProb.DACE.y = z;
         [theta p my sigma2 invR fk] = fitmodel(...
                X, z, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
         [valid3,maxv3,sumv3] = crossval(X,z,my,sigma2,theta,p,epsRank);
      else
         valid3 = 0; maxv3 = inf; sumv3 = inf;
      end
      %if 0
      %   if maxv1 <= min(maxv2,maxv3) 
      %      TRANSFORM = 0;
      %   elseif maxv2 <= maxv3
      %      if TRANSFORM == 1
      %         y = log(F);
      %      else
      %         y = -log(-F);
      %      end
      %   else
      %      TRANSFORM = 3;
      %      y = -1./F;
      %   end
      %else
         if sumv1 <= min(sumv2,sumv3) | isnan(sumv1)
            TRANSFORM = 0;
            y = F;
         elseif sumv2 <= sumv3
            if TRANSFORM == 1
               y = log(F);
            else
               y = -log(-F);
            end
         else
            TRANSFORM = 3;
            y = -1./F;
         end
      %end
      if TRANSFORM < 3
         % Must recompute
         [theta p my sigma2 invR fk] = fitmodel(...
                X, y, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
      end
      if PriLev > 0
         fprintf('EGO - TRANSFORM = %d\n',TRANSFORM); 
      end
      yMin                = y(fIdx);
   elseif Iter == 1 & ~isempty(TRANSFORM)

      TRANSFORM = min(3,max(0,TRANSFORM));
      if TRANSFORM == 0  % No transformation made.
         y = F;
      elseif TRANSFORM == 1  % log(y) transformation made.
         y = log(F);
      elseif TRANSFORM == 2  % -log(-y) transformation made.
         y = -log(-F);
      elseif TRANSFORM == 3  % -1/y transformation made.
         y = -1./F;
      end
      if TRANSFORM > 0
         % Must recompute
         [theta p my sigma2 invR fk] = fitmodel(...
                X, y, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
      end
      yMin                = y(fIdx);
   
   end
   %if 0 & Iter == 1
   %   if PriLev > 0
   %      fprintf('Cross-validation of the DACE model\n');
   %   end
   %   
   %   VALID = crossval(X,y,my,sigma2,theta,p,epsRank);
   
   %   if VALID  % DACE-model valid
   %      if PriLev > 0
   %         fprintf('DACE model valid without transformation of y.\n')
   %      end
   %      TRANSFORM = 0; % No transformation of y      
   %   else
   %      if all(y > 0)
   %         % Refit the DACE model after applying a log(y) transform
   %         if PriLev > 0
   %            fprintf('DACE model not valid. ')
   %            fprintf('Trying log(y) transformation.\n')
   %        end
   %        z = log(y);
   %        DACEProb.DACE.y = z;
   %        [theta p my sigma2 invR fk] = fitmodel(...
   %               X, z, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
   %        VALID = crossval(X,z,my,sigma2,theta,p,epsRank);
   %        if VALID
   %           if PriLev > 0
   %              fprintf('DACE model valid after log(y) transformation.\n')
   %           end
   %           y = z;
   %           TRANSFORM = 1; % log(y) transformation
   %        end
   %     elseif all(y < 0)
   %        % Refit the DACE model after applying a -log(-y) transform
   %        if PriLev > 0
   %           fprintf('DACE model not valid. \n')
   %           fprintf('Trying -log(-y) transformation.\n')
   %        end
   %        z = -log(-y);
   %        DACEProb.DACE.y = z;
   %        [theta p my sigma2 invR fk] = fitmodel(...
   %               X, z, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
   %        VALID = crossval(X,z,my,sigma2,theta,p,epsRank);
   %        if VALID
   %           if PriLev > 0
   %              fprintf('DACE model valid after -log(-y) transformation.')
   %              fprintf('\n')
   %           end
   %           y = z;
   %           TRANSFORM = 2; % -log(-y) transformation
   %        end
   %     end   
   %     if ~VALID         
   %        % Refit the DACE model after applying a -1/y transform
   %        if PriLev > 0
   %           fprintf('DACE model not valid after log transformation. ')
   %           fprintf('Trying -1/y transformation.\n')
   %        end
   %        z = -1./y;
   %        DACEProb.DACE.y = z;
   %        [theta p my sigma2 invR fk] = fitmodel(...
   %               X, z, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
   %        VALID = crossval(X,z,my,sigma2,theta,p,epsRank);
   %        if VALID
   %           if PriLev > 0
   %              fprintf('DACE model valid after -1/y transformation\n')
   %           end
   %           y = z;
   %           TRANSFORM = 3; % -1/y transformation
   %        else
   %           if PriLev > 0
   %              fprintf('DACE model not valid after -1/y transformation. '); 
   %              %fprintf('Problem skipped!\n')
   %           end
   %           if 0
   %           [yMin fIdx] = min(y);
   %           Result.f_k      = yMin;
   %           Result.x_k      = X(:,fIdx);
   %           Result.FuncEv   = n;
   %           Result.Iter     = Iter;
   %           Result.ExitFlag = 1;
   %           Result.ExitText = ...
   %                 'DACE model not valid after -1/y transformation';
 
   %           Result=endSolve(Prob,Result);
   %           return;
   %           end
   %           disp('TRY THIS MODEL NEVERTHELESS')
   %           y = z;
   %           TRANSFORM = 3; % -1/y transformation
   %        end
   %     end
   %  end 
   %end


   % -----------------------------------------------------------------------
   %   STEP 3    % Maximization of the expected improvement 
   % -----------------------------------------------------------------------
   if PriLev > 1
      disp('Maximization of the expected improvement')
   end
   
   % Used in ego_f:  
   EGOProb.EGO.my      = my;     
   EGOProb.EGO.theta   = theta;
   EGOProb.EGO.p       = p;
   EGOProb.EGO.y       = y;
   EGOProb.EGO.X       = X;
   EGOProb.EGO.sigma2  = sigma2;
   EGOProb.EGO.invR    = invR;
   EGOProb.EGO.yMin    = yMin;
   EGOProb.Name        = ['EGO problem ' num2str(Iter)] ;  
   
   if DEBUG2 & d==2  % contour plot of the Expected improvement
      u = linspace(x_L(1),x_U(1),50);
      v = linspace(x_L(2),x_U(2),50);
      for j=1:length(u)
         j
         for i=1:length(v)
            x=[u(j);v(i)];
            ff(i,j)=ego_f(x,EGOProb);
         end
      end
      contour(u,v,ff,75);
      save('EGOvars.mat','u','v','ff') 
      pause(2)
   end

      EGOResult                = tomRun(globalSolver,EGOProb,PriLev-4);
      ExitFlag = EGOResult.ExitFlag;
      %if isempty(EGOResult.x_k)
      if ExitFlag == 7
         % Infeasible solution, try warm start to get feasible
         EGOProb.WarmStart = 1;
         % keyboard
         if IterPrint | PriLev > 0
            fprintf('No feasible point, ');
            fprintf('Warm Start %s',globalSolver);
            fprintf(' with MaxFunc %d',EGOProb.optParam.MaxFunc);
            fprintf(' and MaxIter %d\n',EGOProb.optParam.MaxIter);
            if strcmpi(globalSolver,'glcCluster')
               fprintf('maxFunc1 %d ',EGOProb.GO.maxFunc1);
               fprintf('maxFunc2 %d ',EGOProb.GO.maxFunc2);
               fprintf('maxFunc3 %d ',EGOProb.GO.maxFunc3);
               fprintf('\n');
            end
         end
         %EGOResult                = tomRun(globalSolver,EGOProb,PriLev-4);
         EGOResult                = tomRun(globalSolver,EGOProb,0);
         PrintResult(EGOResult,double(PriLev > 0))
         EGOProb.WarmStart = 0;
      end
      if isempty(EGOResult.x_k)
         % Should not occur
         disp('Failure in global solver. Could be severe ill-conditioning')
         disp('(or difficult nonlinear constraints).')
         disp('Try a local solution instead)')

         EGOProb.optParam.MaxIter = 10000; 
         if isinf(fMin)
            EGOProb.x_0 = 0.5*(x_UU-x_LL);
         else
            EGOProb.x_0 = x_min;
         end
         EGOResult                = tomRun(localSolver,EGOProb,PriLev-4);
      end
      xNew = EGOResult.x_k;
      f_k  = EGOResult.f_k;
      if LOCAL
         LocSteps = min(20,size(xNew,2));
         if LocSteps > 1 & PriLev > 0 
            fprintf('Do %d local search steps ',LocSteps);
            if size(xNew,2) > LocSteps
               fprintf('out of %d\n',size(xNew,2));
            end
            fprintf('\n');
         end
         % Do local search from globally best points
         EGOProb.snP  = Iter;
         % Max 1000 iterations in local solver
         EGOProb.optParam.MaxFunc  = MaxIter*max(d,10); 
         EGOProb.optParam.MaxIter  = MaxIter; 

         %snProb.CGO  = EGOProb.CGO;

         fNew = inf;
         iBest = 0;
         xBest = [];
         for i = 1:LocSteps
             EGOProb.x_0  = xNew(:,i);
             if PriLev == 4
                xprint(EGOProb.x_0,'x0:   ');
             end
             if PriLev > 4
                fprintf('Try local search #%d\n',i);
             end
             gnR = tomRun(localSolver,EGOProb,max(PriLev-4,DebugPriLev*2));
             OK  = 1;
             if dLin > 0
                % Check linear constraints
                Ax = EGOProb.A*gnR.x_k;
                AxMax = max(abs(Ax));
                if any(EGOProb.b_L - Ax > max(1E-5,1E-5*AxMax)  ...
                     | Ax - EGOProb.b_U > max(1E-5,1E-5*AxMax))
                   %'DO NOT ACCEPT LOCAL POINT'
                   %'LINEAR CONSTRAINTS NOT FULFILLED'
                   %ExitFlag = gnR.ExitFlag
                   %Inform   = gnR.Inform
                   OK = 0;
                end
             end
             if gnR.f_k < fNew  & OK
                iBest     = i;
                EGOResult = gnR;
                xBest     = gnR.x_k;
                fNew      = gnR.f_k;
                if PriLev == 4
                   xprint(xBest,'xOpt: ');
                end
             end
         end
         xNew        = xBest;

         if PriLev > 2
            fprintf('Global f(x) %30.20f       with %s\n',f_k,globalSolver);
            fprintf('Local  f(x) %30.20f (#%2d) with %s\n',...
                     fNew,iBest,localSolver);
         end
      else
         xNew       = xNew(:,1);
         fNew       = f_k;
         if PriLev > 2
            fprintf('Global f(x) %30.20f with %s\n',f_k,globalSolver);
         end
      end
      if PriLev > 3
         xprint(xNew,'x: ');
      end
   
   ExpI = -EGOResult.f_k(1);  % Optimal value of expected improvement
   
   % Plot result of Maximization of the expected improvement
   if DEBUG2 & d==2
      % Must read result from saved file to get C
      CC=load('glbSave.mat','C');
      % Transform from [0,1] to original coordinates
      C = tomsol(9, x_L, CC.C, x_U-x_L); 
      plot(C(1,:),C(2,:),'.');
      hold on;
      plot(xNew(1),xNew(2),'*y');
      pause(2)
      hold off
   end

   % *************** UPDATE ***************
   % Remove information for response surface problem

   % New point in original space
   if SCALE
      O_new   = tomsol(9, x_L, xNew, x_D); 
   else
      O_new   = xNew;
   end

   % Clean Tomlab global variables, because switch of problem
   global NLP_x NLP_f NLP_c NLP_xc 
   NLP_x=[]; NLP_f=[]; NLP_c=[]; NLP_xc=[];
   NARG = NARGSAVE;

   if isempty(O_new) % Infeasibility problem
      fNew    = Inf;
      fPen    = Inf;
      SAME1   = 0;
      SAME2   = SAME2+1;
      Update  = -1;
      if size(X,2) > d+1
         % It is possible to remove a point from X
         minDist = -1;
      else
         minDist = 0;
      end
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
         Update = -1;
         minDist = 0;
         fNew   = NaN;
      end
   end

   time  = fix(clock);

   %if REPLACE
   %   fnewV=min(median([F;fNew]),[F;fNew]);
   %   Update = tomsol(24, xNew, fnewV);    
   %else
   %   Update = tomsol(24, xNew, fNew);    
   %end

   %Update = 1;  % Currently no check if points are too close

   if Update == 1
      Dist = tomsol(30,xNew,X);
      [minDist ixDist] = min(Dist);
      if PriLev > 0
         fprintf('Minimal distance to X set');
         fprintf(' %s',minDist);
         fprintf('\n');
      end
   end
     
   if minDist == 0
      if PriLev > -1
         fprintf('Convergence: New point identical to previous point');
         fprintf(' fMin %f,',fMin);
         fprintf(' %d iterations\n',Iter);
         fprintf('Maybe global optimum found?');
      end
      ExitFlag = 0;
      ExitText = 'Algorithm generates same points, maybe global optimum found';
      break;
   end
   if minDist > 1E-5
      F   = [F;fNew];
      Fpen= [Fpen;fPen];
      X   = [X xNew];
      O   = [O O_new];
      O_pre = O_new;
      n   = n + 1;
      Update = 1;
   elseif Update == 1   %if F(ixDist) > fNew
      % Always add new point, but delete close point
      F(ixDist)   = fNew;
      Fpen(ixDist)= fPen;
      X(:,ixDist) = xNew;
      O(:,ixDist) = O_new;
      Update = 2;
   else
      % Infeasible problem, ill-conditioning?
      % Remove most infeasible point
      [maxPen ixDist] = max(Fpen-F);
      if maxPen == 0
         % If all points feasible, remove point with largest f(x)
         [maxPen ixDist] = max(Fpen);
      end
      ix =  [1:ixDist-1,ixDist+1:n];
      F    = F(ix);
      Fpen = Fpen(ix);
      X    = X(:,ix);
      O    = O(:,ix);
      n   = n - 1;
      if ixDist < fIdx
         fIdx = fIdx -1;
      end
      Update = 3;
   end
      
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
      fIdx     = length(Fpen);
      x_min    = xNew;
      O_min    = O_new;
      FLOWER   = nFunc;
   end
   NOUPDATE = 0;

      %if yMin >= 0
      %   yPred = (1-ExpI)*yMin;
      %else
      %   yPred = (1+ExpI)*yMin;
      %end

   yPred = yMin-ExpI;

   if TRANSFORM == 0  % No transformation made.
      yNew = fNew;
      fPred = yPred;
   elseif TRANSFORM == 1  % log(y) transformation made.
      if fNew == 0
         fGoal = 0; % Stop iterations
         yNew  = 1E20;
         fPred = 1E20;
      else
         yNew  = log(fNew);
         fPred = exp(yPred);
      end
   elseif TRANSFORM == 2  % -log(-y) transformation made.
      if fNew == 0
         fGoal = 0; % Stop iterations
         yNew  = 1E20;
         fPred = 1E20;
      else
         yNew  = -log(-fNew);
         fPred = -exp(-yPred);
      end
   elseif TRANSFORM == 3  % -1/y transformation made.
      if fNew == 0
         fGoal = 0; % Stop iterations
         yNew  = 1E20;
         fPred = 1E20;
      else
         yNew  = -1/fNew;
         fPred = -1/yPred;
      end
   end
   if Update == 1
      y(n)      = yNew;
   elseif Update == 2
      y(ixDist) = yNew;
   else
      y         = y(ix);
   end
   yMin = y(fIdx);

   %if PriLev > 1
   %   fprintf('\n\n Relative expected improvement : %e\n',...
   %       100*abs(EGOResult.f_k/yMin));
   %   fprintf('\n min(y) = %f ',yMin);
   %end
   
   %xBest = X(:,fIdx);
   
   if PriLev > 1 | IterPrint
      fprintf('Iter %3d n %3d ', Iter, n);
      fprintf('nFunc %3d ', nFunc);
      fprintf('%d/%d %d:%d ', time([3 2 4 5]));
      if ~isempty(fGoal) & ~isinf(fGoal)
         fprintf(' fGoal %8.5f', fGoal);
      end
      if Feasible
         fprintf(' fMinF %11.8f ', fMin);
      else
         fprintf(' fMinI %11.8f ', fMin);
      end
      fprintf('at %3d fNew %11.8f ', FLOWER, fNew);
      %NOTfprintf(' fLoc %11.8f', min_sn);
      fprintf('fPred %8.5f ', fPred);

      %NEWHKHfprintf(' %11.8f', fGlob);
      %if TRANSFORM ~= 0
      %   fprintf('yMin %8.5f yNew %8.5f ', yMin, yNew);
      %end
      %fprintf('fE %8.5f ',fk);
      %fprintf('I%%');
      %fprintf(' %6.4f ',100*abs(EGOResult.f_k/yMin));

      fprintf('Best %d ', fIdx);
      fprintf('fE %8.5f ',fk);
      fprintf('I%%');
      fprintf(' %6.4f ',100*ExpI);
      fprintf('Q%%');
      if yMin ~= 0, fprintf(' %6.4f ',100*abs(ExpI/yMin)); end
      fprintf('\n');
      if TRANSFORM ~= 0
         fprintf('yMin %8.5f yNew %8.5f ', yMin, yNew);
         fprintf('yPred %8.5f ', yPred);
         fprintf('\n');
      end
      fprintf('Cond invR %e ',cond(invR));
      fprintf('\n');
      xprint(xNew,'xNew:',' %12.8f',8)
      if SCALE > 0
         xprint(O_new,'oNew:',' %12.8f',8)
      end
   end
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
      plot(O_new(1),fNew,'*g');
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
   convflag = isClose(fGoal,fMin,fTol,nFunc,Iter,PriLev);
   if convflag
      if PriLev > -1
         fprintf('Convergence: fMin %f close to optimal fOpt %f ',fMin,fGoal);
         fprintf('after %d iterations\n',Iter);
      end
      ExitFlag = 0;
      ExitText = 'Close to user given optimal function value';
      break;
   end
   if abs(ExpI) < TolExpI & ( ...
      (fPred >= 0.9999*fMin & fMin > 0) | (abs(yPred) < 1E-7 & fMin == 0)|...
      (abs(fPred) >= 1.0001*abs(fMin) & fMin < 0))
      %(fPred >= 0.9999*fMin & fMin > 0) | (abs(yPred) < 1E-7 & fMin == 0)|...
      %(abs(fPred) >= 1.0001*abs(fMin) & fMin < 0))
      LowExp = LowExp + 1;
      if LowExp >= 3
      if PriLev > -1
         fprintf('Convergence: ');
         fprintf('ExpI %e ',ExpI);
         if fMin >= 0
               fprintf(' yPred %f <= 0.999*fMin %f',...
                 fPred,0.999*fMin);
         else
               fprintf(' |yPred| %f <= 1.001*|fMin| %f',...
                 fPred,1.001*fMin);
         end
         fprintf('. Iterations %d\n',Iter);
      end
      ExitFlag = 0;
      ExitText = 'Expected improvement low for three iterations';
      break;
      end
   else
      LowExp = 0;
   end
   %if abs(ExpI) <= abs(0.001*fMin) | (abs(ExpI) < 0.0001 & fMin == 0)
   %   if PriLev > -1
   %      fprintf('Convergence: Expected Improvement %f <= 0.001*fMin %f',...
   %              abs(ExpI),0.001*fMin);
   %      fprintf('. Iterations %d\n',Iter);
   %   end
   %   break;
   %end
   if fk > 1E8
      fprintf('Failure! Ill-conditioning in covariance matrix. ');
      fprintf('Cannot proceed!\n');
      ExitFlag = 3;
      ExitText = 'Cannot proceed! Ill-conditioning in covariance matrix.';
      break;
   end
   % Use previous optimal point as starting point
   %DACEProb.x_0 = DACEResult.x_k;

   if pEst
      DACEProb.x_0 = [log(theta);p];
   else
      DACEProb.x_0 = log(theta);
   end
   
   
end % Main iteration loop

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
   fprintf('constraints are violated')
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
else
   Result.Inform   = convflag;
end
Result          = endSolve(Prob,Result);

%tomsol(25) % Deallocates memory

%Result.x_k       = xBest;
%Result.ExitFlag  = ExitFlag;
%Result.ExitText  = ExitText;


%-----------------------------------------------------------
% function [theta, p, my, sigma2, invR, fk] = fitmodel(...
%           X, y, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);
%-----------------------------------------------------------

function [theta, p, my, sigma2, invR, fk] = fitmodel(...
          X, y, DACEProb, Iter, epsRank, pEst, localSolver,PriLev);

%X            = DACEProb.DACE.X;  % Sample points
%y            = DACEProb.DACE.y;  % Function values at sample points

[d,n]        = size(X);

[xLL,xUU]    = dace_bounds(X,y); % Compute variable bounds

DACEProb.Name= ['DACE problem ' num2str(Iter)];

if pEst
   DACEProb.x_L = [xLL';ones(d,1)];
   DACEProb.x_U = [xUU';2*ones(d,1)];

else
   DACEProb.x_L = xLL';
   DACEProb.x_U = xUU';
end

global GOMaxFunc GOMaxIter globalSolver maxFunc1 maxFunc2 maxFunc3
%if 1 
   [DACEResult] = tomRun(localSolver, DACEProb, PriLev-3);
   %DACEProb.PriLevOpt = 11;
   %[DACEResult] = tomRun(localSolver, DACEProb, 2);
   %pause

   %DACEResult = tomSolve(localSolver, DACEProb);
   ExitFlag = DACEResult.ExitFlag;

   fk    = DACEResult.f_k;
   if ExitFlag == 4 | ExitFlag == 7
      if strcmpi(localSolver,'glcCluster') | strcmpi(localSolver,'glcFast')
         fprintf('Infeasible. TRY warm start %s\n',localSolver);
         if PriLev > 3
            PrintResult(DACEResult,2);
         end
         % keyboard
         DACEProb.WarmStart = 1;
         if PriLev > 0
            fprintf('DACE with Global solver: No feasible point, ');
            fprintf('Warm Start %s',localSolver);
            fprintf(' with MaxFunc %d',DACEProb.optParam.MaxFunc);
            fprintf(' and MaxIter %d\n',DACEProb.optParam.MaxIter);
            if strcmpi(globalSolver,'glcCluster')
               fprintf('maxFunc1 %d ',DACEProb.GO.maxFunc1);
               fprintf('maxFunc2 %d ',DACEProb.GO.maxFunc2);
               fprintf('maxFunc3 %d ',DACEProb.GO.maxFunc3);
               fprintf('\n');
            end
         end
         % DACEResult = tomRun(localSolver, DACEProb,PriLev-2);
         DACEResult = tomRun(localSolver, DACEProb);
         PrintResult(DACEResult,double(PriLev > 0))
         DACEProb.WarmStart = 0;
      elseif pEst
         if PriLev > 3
            PrintResult(DACEResult,2);
         end
         if any(DACEResult.x_k(d+1:d+d) == 2)
            DACEProb.x_U = [DACEProb.x_U(1:d);ones(d,1)];
         else
            DACEProb.x_L = [DACEProb.x_L(1:d);2*ones(d,1)];
         end
         DACEProb.x_0 = [];
         %DACEResult = tomRun(localSolver, DACEProb,PriLev-2);
         DACEResult = tomRun(globalSolver, DACEProb);
         PrintResult(DACEResult,double(PriLev > 0))
      else
         DACEResult = tomRun(globalSolver, DACEProb);
         PrintResult(DACEResult,double(PriLev > 0))
      end
      %if isinf(fk)
      %   %if GOMaxFunc >= 2500
      %   %   DACEProb.optParam.MaxFunc = GOMaxFunc;
      %   %   DACEProb.optParam.MaxIter = GOMaxIter;
      %   %else
      %   %   DACEProb.optParam.MaxFunc = 2*GOMaxFunc;
      %   %   DACEProb.optParam.MaxIter = 2*GOMaxIter;
      %   %end
      %   %if strcmpi(globalSolver,'glcCluster')
      %   %   DACEProb.GO.maxFunc1    = 2*maxFunc1;
      %   %end
      %   DACEProb.WarmStart = 1;
      %   %if IterPrint | PriLev > 0
      %   if PriLev > 0
      %      fprintf('Local min_sn with global solver: No feasible point, ');
      %      fprintf('Warm Start %s',localSolver);
      %      fprintf(' with MaxFunc %d',DACEProb.optParam.MaxFunc);
      %      fprintf(' and MaxIter %d\n',DACEProb.optParam.MaxIter);
      %      if strcmpi(globalSolver,'glcCluster')
      %         fprintf('maxFunc1 %d ',gnProb.GO.maxFunc1);
      %         fprintf('maxFunc2 %d ',gnProb.GO.maxFunc2);
      %         fprintf('maxFunc3 %d ',gnProb.GO.maxFunc3);
      %         fprintf('\n');
      %      end
      %   end
      %   % DACEResult = tomRun(localSolver, DACEProb,PriLev-2);
      %   DACEResult = tomRun(localSolver, DACEProb);
      %   PrintResult(DACEResult,double(PriLev > 0))
      %   DACEProb.optParam.MaxFunc = GOMaxFunc;
      %   DACEProb.optParam.MaxIter = GOMaxIter;
      %   DACEProb.GO.maxFunc1      = maxFunc1;
      %   DACEProb.WarmStart = 0;
      %else
      %   if any(DACEResult.x_k(d+1:d+d) == 2)
      %      DACEProb.x_U = [DACEProb.x_U(1:d);ones(d,1)];
      %   else
      %      DACEProb.x_L = [DACEProb.x_L(1:d);2*ones(d,1)];
      %   end
      %   DACEProb.x_0 = [];
      %   %DACEResult = tomRun(localSolver, DACEProb,PriLev-2);
      %   DACEResult = tomRun(localSolver, DACEProb);
      %   PrintResult(DACEResult,double(PriLev > 0))
      %end

      %zz=DACEProb.optParam.MaxFunc;
      %DACEProb.optParam.MaxFunc = 20000;
      %DACEResult = tomRun('glcFast', DACEProb);
      %PrintResult(DACEResult,1);
      %DACEProb.optParam.MaxFunc = zz;

      fk    = DACEResult.f_k;
   elseif fk > 1E8 & pEst
      fprintf('TRY ANOTHER CALL TO %s\n',localSolver);
      if PriLev > 3
         PrintResult(DACEResult,2);
      end
      if any(DACEResult.x_k(d+1:d+d) == 2)
         DACEProb.x_U = [DACEProb.x_U(1:d);ones(d,1)];
      else
         DACEProb.x_L = [DACEProb.x_L(1:d);2*ones(d,1)];
      end
      DACEProb.x_0 = [];
      DACEResult = tomRun(localSolver, DACEProb,PriLev-2);
      fk    = DACEResult.f_k;
   end
% OLD CODE
%elseif 1 
%   DACEProb.Solver.Alg = 2;
%   DACEProb.NumDiff    = 1;
%   %[DACEResult] = tomRun('ucSolve', DACEProb, PriLev-3);
%   %[DACEResult] = tomSolve('ucSolve',DACEProb);
%   DACEResult = tomRun(localSolver, DACEProb,PriLev-2);
%   fk    = DACEResult.f_k;
%else % Use glbSolve (glbFast)
%   DACEProb.x_L = floor(xLL);
%   DACEProb.x_U = floor(xUU)+1;
%   DACEProb.optParam.MaxIter = 50; % Iterations in glbSolve
%   SolverGLB = 'glbFast';
%   %[DACEResult] = tomRun(SolverGLB, DACEProb, [], [], 'dace_prob',PriLev-3);
%   %[DACEResult] = tomSolve(SolverGLB,DACEProb);
%   DACEResult = tomRun(SolverGLB, DACEProb,PriLev-2);
%   if d==2
%      load('glbSave.mat','C');
%      % Transform from [0,1] to original coordinates
%      C = tomsol(9,DACEProb.x_L, C,DACEProb.x_U-DACEProb.x_L); 
%      plot(C(1,:),C(2,:),'.g');
%      pause
%   end
%   fk    = DACEResult.f_k;
%end  

alpha = DACEResult.x_k(1:d,1);
if fk > 1E8
   PrintResult(DACEResult,3);
end

if size(DACEResult.x_k,1) > d
   p = DACEResult.x_k(d+1:d+d,1);
   if PriLev > 0
      xprint(p,'p:  ');
   end
else
   p = 2;
end

theta = exp(alpha);

%if 0
%R = eye(n); % Compute the correlation matrix R
%if length(p) > 1
%   for i = 1:n-1
%      R(i,i+1:n) = exp(-( theta'*(abs( ...
%          X(:,i)*ones(1,n-i)-X(:,i+1:n) ).^(p*ones(1,n-i)) )));
%      R(i+1:n,i) = R(i,i+1:n)';
%   end
%else
%   for i = 1:n-1
%      R(i,i+1:n) = exp(-( theta'*(abs( X(:,i)*ones(1,n-i)-X(:,i+1:n) ).^p) ));
%      R(i+1:n,i) = R(i,i+1:n)';
%      %for j = i+1:n
%      %   R(i,j) = exp(-(  theta'*(abs( X(:,i)-X(:,j) ).^p)  ));
%      %   R(j,i) = R(i,j);
%      %end
%   end
%end
%end

R = tomsol(28,alpha,X,p);
if R(1,1) > 1E99
   my     = Inf;
   sigma2 = Inf;
   invR   = Inf; % Wrong size, but what to do if this occurs?
   'EGO, fitmodel: Should not occur'
   keyboard
end

[invR, detR, pRank] = tomsol(12, R, epsRank);
if PriLev > 0
   fprintf('pRank %d detR %e \n',pRank,detR);
end

%sum1 = sum(sum(abs(invR*R-eye(n))));
%sum2 = sum(sum(abs(R*invR-eye(n))));

%if (sum1 > 1E-9 | sum2 > 1E-9 ) & detR > 1E-15
%   sum1
%   sum2
%
%   keyboard
%end

% Estimate of mean my: (1'*invR*y)/(1'*invR*1);

my = sum(invR*y)/sum(sum(invR));

% Estimate of the variance s2 = ((y-my*1)'*invR*(y-my*1))/n;

sigma2 = max(1E-300,((y-my)'*invR*(y-my))/n);

%if sigma2 <= 1E-300
%   sigma2
%   keyboard
%end

%if sigma2 < 0
%   fprintf('Negative variance term sigma2 %15.7f\n',sigma2);
%   disp('sigma2 < 0 ');
%   keyboard
%   pause
%end


%-----------------------------------------------------------
%  function [VALID,maxv,sumv] = crossval(X,y,my,sigma2,theta,p,epsRank); 
%-----------------------------------------------------------

function [VALID, maxv,sumv] = crossval(X,y,my,sigma2,theta,p,epsRank); 

PriLev = 0;
if PriLev > 0
   disp('Subfunction CROSSVAL called');
end
n = length(y);
for i = 1:n
   x   = X(:,i);
   XX = [ X(:,1:i-1) X(:,i+1:n) ];
   y_i = [ y(1:i-1);y(i+1:n) ];
   R = zeros(n-1,n-1);
 
   if length(p) > 1
      for ii = 1:n-1
         R(ii,:) = exp(-( theta'*(abs( ...
                   XX(:,ii)*ones(1,n-1)-XX ).^(p*ones(1,n-1))) ));
      end
   else
      for ii = 1:n-1
         R(ii,:) = exp(-( theta'*(abs( XX(:,ii)*ones(1,n-1)-XX ).^p) ));
      end
   end
   %for ii = 1:n-1
   %   for jj = 1:n-1
   %      R(ii,jj) = exp(-( sum( theta.*(abs( XX(:,ii)-XX(:,jj) )).^p) ));
   %   end
   %end

   %[U S V] = svd(R);
   %S_inv = zeros(n-1);
   %S_inv(1,1) = 1/S(1,1);
   %for ii = 2:n-1
   %   if S(ii,ii) > epsRank*S(1,1)
   %      S_inv(ii,ii) = 1/S(ii,ii);
   %   end
   %end
   %invR = V * (S_inv * (U' ));

[invR, detR, pRank] = tomsol(12, R, epsRank);

%sum1 = sum(sum(abs(invR*R-eye(n-1))));
%sum2 = sum(sum(abs(R*invR-eye(n-1))));

%if (sum1 > 1E-9 | sum2 > 1E-9 ) & detR > 1E-15
%   sum1
%   sum2
%
%   keyboard
%end

   if length(p) > 1
      r = exp(-( theta'*(abs( x*ones(1,n-1)-XX ).^(p*ones(1,n-1))) ))';
   else
      r = exp(-( theta'*(abs( x*ones(1,n-1)-XX ).^p) ))';
   end

   %r = zeros(n-1,1);
   %for kk = 1:n-1
   %   r(kk) = exp(-( sum( theta.*(abs( x-XX(:,kk) )).^p) ));
   %end

   Rr = invR*r;

   %f(i) = my + r'*invR*(y_i-my);
   yHat = my + Rr'*(y_i-my);
   
   s = sqrt( sigma2*( 1-r'*Rr + (1-sum(Rr))^2/sum(sum(invR))  ) );
   
   value(i)=(y(i)-yHat)/s;
end % for

if any(abs(value)>3)
   VALID = 0;
else
   VALID = 1;
end
maxv = max(abs(value));
sumv = sum(abs(value));
if PriLev > 0
   fprintf('Standardized cross-validation residual maximum %f\n',maxv);
   fprintf('Standardized cross-validation residual sum %f\n',sumv);
end


% -------------
% dace_bounds.m
% -------------

function [x_L,x_U] = dace_bounds(X,y);

[k,n] = size(X);

for i = 1:k
   xL     =  min(X(i,:));
   xU     =  max(X(i,:));
   r      =  xU - xL;
   denom  =  max(r,r^2);
   x_L(i) =  log(-log(0.99)/max(denom,1E-30));
   r      =  0.05*r;
   denom  =  min(r,r^2);
   x_U(i) =  log(-log(0.01)/max(denom,1E-30));
end


%************************************
% The rest is just a saving of code
%************************************

% % contour plot of the Likelihood function
%   
% u = linspace(1E-5,0.05,50);
% v = linspace(1E-5,0.005,50);
% for j=1:length(u)
%    j
%    for i=1:length(v)
%       x=[u(j);v(i)];
%       ff(i,j)=DACE_f(x,DACEProb);
%       %ff(i,j)=feval(DACEProb.USER.f,x,DACEProb);
%   end
%end
%contour(u,v,ff,50);
%save 'vars' u v ff 
%end
%
%% contour plot of model function
%
%o = ones(n-1,1);   
%y = y(1:length(y)-1);
%u = linspace(x_L(1),x_U(1),100);
%v = linspace(x_L(2),x_U(2),100);
%%for j=1:length(u)
%   for i=1:length(v)
%      x=[u(j);v(i)];
%      r = zeros(n-1,1);
%      for kk = 1:n-1
%         r(kk) = exp(-( sum( theta.*(abs( x-X(:,kk) )).^p) ));
%      end
%      %ff(i,j)=nlp_f(x,Prob, varargin{:});
%      ff(i,j) = my_hat + r'*invR*(y-my_hat*o);
%   end
%end
%contour(u,v,ff,75);
%end
%
%plot(DACEResult.GLOBAL.C(1,:),DACEResult.GLOBAL.C(2,:),'*');
%pause(2)

%alpha = DACEResult.x_k(1:k,1);
%theta  = exp(alpha);

% ====================================================================
function convflag = isClose(fGoal,f,fTol,nFunc,Iter,PriLev)
% ====================================================================

convflag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convflag = 1;
elseif fGoal == 0
   if abs(f-fGoal) < fTol
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
    [t,p]=size(X);
    while min(sqrt(sum((x_n*ones(1,p)-X).^2)))<minDist
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
    X(j,:) = x_L(j) + floor(rand(1,n)*(x_U(j)-x_L(j)+1));
end


% MODIFICATION LOG
%
% 981127  mbk  DACEProb is defined without dace_prob.m.
%              DACEProb.x_0 is set in glb_prob.
% 981127  hkh  Add calls to ini/endSolve, globalSave/Get.
% 981208  hkh  Init of yMin = Inf; xBest= Prob.x_0;
%              Add Result.ExitFlag. Fix error return of Result.
% 981208  mbk  Init of yMin = Inf; xBest= Prob.x_0; ABORTED
%              Return the smallest function value of those initially
%              computed and corresponding x-value.
% 001118  hkh  Major revision
% 010304  hkh  Use tomSolve
% 010718  hkh  Use glbFast, revised similar to glbSolve
% 010726  hkh  Change local / solver handling and printing.
% 010903  hkh  General initialization of optPar vector to -999.
% 020103  hkh  Major revision taking many elements from rbfSolve
% 020616  hkh  Use input TolExpI and TRANSFORM. 
% 020616  hkh  Fix bug calling SOL solvers, estimating internal derivatives
% 020820  hkh  Print sub problem results dependent on PriLev
% 020820  hkh  Try local search for all possible global solutions
% 020823  hkh  Error in transforming linear constraints if SCALE = 1
% 020824  hkh  Accept local solution only if linear constraints fulfilled
% 020828  hkh  Use GetSolver to select default local solver
% 021020  hkh  Check dCon > 0 before using xnargin
% 021020  hkh  Use Fpen=F+sum(max(0,dev)) to handle infeasible initial points
% 030222  hkh  Add check for empty Result.x_k, use best infeasible solution
% 030514  hkh  Avoid size as variable name in randomtest, use p instead
% 030514  hkh  Use rbf_c interface if pEst set, because x is expanded
% 030514  hkh  Also do xnargin on rbf_d2c (normally not used)
% 040104  hkh  Use NARGSAVE to set NARG, must set nlp_c=nlp_xc = [];
% 040104  hkh  sum2 changed to sumv2, sum3 to sumv3
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define fields mLin and mNonLin in DACEProb and EGOProb
% 040309  hkh  Revised functions gutmann and corners, added new sampleInts
% 040309  hkh  Set field Prob.MIP in DACEProb.MIP and EGOProb.MIP
% 040309  hkh  Set globalSolver='glcCluster' if Prob.MIP.IntVars is non-empty
% 040309  hkh  Making EGO handle mixed-integer problems
% 040309  hkh  Add CGO.RandState, to initialize the random generator
% 040309  hkh  Check for cycle of same points
% 040309  hkh  Set Result.Inform to convflag
% 040309  hkh  Revise feasiblity handling, print fMinI, fMinF, use Feasible
% 040309  hkh  Special treatment of pure IP, with glcSolve and glcFast
% 040310  hkh  Add new CGO.pEst parameter, if true estimate || ||_p (default)
% 040323  hkh  Revised with similar initial point features as rbfSolve
% 040323  hkh  Set cTol and bTol in suboptimizations, avoiding infeasibility
% 040323  hkh  Default for Percent if dLin+dCon > 0, use -5000 instead of -d 
% 040404  hkh  Default for pEst=0, for SCALE = 1
% 040404  hkh  Major revision, global solvers do not return empty
% 040406  hkh  Save fMinIdx (=fIdx(1)) on cgoSave.mat
% 040414  hkh  Do not set MENU or LargeScale in EGOProb and DACEProb
% 040414  hkh  Return Result.c_k and Result.Ax
% 040415  hkh  Safe guard if localSolver return > 1 solutions
% 040415  hkh  Avoid LOCAL=1 if pure IP problems
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 040923  hkh  Wrong XL XU when pEst=0, DACEProb wrong.
% 040923  hkh  Changed text to "Cannot proceed", not Failure, if ill-conditioned
% 040929  hkh  sampleInts, x_L(j), offset from 0, must be added
% 041005  hkh  Add localSolver and globalSolver to SolverAlgorithm output
% 041123  hkh  Change call to tomRun


