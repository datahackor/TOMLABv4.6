% TOMLAB NLPJOB NLP Solver
%
% function Result = nlpjobTL(Prob)
%
% INPUT:  
%
% Prob   Problem structure in TOMLAB format
%
% -----------------------------------------------
% Fields used in input structure Prob (call Prob=ProbDef; to define Prob)
% -----------------------------------------------
%
% x_L, x_U  Bounds on variables. 
% b_L, b_U  Bounds on linear constraints. 
% c_L, c_U  Bounds on nonlinear constraints. 
% A         Linear constraint matrix
% PriLevOpt Print level in solver.
%
% -----------------------------------------------
% Fields used in Prob.NLPJOB:
% -----------------------------------------------
% maxfun    The integer variable defines an upper bound for the number
%           of function calls during the line search.
%
% maxit     Maximum number of outer iterations, where one iteration
%           corresponds to one formulation and solution of the quadratic
%           programming subproblem, or, alternatively, one evaluation of
%           gradients.
%
% acc       The user has to specify the desired final accuracy (e.g. 1.0e-7).
%           The termination accuracy should not be smaller
%           than the accuracy by which gradients are computed.
%
% scbou     The real variable allows an automatic scaling of the problem functions. If at the
%           starting point x 0, a function value is greater than SCBOU (e.g. E+3), then
%           the function is divided by the square root. If SCBOU is set to any negative
%           number, then the objective function will be multiplied by the value stored in
%           WA(MMAX+1) and the Jth constraint function by the value stored in WA(J),
%           J=1,...,M.
%
% w         Weight vector of dimension L, to be filled with suitable values when calling
%           NLPJOB depending on the transformation model: MODEL=1,10,12,13,14,15
%           - weights, MODEL=2 - bounds, MODEL=3 - bounds for objective functions,
%           MODEL=4,5 - goal values.
%
% fk        For MODEL=2,6,7,11,12,14,15, FK has to contain the optimal values of the
%           individual scalar subproblems when calling NLPJOB.
%
% model     1-15. See the user's guide for more information.
%
% imin      If necessary (model = 2 or 3), imin defines the index of the objective function
%           to be take into account for the desired scalar transformation.
%
% PrintFile Name of NLPJOB Print file. Amount and type of printing determined 
%           by PriLevOpt.
%
% -----------------------------------------------------------------------
%
% OUTPUT: 
%
% Result   Structure with results (see ResultDef.m):
%
% x_k      Solution vector.
% x_0      Initial solution vector.
%
% f_k      Function value at optimum.
% g_k      Gradient of the objective function.
%
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
%
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
%
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
%
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
%
% ExitFlag Exit status.
% ExitText Exit text from solver.
%
% Inform   NLPJOB information parameter.
%
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% Solver   Name of the solver (NLPJOB).
% SolverAlgorithm  Description of the solver.
%
% NLPJOB.act  The logical array indicates constraints, which NLPQL considers to be 
%             active at the last computed iterate.
%
% NLPJOB.u    Contains the multipliers with respect to the actual iterate 
%             stored in x_k. See user's guide.
%
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Revision: 1.8 $
% Written Jun 1, 2003.    Last modified Dec 2, 2004.
%

function Result = nlpjobTL(Prob)

%#function nlJac

if nargin < 1, error('nlpjobTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

if isempty(Prob.LS)
   fprintf('Prob.LS not defined\n')
   error('No least squares or approximation problem')
end

Prob.solvType = 6;          % Constrained LS solver, cls
Prob = iniSolve(Prob,3,1,1);% Init globals

Prob = ksDef(Prob);         % Klaus Schittkowski problem definitions

if ~isfield(Prob,'NLPJOB')  % Make sure NLPJOB field exists
   Prob.NLPJOB = [];
end

Result=ResultDef(Prob);     % Define the Result struct
Result.Solver='NLPJOB';     % Solver name
Result.SolverAlgorithm='Dense SQP';

PriLev=Prob.PriLevOpt;      % Printing level in solver

%
% Define lower and upper bound arrays for NPSOL
%
% Inf are changed to BIG (=1E12), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints
%

BIG=1E12;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

m = m1+m2;

% Initial checks on the inputs
xl = bl(1:n); % Prob.x_L
xu = bu(1:n); % Prob.x_U

me     = Prob.mEQ;
mi     = Prob.mIN;

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( xl,min(xu,x_0(:) ) ); 

[f, g ] = nlresid(x_0, Prob);
L = length(f);  % Number of objective functions (or residuals)

% Define default solver options.
maxfun    = DefPar(Prob.NLPJOB,'maxfun',400);
maxit     = DefPar(Prob.NLPJOB,'maxit',2000);
acc       = DefPar(Prob.NLPJOB,'acc',1e-6);
scbou     = DefPar(Prob.NLPJOB,'scbou',1e14);
PrintFile = DefPar(Prob.NLPJOB,'PrintFile','nlpjob.txt');
options   = [maxfun, maxit, acc, scbou, eps];

% Default objective functions
model     = DefPar(Prob.NLPJOB,'model',1);
imin      = DefPar(Prob.NLPJOB,'imin',1);
w         = DefPar(Prob.NLPJOB,'w',ones(L,1));
fk        = DefPar(Prob.NLPJOB,'fk',zeros(L,1));

[x,Inform,fw,f_k,df_k,g_k,dg_k,u,act] =  ...
    nlpjob(Prob, me, mi, n, L, x_0, bl(1:n), bu(1:n), w, fk, ...
    imin, PriLev, PrintFile, options, model);

Result.f_k  = f_k;  %optimal value
Result.r_k  = fw;   %residuals
Result.J_k  = df_k;

c_k   = zeros(m2,1);
dc_k  = zeros(m2,n);
cixEQ = Prob.cixEQ;
if ~isempty(cixEQ)
   c_k(cixEQ)    = g_k(cixEQ) + Prob.c_L(cixEQ);
   dc_k(cixEQ,:) = dg_k(cixEQ,:);
end
cixLow = Prob.cixLow;
if ~isempty(cixLow)
   c_k(cixLow)    = g_k(cixLow) + Prob.c_L(cixLow);
   dc_k(cixLow,:) = dg_k(cixLow,:);
end
cixUpp = Prob.cixUpp;
if ~isempty(cixUpp)
   c_k(cixUpp)    = Prob.c_U(cixUpp) - g_k(cixUpp);
   dc_k(cixUpp,:) = -dg_k(cixUpp,:);
end

if m1>0, Result.Ax  = Prob.A*x; else Result.Ax = [];  end

Result.c_k  = c_k;
Result.cJac = dc_k;
Result.x_k  = x;

optParam = Prob.optParam;

Result = StateDef(Result, x(1:n), Result.Ax, Result.c_k, ...
   optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 1);

Result.NLPJOB.u = u;
Result.NLPJOB.act = act;

ExitFlag = 99;
switch(Inform)
    case -2        
        ExitText = 'Unknown Error: -2';
    case -1
        ExitText = 'Unknown Error: -1';
    case 0
        ExitFlag = 0;
        ExitText = 'The optimality conditions are satisfied.';
    case 1 % iterlimit
        ExitFlag = 1;
        ExitText = 'The algorithm has been stopped after MAXIT iterations.';
    case 2
        ExitText = 'The algorithm computed an uphill search direction.';
    case 3
        ExitText = 'Underflow occurred.';
    case 4
        ExitText = 'The line search could not be terminated successfully.';
    case 5
        ExitText = 'Length of a working array is too short. Increase print level.';
    case 6
        ExitText = 'There are false dimensions.';
    case 7
        ExitText = 'The search direction is close to zero. Still infeasible';
    case {Inform > 10}
        ExitText = 'QP sub problem error.';
    otherwise
        ExitText = 'Unknown Error';
end

Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 030601 hel Created NLPJOB interface.
% 040113 med Setting up interface.
% 041202 hkh Revise calls to defblbu and StateDef, avoid vector reshuffling
% 041202 hkh Check if NLLS or approximation problem
% 041202 hkh SolvType should be 6, not 4
