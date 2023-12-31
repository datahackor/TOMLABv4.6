% TOMLAB NLPQL Solver
%
% function Result = nlpqlTL(Prob)
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
% A         Linear constraint matrix.
% PriLevOpt Print level in solver.
% WarmStart If true, use warm start, otherwise cold start.
%
% -----------------------------------------------
% Fields used in Prob.NLPQL:
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
% accqp     The tolerance is needed for the QP solver to perform several 
%           tests, for example whether optimality conditions are satisfied
%           or whether a number is considered as zero or not. If ACCQP is
%           less or equal to zero, then the machine precision is computed 
%           by NLPQL and subsequently multiplied by 1.0e+4.
%
% PrintFile Name of NLPQL Print file. Amount and type of printing determined
%           by PriLevOpt.
%
% The following parameters are required if WarmStart is used.
%
% u         Contains the multipliers with respect to the actual iterate
%           stored in the first column of X. The first M locations contain
%           the multipliers of the M nonlinear constraints, the subsequent
%           N locations the multipliers of the lower bounds, and the
%           final N locations the multipliers of the upper bounds.
%           At an optimal solution, all multipliers with respect to
%           inequality constraints should be nonnegative.
%
% c         On return, C contains the last computed approximation
%           of the Hessian matrix of the Lagrangian function stored in
%           form of an LDL decomposition. C contains the lower triangular
%           factor of an LDL factorization of the final quasi-Newton matrix
%           (without diagonal elements, which are always one).
%           In the driving program, the row dimension of C has to be equal
%           to NMAX.
%
% d         The elements of the diagonal matrix of the LDL decomposition
%           of the quasi-Newton matrix are stored in the one-dimensional
%           array D.
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
% Inform   NLPQL information parameter.
%
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% Solver   Name of the solver (NLPQL).
% SolverAlgorithm  Description of the solver.
%
% NLPQL.act  The logical array indicates constraints, which NLPQL considers to be 
%            active at the last computed iterate.
%
% NLPQL.u    See inputs.
% NLPQL.c    See inputs.
% NLPQL.d    See inputs.
%
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0 $
% Written Jun 1, 2003.    Last modified Dec 2, 2004.
%

function Result = nlpqlTL(Prob)

if nargin < 1, error('nlpqlTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

Prob.solvType = 3;          % NLP (CON) solver
Prob = iniSolve(Prob,3,1,1);% Init globals

Prob = ksDef(Prob);         % Klaus Schittkowski problem definitions

if ~isfield(Prob,'NLPQL')   % Make sure NLPQL field exists
   Prob.NLPQL = [];
end

Result=ResultDef(Prob);     % Define the Result struct
Result.Solver='NLPQL';      % Solver name
Result.SolverAlgorithm='Dense SQP';

PriLev=Prob.PriLevOpt;      % Printing level in solver

%
% Define lower and upper bound arrays for NLPQL
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

% Initial checks on the inputs

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

m  = m1+m2;

xl = bl(1:n); % Prob.x_L
xu = bu(1:n); % Prob.x_U

% Check if Warm Start, then set U, C and D

u = DefPar(Prob.NLPQL,'u',[]);
c = DefPar(Prob.NLPQL,'c',[]);
d = DefPar(Prob.NLPQL,'d',[]);

WarmStart = DefPar(Prob,'WarmStart',0);
if WarmStart
    % Warm start for NLPQL solver
    u    = Prob.NLPQL.u;
    c    = Prob.NLPQL.c;
    d    = Prob.NLPQL.d;
end

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( xl,min(xu,x_0(:) ) ); 

Result.x_0 = x_0;

mTot = Prob.m;
me  = Prob.mEQ;
n   = size(x_0,1);
xl  = Prob.x_L;
xu  = Prob.x_U;

[f, g ] = nlfunc(x_0, Prob);
[df,dg] = nlgrad(x_0, Prob);
Result.f_0 = f;

% Define default solver options.
maxfun  = DefPar(Prob.NLPQL,'maxfun',400);
maxit   = DefPar(Prob.NLPQL,'maxit',2000);
acc     = DefPar(Prob.NLPQL,'acc',1e-6);
accqp   = DefPar(Prob.NLPQL,'accqp',1e-6);
PrintFile = DefPar(Prob.NLPQL,'PrintFile','nlpql.txt');

options = [maxfun, maxit, acc, accqp];

[x,Inform,f_k,g_k,df_k,dg_k,u,c,d,act] =  ...
    nlpql(Prob, mTot, n, me, f, g, df, dg, x_0, xl, xu, PriLev, ...
	  PrintFile, options, WarmStart, u, c, d);

Result.f_k  = f_k;
Result.g_k  = df_k;

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
   optParam.xTol, optParam.bTol, optParam.cTol, bl, bu,1);

Result.NLPQL.act = act;

% Outputs for Warm Start

Result.NLPQL.u = u;
Result.NLPQL.c = c;
Result.NLPQL.d = d;

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
    case 8
        ExitText = 'The starting point violates a lower or upper bound.';
    case 9
        ExitFlag = 9;
        ExitText = 'Wrong input parameter.';
    case {Inform > 9}
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
% 030601 hel Created NLPQL interface.
% 040103 med Setting up interface.
% 040109 hkh Clean up. Remove ProbCheck call.
% 040117 med New mTot for total number of constraints
% 040119 med New default settings.
% 041202 hkh Revise calls to defblbu and StateDef, avoid vector reshuffling
