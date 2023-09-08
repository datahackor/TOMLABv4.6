% TOMLAB SNOPT NLP Solver
%
% function Result = snoptTL(Prob)
%
% INPUT:  
% Prob   Problem structure in TOMLAB format.
%
% -----------------------------------------------
% Fields used in input structure Prob (call Prob=ProbDef; to define Prob)
% -----------------------------------------------
%
% x_L, x_U  Bounds on variables. 
% b_L, b_U  Bounds on linear constraints. 
% c_L, c_U  Bounds on nonlinear constraints. 
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% PriLevOpt Print level.
% WarmStart If true, use warm start, otherwise cold start.
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution and slacks from previous run.
% hs        State for solution and slacks from previous run.
% nS        Number of superbasics from previous run.
% hElastic  Defines which variables are elastic in elastic mode. hElastic(j):
%           0 = variable j is non-elastic and cannot be infeasible.
%           1 = variable j can violate its lower bound.
%           2 = variable j can violate its upper bound.
%           3 = variable j can violate either its lower or upper bound.
% moremem   Add more memory if SNOPT stops with not enough storage message.
% SpecsFile Name of user defined SPECS file, read BEFORE optPar() is used.
% PrintFile Name of SOL Print file. Amount and type of printing determined 
%           by SPECS parameters or optPar parameters.
%           Output is written on file snoptpri.txt, if not given.
% SummFile  Name of SOL Summary File.
%           Output is written on file snoptsum.txt, if not given.
% To make snopt to not open and not write anything to file:
%    Set SpecsFile and PrintFile empty.
%    Set optPar(1) = 0 AND optPar(3) = 0.
% optPar    Elements > -999 takes precedence over corresponding TOMLAB
%           params.
%
% -----------------------------------------------
% How optPar is used for setting SPECS parameters:
% -----------------------------------------------
%
% optPar  Structure with optimization parameters. 
%
% SNOPT keywords in optPar(#) - and the fields used in optParam
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
% Printing
% 1.  MAJOR PRINT LEVEL             0        1         111111
%     -PriLev=0 to PRINT LEVEL 0, 1 to 1, 2 to 11111, 3 to 111111
% 2.  MINOR PRINT LEVEL             0        0         10      0, 1 or 10
% 3.  PRINT FILE                    0        0                 Fortran Unit #
%     -Set to 9 by default.
% 4.  SUMMARY FILE                  0        0                 Fortran Unit #
%     -Set to 8 by default
% 5.  PRINT FREQUENCY               0        100
%     -PriFreq
% 6.  SUMMARY FREQUENCY             0        100
%     -SummFreq
% 7.  SOLUTION YES/NO               0        1         1       1 = YES; 0 = NO
% 8.  SUPPRESS OPTIONS LISTING      0        0         1       1 = True
%     Also called SUPPRESS PARAMETERS
%
% Convergence Tolerances
% 9.  MAJOR FEASIBILITY TOLERANCE   >0       1E-6
%     -cTol     Feasibility tol on nonlinear constraints
% 10. MAJOR OPTIMALITY TOLERANCE    >0       max(1E-6,(10eps_R)^0.5)=1.73E-6
%     -eps_x    Optimality tolerance
% 11. MINOR FEASIBILITY TOLERANCE   >0       1E-6  
%     -bTol     Feasibility tol on linear constraints 
% 12. MINOR OPTIMALITY TOLERANCE    >0       1E-6  
%     -MinorTolX  Minor Optimality tolerance
%
% Derivative checking
% 13. VERIFY LEVEL                  -1       0         3       {-1,0,1,2,3}
% 14. START OBJECTIVE CHECK AT COL  0        1         nnObj         
% 15. STOP OBJECTIVE CHECK AT COL   0        nnObj     nnObj
% 16. START CONSTRAINT CHECK AT COL 0        1         nnJac          
% 17. STOP CONSTRAINT CHECK AT COL  0        nnJac     nnJac
%
% Scaling
% 18. SCALE OPTION                  0        1 or 2    2       2 if LP,1 if NLP
%     Option 1 changed to QP Scaling 0 by SNOPT if no linear constraints
% 19. SCALE TOLERANCE               >0       0.9       <1
% 20. SCALE PRINT                   0        0         1       1 = True
%
% Other Tolerances
% 21. CRASH TOLERANCE               0        0.1       <1
% 22. LINESEARCH TOLERANCE          >0       0.9       <1
%     -sigma    Line search accuracy, LineSearch.sigma
% 23. LU FACTOR TOLERANCE           1        100 or 3.99       100 if LP
% 24. LU UPDATE TOLERANCE           1        10  or 3.99       10  if LP
% 25  LU DENSITY TOLERANCE          >0       0.6                            
% 26. LU SINGULARITY TOLERANCE      >0       3.25E-11          eps^(0.67)
% 27. PIVOT TOLERANCE               >0       3.25E-11          eps^(0.67)
%     -eps_Rank Rank tolerance     
%
% QP subproblems
% 28. CRASH OPTION                  0        0         3       {0,1,2,3}
% 29. ELASTIC WEIGHT                0        100.0      
% 30. ITERATION LIMIT               0        10000             or 20m, if more
%     -MaxIter
%     Maximal sum of minor iterations
% 31. PARTIAL PRICE                 0        10 or 1           10 for LP
%
% SQP method
% 32. MAXIMIZE                      0        0         1       1=maximize
% 33. FEASIBLE POINT                0        0         1       1=feasible pnt
% 34. VIOLATION LIMIT               >0       10
% 35. MAJOR ITERATIONS LIMIT        >0       1000              or m,  if more
%     -MajorIter  Iterations (major its)       
%     Maximal number of major iterations
% 36. MINOR ITERATIONS LIMIT        >0       1000              or 5m, if more
%     -MinorIter  Iterations (minor its)       
%     Maximal number of minor iterations, i.e. in the solution of QP or simplex
% 37. MAJOR STEP LIMIT              >0       2
% 38. HESSIAN FREQUENCY             >0       99999999
% 39. DERIVATIVE LEVEL              0        3         3       {0,1,2,3}
%     -DerLevel - Is set by snoptTL dependent on Prob.ConsDiff, Prob.NumDiff
% 40. DERIVATIVE LINESEARCH         0        1         1       0=NONDERIVATIVE
%     -LineAlg  Type of line search algorithm
%     -  LineAlg=0 is quadratic - for SNOPT gives quadratic, without gradient values
%     -  LineAlg=1 is cubic - for SNOPT gives cubic, always using gradient values
%     Default: 0 if numerical derivatives, otherwise 1
% 41. FUNCTION PRECISION            >0       3.0E-13           eps^0.8=eps_R
%     -fTol     Relative tolerance on function 
% 42. DIFFERENCE INTERVAL           >0       5.48E-8           eps^0.4
%     -DiffInt
% 43. CENTRAL DIFFERENCE INTERVAL   >0       6.70E-5           eps^{0.8/3}
%     -CentralDiff
% 44. FEASIBLE EXIT                 0        0         1       1=feasible exit
% 45. UNBOUNDED STEP SIZE           >0       1E20
% 46. UNBOUNDED OBJECTIVE           >0       1E15
%
% Hessian approximation (Superbasics as #48 as MINOS, Hessian freq as 38)
% 47. HESSIAN FULL MEMORY           0        1         1       =1 if nnL <= 75
% or  HESSIAN LIMITED MEMORY                                   =0 if nnL >  75
% 48. SUPERBASICS LIMIT             >0       min(500,1+nnL)
%     Default: 500, but n if n-size(A,1)-length(c_L) > 450 & n <= 5000
%     If n > 5000: max(500,n-size(A,1)-length(c_L)) 
% 49. HESSIAN UPDATES               >0       20
% 50. HESSIAN FLUSH                 >0       99999999
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
% 63. LU COMPLETE PIVOTING          0        0         3    0=partial, 1=complete,
%     or LU PARTIAL PIVOTING                                2=rook,    3=diagonal     
%     or LU ROOK PIVOTING      
%     or LU DIAGONAL PIVOTING  
%
% OUTPUT: 
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial  solution vector.
% g_k      Gradient of the function.
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from snopt.m (similar to TOMLAB).
% Inform   SNOPT information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
% Solver   Name of the solver (snopt).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field:
% xs       Solution and slack variables.
% hs       State for variables and slacks in xs.
% nS       #   of superbasics.
% nInf     #   of infeasibilities.
% sInf     Sum of infeasibilities.
%
% -----------------------------------------------------------------------
%
% For a problem description, see conAssign.m
%
% -------------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written July 4, 2000.   Last modified Dec 22, 2004.
%

function Result = snoptTL(Prob,rmode)

%#function nlp_cdcS nlp_fg

if nargin < 2
   rmode = 0;
   if nargin < 1 
      error('snoptTL needs the Prob structure as input');
   end
end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

% funfdf Name of routine [f,gradf] = funfdf(x, Prob, mode, nstate)
%        funfdf=nlp_fg, included in TOMLAB.
% funcdc Name of routine [g,gJac]  = funcdc(x, Prob, mode, nstate)
%        funcdc=nlp_cdcS, included in TOMLAB.

Prob.solvType = 3; % NLP (CON) solver


Prob = iniSolve(Prob,3,Prob.NumDiff~=6,Prob.mNonLin>0 & Prob.ConsDiff~=6);
%Prob = iniSolve(Prob,3,1,1);

Result=ResultDef(Prob);
Result.Solver='SNOPT';
Result.SolverAlgorithm='SNOPT 6.2-2 NLP code';

PriLev=Prob.PriLevOpt;

%
% Define lower and upper bound arrays for SNOPT
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints
%
BIG=1E20;

% Request bl/bu in order (x,nonlinear constraints, linear constraints)
[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 0);

% Initial checks on the inputs
ObjAdd  = 0;

if isfield(Prob.SOL,'moremem')
   moremem = Prob.SOL.moremem;
else
   moremem = [];
end
if isempty(moremem), moremem = 0; end

m       = m1+m2;

% Check if any linear part
probType = Prob.probType;
%if ~any(probType==[2 7 8]) & isfield(Prob.QP,'c')
if isfield(Prob.QP,'c')
   if probType == 2
      c  = [];
      m3 = 0;
      % Does not speed up much to avoid linear part in callback
      %c  = Prob.QP.c;
      %m3 = length(c) > 0;
      %if m3
      %   Prob.USER.f='qp_fX';
      %   Prob.USER.g='qp_gX';
      %end
   elseif any(probType==[7 8])
      c  = Prob.QP.c;
      m3 = length(c) > 0;
   else
      c  = [];
      m3 = 0;
   end
else
   c  = [];
   m3 = 0;
end

nb     = n+m1+m2+m3;

% Check on Warm start, then set hs,xs,nS

if Prob.WarmStart
   % Warm start for SOL solver
   Warm = 1;
   xs   = Prob.SOL.xs;
   hs   = Prob.SOL.hs;
   %xs   = Prob.SOL.xs(1:nb);
   %hs   = Prob.SOL.hs(1:nb);
   nS   = Prob.SOL.nS;
   x_0  = xs(1:n);
else
   % Initial values
   x_0  = Prob.x_0(:);
   if length(x_0) < n, x_0=zeros(n,1); end
   % Safe guard x_0 
   x_0  = max(bl(1:n),min(x_0,bu(1:n)));
   Warm = 0;
   nS   = Prob.SOL.nS;
   hs   = Prob.SOL.hs;
   xs   = x_0;
end

hElast = Prob.SOL.hElastic;

%if isempty(c), c=zeros(n,1); end

nnObj = n; % number of nonlinear variables
% Avoid putting c both in objective, and calling lp_f to get c'x
if strcmpi('lp',checkType([],Prob.probType))
   nnObj = 0;
end

nEqual = 0;

if m2 > 0
   % Determine the sparse problem structure

   % number of nonlinear constraints   (nnCon)
   nnCon=m2;
   if isempty(Prob.ConsPattern)
      gJac=ones(m2,n);
      nnJac=n;
   else
      gJac=Prob.ConsPattern;
      [ix,iy]=find(gJac);
      %[i,j]=find(Prob.ConsPattern);

      % Send linear index from multiple subscripts for nonzero pattern
      Prob.ConsIdx = sub2ind(size(Prob.ConsPattern),ix,iy);

      % Number of nonlinear Jacobian variables (nnJac)
      nnJac=size(gJac,2);
   end
   % How many nonlinear equality constraints 
   %nEqual=sum(bl(n+m1+1:n+m)==bu(n+m1+1:n+m));
   nEqual=sum(bl(n+1:n+m2)==bu(n+1:n+m2));
   % Must call function, in case global values are needed for constraints
   % Because snopt calls constraints first, then functions
   Result.f_0= nlp_f(x_0, Prob);
   % Sometimes maybe also a call to the gradient might be needed.
   % Do not make such a call now, it is highly unlikely.
else
   nnCon      = 0;
   nnJac      = 0;
   gJac       = [];
   Result.f_0= nlp_f(x_0, Prob);
   %Result.f_0 = c'*x_0;
end

% Move nonlinear constraints BEFORE linear constraints

%bl(n+1:n+m)=bl([n+m1+1:n+m,n+1:n+m1]);
%bu(n+1:n+m)=bu([n+m1+1:n+m,n+1:n+m1]);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end 
end 

% Set up the constraint matrix A 
if isempty(gJac)
   A=sparse(Prob.A);
else
   if isempty(Prob.A)
      A = sparse(gJac);
   else
      if nnJac < n
         A = sparse([[gJac,zeros(m2,n-nnJac)];Prob.A]);
      else
         A = sparse([gJac;Prob.A]);
      end
   end
end
if m==0
   % Construct one linear constraint from some bounded variable
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);

   %bl=[bl;max(bl(i),-BIG)];
   %bu=[bu;min(bu(i),BIG)];
   bl(n+1)=max(bl(i),-BIG);
   bu(n+1)=min(bu(i),BIG);
   m=1;
   A = sparse(zeros(1,n));
   A(1,i)=1;
end

if m3 > 0
   A = [A;c'];
   iObj=m+1;  % The constraint row with linear obj term c is put last
   % Add the bounds for the objective in bl and bu
   bl=[bl;-BIG];
   bu=[bu; BIG];
else
   iObj=0;
end

[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('snopt',3,...
         nnObj, nnJac, size(A,1), Prob);

if optPar(40) < 0 & (Prob.NumDiff ~= 0 | (Prob.ConsDiff ~= 0 & m2 > 0) | ...
   isempty(Prob.USER.g) | (isempty(Prob.USER.dc) & m2 > 0))
   % Change default to non-Derivative line search if numerical derivatives
   optPar(40) = 0;
end
if optPar(48) < 0 & (n-m1-m2) > 450
   % Increase number of superbasics
   if n > 5000
      optPar(48) = max(500,n-m1-m2+200);
   else
      optPar(48) = max(500,n);
   end
end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

% Lagrange multipliers, known info
pi(1:m+m3)=0;

if Prob.NumDiff == 6
   if Prob.ConsDiff == 6
      optPar(39)=0;
   else
      optPar(39)=2;
   end
elseif Prob.ConsDiff == 6
   optPar(39)=1;
else
   optPar(39)=3;
end

if rmode > 0
[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
     rsnopt( A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
            Warm, hs, xs, pi, nS, ...
            SpecsFile, PrintFile, SummFile, ...
            PriLev, ObjAdd, moremem, Prob.Name,rmode );
else
[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
     snopt( A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
            Warm, hs, xs, pi, nS, ...
            SpecsFile, PrintFile, SummFile, ...
            PriLev, ObjAdd, moremem, Prob.Name );
end

switch Inform
   case 3
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 1
     ExitFlag=4;  % Infeasible
   case {10,22}
     ExitFlag=3;  % Rank problem
   case {5,7,8,20,21,30,32,40,42,43,44}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k=Obj;
Result.g_k=gObj;
Result.c_k=fCon;
Result.x_k=xs(1:n);
Result.x_0=x_0;
% Multipliers for bounds-linear-nonlinear
Result.v_k=[rc(1:n);pi(m2+1:m2+m1);pi(1:m2)];

if ~isempty(gCon)
   [ix,iy] = find(gJac);
   % TOMLAB now has the format one constraint / row, same as SOL solvers
   Result.cJac=sparse(ix,iy,gCon,size(gJac,1),size(gJac,2));
end

% Saved for a warm start
%if m3 > 0
   Result.SOL.xs=xs;
   Result.SOL.hs=hs;
%else
%   % The dummy objective row is last in xs, value -Obj, last in hs is 3 (basic)
%   Result.SOL.xs=[xs;-Obj];
%   Result.SOL.hs=[hs;3];
%end
Result.SOL.nS=nS;
Result.SOL.nInf=nInf;
Result.SOL.sInf=sInf;
Result.SOL.optPar=optPar;

Result.FuncEv    = iwCount(3);
Result.GradEv    = sum(iwCount(4:6));
Result.ConstrEv  = iwCount(7);
Result.Iter      = iwCount(2);
Result.MinorIter = iwCount(1);
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

optParam = Prob.optParam;
if m1 > 0
   Result.Ax = A(m2+1:m2+m1,:)*xs(1:n);
   Result    = StateDef(Result, xs(1:n), Result.Ax, fCon, ...
                        optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 0);
else
   Result    = StateDef(Result, xs(1:n), [], fCon, ...
                        optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 0);
end

switch Inform
   case 0
      Text = 'Optimal solution found';
   case 1
      Text = 'Problem Infeasible';
   case 2
      Text = 'Problem Unbounded (or badly scaled)';
   case 3
      Text = 'Too many iterations';
   case 4
      Text = str2mat('Feasible solution, but the requested accuracy'...
             ,'in the dual feasibilities could not be achieved');
   case 5
      Text = 'The Superbasics limit is too small';
   case 6
      Text = str2mat('User requested termination by returning' ...
             ,'mode <= -2 from funobj or funcon');
   case 7
      Text = 'funobj seems to give incorrect gradients';
   case 8
      Text = 'funcon seems to give incorrect gradients';
   case 9
      Text = 'The current point cannot be improved';
   case 10
      Text = str2mat('Numerical error in trying to satisfy the linear ' ...
             ,'constraints. The basis is very ill-conditioned.');
   case 12
      Text = str2mat('Basis factorization requested twice in a row' ...
             ,'Like case Inform = 9. Possibly convergence?');
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
   case 42
      Text = 'Not enough 8-character workspace to solve the problem';
   case 43
      Text = 'Not enough integer workspace to solve the problem';
   case 44
      Text = 'Not enough real workspace to solve the problem';
   otherwise
      Text = 'NOTE: UNKNOWN SNOPT Inform value.';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nSNOPT solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('SNOPT: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Major       iterations%7d. ',iwCount(2));
   fprintf('Total minor iterations%7d. ',iwCount(1));
   fprintf('\n');
   
   fprintf('fObj and gObj evaluations%7d %d %d %d\n',iwCount(3:6));
   if m2 > 0
      fprintf('fCon and gCon evaluations%7d %d %d %d\n',iwCount(7:10));
   end
   fprintf('nInf (# of infeasible constraints) %7d. ',nInf);
   fprintf('nS (# of superbasics) %7d. ',nS);
   fprintf('sInf (Sum of infeasibilities outside bounds) %14.7e\n',sInf);

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
% 000704 hkh  New snoptMex interface.
% 000916 hkh  Return string ExitText with interpretation of Inform flag 
% 010901 hkh  Give snopt correct number of nonlinear Jacobian variables 
%             using size of Prob.ConsPattern
% 010903 hkh  Add complete pivoting option. Allow cold start change of hs,ns
% 011212 hkh  Improve comments on iterations
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 020210 hkh  Send linear index for sparse gJac subscripts in Prob.ConsIdx
% 020304 hkh  Set optPar(39) DERLVL dependent on Prob.ConsDiff, Prob.NumDiff
% 020417 hkh  Send optPar back in Result.SOL
% 021105 hkh  Make a comment about SCALE OPTION 1 set to 0 when only linear con
% 030801 hkh  Correct test for LP problem, avoid f==0 results
% 031208 ango Edit comments (optPar(63)), Result.SolverAlgorithm for SNOPT 6.2-2
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040304 hkh  Safe guard x_0 x_0 = max(bl(1:n),min(bu(1:n),x_0));
% 040602 med  Help fix PriLev to PriLevOpt
% 040728 med  Pragmas added for MATLAB Compiler
% 041005 ango Change order in Result.v_k: bounds, linear, nonlinear
% 041129 hkh  Default non-derivative line search if numerical derivatives
% 041129 hkh  Different Number of Superbasics dependent on problem size
% 041202 hkh  Removed unnecessary print statements from old debug
% 041202 hkh  Revise calls to defblbu and StateDef, use different Order=0
% 041202 hkh  If internal differentiation, different call to iniSolve 
% 041222 med  Safeguard added for x_0
