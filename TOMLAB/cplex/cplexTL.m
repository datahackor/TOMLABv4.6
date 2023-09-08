% TOMLAB /CPLEX LP, QP, MILP, MIQP and MIQQ Solver
%
% cplexTL converts the problem from the TOMLAB structure format and 
% calls cplex.m. On return converts the result to the TOMLAB structure format.
% Also see the help for cplex.m
% 
% cplexTL.m solves the following mixed integer (linear or quadratic) 
% programming problem (LP, QP, MILP, MIQP, MIQQ):
%
%   minimize   0.5 * x'*F*x + c'x    
%      x            
% 
%   subject to  x_L <=    x   <= x_U
%               b_L <=   Ax   <= b_U
%
%   and also, optionally, subject to the quadratic constraints
%
%       x'*Q_i*x + a_i'*x <= r_U(i), i = 1,2,...,n_qc 
%
%   where 
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the CPLEX-MP sparse matrix format.
%   c, x_L, x_U has dimension n
%   b_L, b_U has dimension m
%   F is a n x n symmetric matrix, sparse or dense. 
%   If F is empty, an LP or MILP problem is solved
%   Q_i are sparse n x n Matlab matrices, a_i has dimension n and r_U(i) is a
%   scalar upper bound for the i:th quadratic constraint
%
%   Some or all x may be integer valued as specified by other input
%   variables.
%
% ---------------------------------------------------------------------------
%
% function Result = cplexTL(Prob)
%
% INPUT:    
%
% Prob        Problem structure in TOMLAB format.
%             Use lpAssign, mipAssign, qpAssign, miqpAssign to 
%             define the Prob structure.
%
%             Fields used in input structure Prob: 
%
%
%   x_L, x_U  Lower and upper bounds on variables, size n x 1
%   b_L, b_U  Lower and upper bounds on linear constraints, size m x 1 
%   A         Linear constraint matrix, dense or sparse m x n matrix
%   
%             NOTE - all bounds vectors - if [], +/- Inf is assumed
%   
%   QP.c      Linear objective function coefficients, size n x 1
%   QP.F      Quadratic matrix of size n x n 
%
%   QP.qc     Structure array defining quadratic constraints ("qc").
%
%             Please note that CPLEX 9.0 only handles single-sided bounds 
%             on qc's. An arbitrary number of qc's is set using the Prob.QP.qc
%             structure array:
%
%          qc(1).Q   = sparse( <quadratic coefficient nxn matrix> );
%          qc(1).a   = full  ( <linear coefficient nx1 vector   > );
%          qc(1).r_U = <scalar upper bound> ;  
%
%             And similarly for qc(2), ... , qc(n_qc).
%
%             The standard interpretation is x'*Q*x + c'*x <= r_U, but it is 
%             possible to define an alternative sense x'*Q*x + c'*x >= r_L 
%             by setting qc(i).sense to a nonzero value and specifying a
%             lower bound in qc(i).r_L.
%
%             Observe that the Q matrix must be sparse, non-empty and positive 
%             semi-definite for all qc's. The linear coefficient vector qc(i).a 
%             may be omitted or set empty, in which case all zeros are assumed. 
%
%             Likewise, if a bound r_U or r_L is empty or not present, it
%             is assumed to be 0.0. Note that this is contrary to the usual
%             Tomlab standard, where an empty or omitted bound is assumed
%             to be +/- Inf. The reason is that a single-sided constraint with 
%             an infinite bound would have no meaning. 
%
%   PriLevOpt Print level in cplexTL, the cplex m-file and cplexmex C-interface.
%             = 0  Silent
%             = 1  Warnings and Errors
%             = 2  Summary information
%             = 3  More detailed information
%             > 10 Pause statements, and maximal printing (debug mode)
%
%
% Fields used in Prob.CPLEX (Structure with CPLEX specific parameters)
%
%   LogFile   Name of file to receive the CPLEX iteration and results log. 
%             If empty or not present, no log is written. 
%
%   SaveFile  Name of file for saving the CPLEX problem just prior to calling the 
%             CPLEX solver. If empty, nothing is written. Also see the SaveMode 
%             parameter below.
%
%   SaveMode  Integer flag indicating which format to use for the save file.
%             The following values are possible:
%
%          1: SAV  Binary SAV file                 (default) 
%          2: MPS  MPS format, original format
%          3: LP   CPLEX LP format, original format
%          4: RMP  MPS file with generic names
%          5: REW  MPS file with generic names
%          6: RLP  LP  file with generic names
%
%             The SAV format is a binary format suitable for submission to
%             ILOG help desk.
%
%   iis       Flag indicating whether to compute an IIS and return it to
%             MATLAB. This option can only be set for an LP problem using
%             the algorithms primal or dual simplex or barrier with
%             crossover.
%               
%             = 0     Don't return IIS to MATLAB (default).
%             = 1     Compute IIS and return it to MATLAB if an LP problem 
%                     has been proven infeasible.
%
%             The IIS is returned through the output parameter 'iis'.
%
%   iisFile   Name of a file to write the IIS to. No file is written if
%             this input parameter is empty or if no IIS is available.
%
% sa          Structure telling whether and how you want CPLEX to perform a
%             sensitivity analysis (SA). You can complete an SA on the 
%             objective function, right hand side vector, lower and 
%             upper bounds. The sa structure contains four
%             sub structures:
%
%                 .obj, .rhs, .xl, .xu
%               
%             Each one of these contain the field:
%
%                 .index
%
%             .index contain indices to variables or constraints 
%             of which to return possible value ranges.
%
%             The .index array has to be sorted, ascending.
%
%             To get an SA of objective function on the four variables 120 
%             to 123 (included) and variable 19, the sa structure 
%             would look like this:
%
%                 sa.obj.index = [19 120 121 122 123];
%
%             The result is returned through the output parameter 'sa'.
%
%
% Fields used in Prob.MIP:
% 
% See the corresponding variables in cplex.m for an explanation
%
%   MIP.IntVars
%             Defines which variables are integers, of general type I or binary B
%             Variable indices should be in the range [1,...,n].
%             IntVars is a single integer ==> Variable 1:IntVars are integer 
%             IntVars is a logical vector ==> x(find(IntVars > 0)) are integers 
%             IntVars is a vector of indices ==> x(IntVars) are integers 
%             (if [], then no integers of type I or B are defined)
%             cplex checks which variables has x_L=0 and x_U=1, i.e. binary.
%   
%   MIP.PI
%             Integer variables of type Partially Integer (PI), i.e. takes an
%             integer value up to a specified limit, and any value above that
%             limit.
%             PI must be a structure array where:
%             PI.var  Vector of variable indices in the range [1,...,n]
%             PI.lim  A vector of limit values, for each of the variables
%                     specified in PI.var, i.e. for variable i, 
%                     that is the PI variable with index j in PI.var:
%                     x(i) takes integer values in [x_L(i),PI.lim(j)] and
%                     continous values in [PI.lim(j),x_U(i)].
%   
%   MIP.SC    A vector with indices for the integer variables of type 
%             Semi-continuous (SC), i.e. that takes either the value 0 or a 
%             real value in the range [x_L(i),x_U(i)], assuming for some j,
%             i = SC(j), where i is an variable number in the range [1,n].
%   
%   MIP.SI    A vector with indices for the integer variables of type 
%             Semi-integer (SI), i.e. that takes either the value 0 or 
%             an integer value in the range [x_L(i),x_U(i)], assuming for some j,
%             i = SI(j), where i is an variable number in the range [1,n].
%   
%   MIP.sos1  A structure defining the Special Ordered Sets of Type One (sos1). 
%             Assume there are k sets of type sos1, then
%             sos1(1).var is a vector of indices for variables in sos1, set 1.
%             sos1(1).row is the row number for the reference row identifying
%                         the ordering information for the sos1 set, i.e.
%                         A(sos1(1).row,sos1(1).var) identifies this information
%             sos1(2).var is a vector of indices for variables in sos1, set 2.
%             sos1(2).row is the row number for the reference row of sos1 set 2.
%             ...
%             sos1(k).var is a vector of indices for variables in sos1, setk.
%             sos1(2).row is the row number for the reference row of sos1 set k.
%   
%   MIP.sos2  A structure defining the Special Ordered Sets of Type Two (sos2). 
%             Specified exactly as sos1 sets, see sos1 input variable description
%   
%   MIP.basis Vector containing a CPLEX Basis. If re-solving a similar problem 
%             several times, this can be set to the 'basis' output argument of an 
%             earlier call to cplex.m. 
%
%             The length of this vector must be equal to the sum  of the number 
%             of rows (m) and columns (m). 
%
%             Furthermore, please note that if cpxControl.ADVIND is set to zero, 
%             the advanced basis information will not be used.
%
%             The first m elements contain row basis information, with the
%             following possible values for non-ranged rows:
%
%           0 associated slack/surplus/artificial variable nonbasic at value 0.0 
%           1 associated slack/surplus/artificial variable basic 
%
%             and for ranged rows (both upper and lower bounded)
% 
%           0 associated slack/surplus/artificial variable nonbasic at its lower bound  
%           1 associated slack/surplus/artificial variable basic 
%           2 associated slack/surplus/artificial variable nonbasic at its upper bound  
%
%
%            The last n elements, i.e. basis(m+1:m+n) contain column 
%            basis information:
%
%           0 variable at lower bound 
%           1 variable is basic 
%           2 variable at upper bound 
%           3 variable free and nonbasic 
%
%   MIP.cpxControl
%
%           cpxControl Structure, where the fields are set to the CPLEX Parameters of
%           type CPX_PARAM_, see the CPLEX manuals, or the Tomlab /CPLEX 
%           User's Guide.
%           The user only sets the fields corresponding to the parameters
%           to be changed from its default values.
%
%           NOTE - the prefix CPX_PARAM_ is NOT used in the field name
%
%           Examples: 
% 
%             cpxControl.STARTALG = 1 Solve root node with Primal instead of Dual simplex
%             cpxControl.SUBALG   = 4 Solve subnodes with Barrier with crossover
%
% Fields used in Prob.optParam: (Structure with optimization parameters)
% 
%   MaxIter   Limit of iterations  (if not cpxControl.ITLIM is set)
%
%   MIP.callback
%
%           Logical vector defining which callbacks to use in CPLEX
%           If the i:th entry of logical vector callback is set, the corresponding 
%           callback is defined. See Tomlab /CPLEX User's Guide
%           The callback calls the m-file specified in cplex.m. 
%           The user may edit the existing files, or make copies and put
%           them before the predefined files in the Matlab path.
%
% ------------------------------------------------------------------------------
%
% OUTPUTS: 
%
% Result   Structure with results (see ResultDef.m):
%
%   f_k      Function value at optimum
%   x_k      Solution vector
%   x_0      Initial  solution vector not known, set as empty
%   g_k      Exact gradient computed at optimum, computed as c or c + Fx
%
%   xState   State of variables.   Free==0; On lower == 1; On upper == 2; Fixed == 3;
%   bState   State of constraints. Free==0; On lower == 1; On upper == 2; Equality == 3;
%
%   v_k      Lagrangian multipliers (for bounds + dual solution vector)
%            v_k = [rc;v]. rc n-vector of reduced costs. v holds m dual variables
%
%   rc       Reduced costs. If ninf=0, last m == -v_k
%
%   ExitFlag Exit status, TOMLAB standard
%
%   Inform   CPLEX information parameter, see Tomlab /CPLEX User's Guide 
%   
%      LP/QP Inform values
%
%       1 (S,B) Optimal solution found
%       2 (S,B) Model has an unbounded ray
%       3 (S,B) Model has been proved infeasible
%       4 (S,B) Model has been proved either infeasible or unbounded
%       5 (S,B) Optimal solution is available, but with infeasibilities after unscaling
%       6 (S,B) Solution is available, but not proved optimal, due to numeric difficulties 
%      10 (S,B) Stopped due to limit on number of iterations
%      11 (S,B) Stopped due to a time limit
%      12 (S,B) Stopped due to an objective limit
%      13 (S,B) Stopped due to a request from the user
%      
%      20 (B)   Model has an unbounded optimal face
%      21 (B)   Stopped due to a limit on the primal objective
%      22 (B)   Stopped due to a limit on the dual objective
%
%      32 Converged, dual feasible, primal infeasible';
%      33 Converged, primal feasible, dual infeasible';
%      34 Converged, primal and dual infeasible';
%      35 Primal objective limit reached';
%      36 Dual objective limit reached';
%      37 Primal has unbounded optimal face';
%      38 Non-optimal solution found, primal-dual feasible';
%      39 Non-optimal solution found, primal infeasible';
%      40 Non-optimal solution found, dual infeasible';
%      41 Non-optimal solution found, primal-dual infeasible';
%      42 Non-optimal solution found, numerical difficulties';
%      43 Barrier found inconsistent constraints';
%
%     101 Optimal integer solution found
%     102 Optimal sol. within epgap or epagap tolerance found
%     103 Solution is integer infeasible
%     104 The limit on mixed integer solutions has been reached 
%     105 Node limit exceeded, integer solution exists
%     106 Node limit exceeded, no integer solution
%     107 Time limit exceeded, integer solution exists
%     108 Time limit exceeded, no integer solution
%     109 Terminated because of an error, but integer solution exists.
%     110 Terminated because of an error, no integer solution 
%     111 Limit on tree memory has been reached, but an integer solution exists 
%     112 Limit on tree memory has been reached; no integer solution 
%     113 Stopped, but an integer solution exists.
%     114 Stopped; no integer solution.
%     115 Problem is optimal with unscaled infeasibilities 
%     116 Out of memory, no tree available, integer solution exists 
%     117 Out of memory, no tree available, no integer solution 
%     118 Model has an unbounded ray 
%     119 Model has been proved either infeasible or unbounded 
%
%   Iter     Number of iterations / nodes visited
%
%   FuncEv   Number of function evaluations. Set to Iter.
%
%   GradEv   Number of gradient evaluations. Set to Iter if
%            QP/MIQP, otherwise 0. 
%            FuncEv and ConstrEv set to Iter. GradEv=0.
%
%   ConstrEv Number of constraint evaluations. Set to 0.
%
%   QP.B     Basis vector in TOMLAB QP standard ???
%
%   Solver           Name of the solver  (CPLEX)
%   SolverAlgorithm  Description of the solver
%
% -----------------------------------------------
% Output fields in Result.MIP:
% -----------------------------------------------
%   MIP.slack     Slack variables (m x 1 vector)
%   MIP.ninf      Number of infeasibilities
%   MIP.sinf      Sum of infeasibilities
%   MIP.lpiter    Number of LP iterations
%   MIP.glnodes   Number of nodes visited
%   MIP.basis     basis status of constraints + variables, (m + n x 1 vector)
%                 in the CPLEX format, fields xState and bState has the same
%                 information in the Tomlab format.
%
% Also set into the Result.MIP output structure is:
%
%   MIP.cpxRetVec  A vector with information on return from CPLEX, see the 
%                  Tomlab /CPLEX User's Guide for a description of each element
%
% -----------------------------------------------
% Output fields in Result.CPLEX:
% -----------------------------------------------
%
% iis       Structure of information about an IIS if requested. The fields:
%
%               iisStatus   Status information. Positive on success,
%                           negative on failure. Possible negative
%                           values and their meaning:
%                        
%                           -1  Problem is not infeasible.
%                           -2  The algorithm used was barrier with no
%                               crossover. 
%                           -3  The presolve for the barrier algorithm
%                               found the infeasibility. See output
%                               (if PriLev >= 1) to get information 
%                               about the infeasibility. (It is not
%                               possible to turn the barrier presolve
%                               off.)
%                           -4  Internal CPLEX error. See iisMessage.
%                           -5  Problem was infeasible or unbounded.
%                               See the presolver output (PriLev >= 1) 
%                               for a description of the infeasibility
%                               or the unbounded ray. If the presolver
%                               was enabled, disable it and try again
%                               to obtain an IIS if the problem was
%                               infeasible and not unbounded.
%
%                           If iisFile was defined, iisStatus should be
%                           equal to or greater than 2 on return. If not,
%                           the file output failed.
%
%               iisMessage  An error message if an error occurred.
%
%               iisstat     Status information from CPLEX. Possible values:
%                       
%                            1  IIS found.
%                            2  Time limit exceeded, partial IIS found.
%
%               rowind      An array containing the indices of the rows in 
%                           the IIS. 
%
%               rowbdstat   An array that identifies the right-hand side
%                           value for each row in the IIS that is causing
%                           the infeasibility. Possible values:
%
%                            0  Row is at lower bound.
%                            1  Row is fixed.
%                            2  Row is at upper bound.
%
%               colind      An array containing the indices of the columns in 
%                           the IIS. 
%
%               colbdstat   An array that identifies the bound for each column 
%                           in the IIS that is causing the infeasibility. 
%                           Possible values:
%
%                            0  Column is at lower bound.
%                            1  Column is fixed.
%                            2  Column is at upper bound.
%
%   sa      Structure with information about the requested SA, if requested.
%           The fields:
%
%               obj         Ranges for the variables in the objective function.
%
%               rhs         Ranges for the right hand side values.
%   
%               xl          Ranges for the lower bound values.
%
%               xu          Ranges for the upper bound values.
%
%           These fields are structures themselves. All four structures 
%           have identical field names:
%
%               status      Status of the SA operation. Possible values:
%
%                            1  Successful
%                            0  SA not requested.
%                           -1  Error: begin is greater than end.
%                           -2  Error: The selected range (begin...end) stretches 
%                               out of available variables or constraints.
%                           -3  Error: No SA available.
%
%               lower       The lower range.
%
%               upper       The upper range.
%
%
% The return vector information is also available as a global variable 
% after the run, do: global cpxRetVec
%
% -----------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 9.0.2 $
% Written July 8, 1999.      Last modified Feb 1, 2005.
%

function Result = cplexTL(Prob)

%#function lp_f lp_g lp_H cpx2cbinfo cpx2retvec

if nargin < 1 
    error('cplexTL needs the Prob structure as input');
end

% Information on return from CPLEX
global cpxCBInfo cpxRetVec

global MAX_x MAX_c % Max number of variables and constraints to print

Prob.solvType = 11; % MIQP solver

Prob = iniSolve(Prob,11,0,0);

nCALLBACKS=11;
%DEBUG = 0;

PriLev=Prob.PriLevOpt;

Result=ResultDef(Prob);
Result.Solver='CPLEX';
Result.SolverAlgorithm='CPLEX LP/QP/MILP/MIQP/MIQQ solver';

% Initial checks on the inputs
%
% Define lower and upper bound arrays for CPLEX
%
% Inf are changed to BIG (default = 1E12), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%

% The cplex MEX has 1E12 as default, but is changed if BIG is noempty
Prob.BIG=DefPar(Prob,'BIG',1E12);

[bl, bu, n, m, m2] = defblbu(Prob, Prob.BIG);

nTot=n+m;

% Initial values (cplex does not use x_0)

Fzero = isempty(Prob.QP.F);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   %   fprintf('mA = %i, m = %i\n', mA, m);
   if mA~=m, error('Linear constraints A MUST have m rows!'); end 
end 

Result.f_0=0;

% Check if any linear part
c = Prob.QP.c(:);

if isempty(c), c=zeros(n,1); end


Prob.MIP.IntVars=DefPar(Prob.MIP,'IntVars',[]);
Prob.MIP.sos1=DefPar(Prob.MIP,'sos1',[]);
Prob.MIP.sos2=DefPar(Prob.MIP,'sos2',[]);
Prob.MIP.SC=DefPar(Prob.MIP,'SC',[]);
Prob.MIP.SI=DefPar(Prob.MIP,'SI',[]);
Prob.MIP.PI=DefPar(Prob.MIP,'PI',[]);

if isempty(Prob.MIP)
   MIP=0;
   
elseif isempty(Prob.MIP.IntVars)
   if isempty(Prob.MIP.sos1) & isempty(Prob.MIP.sos2) & ...
         isempty(Prob.MIP.SC)  & isempty(Prob.MIP.SI)  & isempty(Prob.MIP.PI)
      MIP=0;
   else
      MIP=1;
   end
else
   MIP=1;
end

cpxControl = DefPar(Prob.MIP,'cpxControl',[]);
callback   = DefPar(Prob.MIP,'callback',zeros(nCALLBACKS,1));

if length(callback) < nCALLBACKS
   callback=[callback(:);zeros(nCALLBACKS-length(callback),1)];
end

if MIP
   IntVars = DefPar(Prob.MIP,'IntVars',[]);
   PI      = DefPar(Prob.MIP,'PI',[]);
   SC      = DefPar(Prob.MIP,'SC',[]);
   SI      = DefPar(Prob.MIP,'SI',[]);
   sos1    = DefPar(Prob.MIP,'sos1',[]);
   sos2    = DefPar(Prob.MIP,'sos2',[]);
   
   %if DEBUG
   %   callback(1)=0;
   %   callback(2)=0;
   %   callback(3)=0;
   %   callback(4)=0;
   %   callback(5)=0;
   %   callback(6)=0;
   %   callback(7)=0;
   %   callback(8)=0;
   %   %callback(9)=1;
   %   %callback(10)=1;
   %   %callback(11)=1;
   %   callback(9)=0;
   %   callback(10)=0;
   %   callback(11)=0;
   %end
   
else
   IntVars = [];
   PI      = [];
   SC      = [];
   SI      = [];
   sos1    = [];
   sos2    = [];
   
   %if DEBUG
   %   %      callback(8)=1;
   %   %      callback(9)=0;
   %   %      callback(1:7)=0;
   %   %      callback(10)=1;
   %   %      callback(11)=1;
   %end
end

%Simplex iteration limit
if ~isfield(Prob.optParam,'MaxIter')
   Prob.optParam.MaxIter=2000;
end
if ~isfield(cpxControl,'ITLIM') & Prob.optParam.MaxIter~=2000
   cpxControl.ITLIM=Prob.optParam.MaxIter;
end

CPLEX = DefPar(Prob,'CPLEX');

SaveFile = DefPar(CPLEX,'SaveFile','');
SaveMode = DefPar(CPLEX,'SaveMode',1);
LogFile  = DefPar(CPLEX,'LogFile','');

if ~isempty(SaveFile) & (SaveMode<1 | SaveMode>6)
   error('SaveFile is nonempty and SaveMode is not in 1,2,...,6')
end

% Advanced starting basis, if available
advbas = DefPar(Prob.MIP,'basis',[]);

% Switch ADVIND on, if not explicitly set by the user
if length(advbas) == nTot
   ADVIND = DefPar(cpxControl,'ADVIND',[]);
   if isempty(ADVIND)
      cpxControl.ADVIND = 1;
   end
end

iisRequest = DefPar(CPLEX,'iis',[]);
iisFile    = DefPar(CPLEX,'iisFile','');
saRequest  = DefPar(CPLEX,'sa',[]);

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('cplex',2,...
%         nnObj, 0, size(Prob.A,1), Prob);

% Quadratic constraints
qc = DefPar(Prob.QP,'qc',[]);


[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, ...
        iis, sa] = ...
   cplex(c, Prob.A, bl(1:n), bu(1:n), bl(n+1:nTot), bu(n+1:nTot), cpxControl,...
   callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2, Prob.QP.F,...
   LogFile, SaveFile, SaveMode, qc, iisRequest, iisFile, saRequest, advbas);


%xprint(x,'x:');
%xprint(slack,'s:');
%xprint(v,'v:');
%xprint(rc,'rc:');
%xprint(basis(1:m),'cbasis:');
%xprint(basis(m+1:m+n),'xbasis:');

% Exit texts, depending on Inform
switch Inform
   case 1, Text='Optimal solution found';
   case 2, Text='Model has an unbounded ray';
   case 3, Text='Model has been proved infeasible';
   case 4, Text='Model has been proved either infeasible or unbounded';
   case 5, Text='Optimal solution is available, but with infeasibilities after unscaling';
   case 6, Text='Solution is available, but not proved optimal, due to numeric difficulties';
   case 10, Text='Stopped due to limit on number of iterations';
   case 11, Text='Stopped due to a time limit';
   case 12, Text='Stopped due to an objective limit';
   case 13, Text='Stopped due to a request from the user';
   case 20, Text='Model has an unbounded optimal face';
   case 21, Text='Stopped due to a limit on the primal objective';
   case 22, Text='Stopped due to a limit on the dual objective';
   case 32, Text='Converged, dual feasible, primal infeasible';
   case 33, Text='Converged, primal feasible, dual infeasible';
   case 34, Text='Converged, primal and dual infeasible';
   case 35, Text='Primal objective limit reached';
   case 36, Text='Dual objective limit reached';
   case 37, Text='Primal has unbounded optimal face';
   case 38, Text='Non-optimal solution found, primal-dual feasible';
   case 39, Text='Non-optimal solution found, primal infeasible';
   case 40, Text='Non-optimal solution found, dual infeasible';
   case 41, Text='Non-optimal solution found, primal-dual infeasible';
   case 42, Text='Non-optimal solution found, numerical difficulties';
   case 43, Text='Barrier found inconsistent constraints';
      
   case 101, Text='Optimal integer solution found';
   case 102, Text='Optimal sol. within epgap or epagap tolerance found';
   case 103, Text='Solution is integer infeasible';
   case 104, Text='The limit on mixed integer solutions has been reached';
   case 105, Text='Node limit exceeded, integer solution exists';
   case 106, Text='Node limit exceeded, no integer solution';
   case 107, Text='Time limit exceeded, integer solution exists';
   case 108, Text='Time limit exceeded, no integer solution';
   case 109, Text='Terminated because of an error, but integer solution exists';
   case 110, Text='Terminated because of an error, no integer solution';
   case 111, Text='Limit on tree memory has been reached, but an integer solution exists';
   case 112, Text='Limit on tree memory has been reached; no integer solution';
   case 113, Text='Stopped, but an integer solution exists';
   case 114, Text='Stopped; no integer solution';
   case 115, Text='Problem is optimal with unscaled infeasibilities';
   case 116, Text='Out of memory, no tree available, integer solution exists';
   case 117, Text='Out of memory, no tree available, no integer solution';
   case 118, Text='Model has an unbounded ray';
   case 119, Text='Model has been proved either infeasible or unbounded';
    
    otherwise, Text='Unknown CPLEX Status value';

end
Result.ExitText = Text;


% Exitflags, depending on Inform
switch(Inform)
   
   case {1,11,32,33,34,101,102,115} % Successful
      ExitFlag = 0;
      
   case {4,5,6,7,8,9} % Time/Iterations limit exceeded
      ExitFlag = 1;
      
   case {3,19,37,} % Unbounded
      ExitFlag = 2;
      
   case {2,14,15,16,103} % Infeasible
      ExitFlag = 4;
      
   case {} % Input errors
      ExitFlag = 10;
      
   case {116,117,118,119} % Memory errors
      ExitFlag = 11;
      
   otherwise % Other Inform values
      ExitFlag = -1;
end

Result.ExitFlag = ExitFlag;

% Result printing
if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / CPLEX solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   
   fprintf('%s',Text);
   fprintf('\n');
   
   fprintf('\nObjective function at x (obj) %25.16f\n\n',f_k);
   if MIP
      fprintf('LP iterations%7d. ',lpiter);
   else
      fprintf('Nodes visited%7d. ',glnodes);
   end
   fprintf('\n');
   fprintf('Number of infeasibilities%7d. ',ninf);
   fprintf('Sum of infeasibilities %e. ',sinf);
   fprintf('\n');
end


if PriLev > 1
   if isempty(MAX_x)
      MAX_x=length(x);
   end
   fprintf('Optimal x = \n');
   xprinte(x(1:min(length(x),MAX_x)),'x:  ');
end


if PriLev > 3
   if isempty(MAX_c)
      MAX_c=20;
   end
   fprintf('Dual variables (Lagrangian multipliers) v_k = \n');
   xprinte(v(1:min(length(v),MAX_c)),'v_k:');
   
   fprintf('Reduced costs rc: Last %d elements should be -v_k\n',length(v));
   xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'rc: ',' %14.9f',5);
end

Result.f_k=f_k;
Result.x_0=[];
Result.x_k=x;
Result.v_k=[rc;v];

%Result.g_k=[];

if ~isempty(c)
   if ~Fzero
      Result.g_k=Prob.QP.F*x+c;
   else
      Result.g_k=c;
   end
elseif Fzero
   Result.g_k=[];
else
   Result.g_k=Prob.QP.F*x;
end

Result.c_k=[];
Result.cJac=[];

if length(basis) ~= nTot, basis = zeros(nTot,1); end

% State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
xbasis=basis(m+1:m+n);
xbasis(xbasis==0)=-1;
xbasis(xbasis==1)=0;
xbasis(xbasis==-1)=1;

Result.xState = xbasis;

cbasis = basis(1:m);
cbasis(cbasis==0)=-1;
cbasis(cbasis==1)=0;
cbasis(cbasis==-1)=1;

Result.bState = cbasis;

if PriLev > 2
   fprintf('State vector xState for x = \n');
   xprinti(xbasis(1:min(length(xbasis),MAX_x)),'xS: ');
end
if MIP
   Result.Iter     = glnodes; % Number of nodes visited
else
   Result.Iter     = lpiter;  % Simplex iterations
end
Result.FuncEv   = lpiter;
if ~Fzero
   Result.GradEv   = lpiter;
else
   Result.GradEv   = 0;
end

Result.ConstrEv = lpiter;


if 0
   Result.Inform   = glstatus;
else
   Result.Inform   = Inform;
end

%xTol=Prob.optParam.xTol;

B       = xbasis;
B(B==2) = -1;
Result.QP.B=B;

Result.MIP.ninf=ninf;
Result.MIP.sinf=sinf;
Result.MIP.slack=slack;
Result.MIP.lpiter=lpiter;
Result.MIP.glnodes=glnodes;
Result.MIP.basis=basis;

Result.CPLEX.iis = iis;
Result.CPLEX.sa = sa;

% Not needed to save callback info
% Result.MIP.cpxCBInfo = cpxCBInfo;

% Return information from CPLEX
Result.MIP.cpxRetVec = cpxRetVec;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 020818 hkh  Last revision for xpress
% 020805 fhe  Modified xpressTL to be used with CPLEX as cplexTL
% 020922 hkh  Revision, using cpxControl for CPLEX parameters
% 030731 ango ExitFlag now conforms with TOMLAB standard, comments revision
% 030930 ango Added LogFile and SaveFile features
% 031014 ango Better check on SaveMode and SaveFile values, also commented.
% 031125 ango Adapt status messages and codes to CPLEX 9.0 
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040113 ango Comments revised for quadratic constraints
% 040528 hkh  Avoid lowering ITLIM, unless user has set optParam.MaxIter
% 040528 hkh  Comment out all DEBUG parts
% 040803 med  Added pragmas for MATLAB Compiler
% 040818 fhe  Added iis and sa and some help text for iis and sa.
% 040825 ango Slight change for sa help text.
% 041213 hkh  Use BIG as Prob.BIG is nonempty, otherwise 1E12
% 041213 hkh  Report glnodes as Result.Iter for MIP problems
% 041213 hkh  Use lpiter for Result.FuncEv, GradEv, ConstrEv 
% 050201 ango Add support for advanced basis
