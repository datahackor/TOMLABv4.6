% TOMLAB /CPLEX MILP, MIQP, QP and LP Solver
%
% -----------------------------------------------------------------------
%
%   cplex.m solves the following 
%   mixed integer (linear or quadratic) programming problem (MILP, MIQP):
%
%   minimize   0.5 * x'*F*x + c'x     subject to:	
%      x             x_L <=    x   <= x_U
%                    b_L <=   Ax   <= b_U
%   where 
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the CPLEX sparse matrix format.
%   c, x_L, x_U has dimension n
%   b_L, b_U has dimension m
%   F is a n x n symmetric matrix, sparse or dense. 
%   If F is empty, an LP or MILP problem is solved
%
%   CPLEX 9.0 also handle quadratic constraints (QC) for the continuous case. 
%   QC's are defined as
%
%               x'*Qj*x + aj'*x <R> rj   , j=1,2,...,nq
%
%   where <R> is one and only one of the relations <= or >=.
%
%   The aj vector may be zero for any given QC, but the Qj matrix must have at
%   least one nonzero element.
%
% ------------------------------------------------------------------------
%
% function [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
%          glnodes, iis, sa] = cplex(c, A, x_L, x_U, b_L, b_U, ...
%          cpxControl, callback, PriLev, Prob, IntVars, PI, SC, SI, ...
%          sos1, sos2, F, logfile, savefile, savemode, qc, ...
%          iisRequest, iisFile, saRequest, basis);
%
% INPUT:  
% c         Linear objective function cost coeffs, n x 1.
% A         Linear constraint matrix, dense or sparse m x n matrix.
%           cplex.m converts the matrix to a sparse format.
% x_L       Lower bounds on x. (if [], then x_L=0 assumed)
% x_U       Upper bounds on x
% b_L       Lower bounds on linear constraints
%
%    The following parameters are optional: 
%
% b_U        Upper bounds on linear constraints (if [], then b_U=b_L assumed)
% cpxControl Structure, where the fields are set to the CPLEX 
%            parameters that the
%            user wants to specify values for. The prefix CPX_PARAM
%            is not used. For example:
%
%  cpxControl.STARTALG = 1 Solve root node with Primal instead of Dual simplex
%  cpxControl.SUBALG   = 4 Solve subnodes with Barrier with crossover
%
% callback  Logical vector defining which callbacks to use in CPLEX
%  If the i:th entry of logical vector callback is set, the corresponding 
%  callback is defined. See CPLEX manuals and the Tomlab /CPLEX User's Guide
%  The callback calls the m-file specified below. The user may edit this file,
%  or make a new copy, which is put before in the Matlab path.
%
% callback(1)  cpxcb_PRIM.m       Primal simplex callback
%         (2)  cpxcb_DUAL.m       Dual simplex callback
%         (3)  cpxcb_PRIMCROSS.m  Primal crossover callback 
%         (4)  cpxcb_DUALCROSS.m  Dual crossover callback
%         (5)  cpxcb_BARRIER.m    Barrier log callback
%         (6)  cpxcb_PRESOLVE.m   Presolve callback 
%         (7)  cpxcb_MIP.m        MIP callback
%         (8)  cpxcb_MIPPROBE.m   MIP probe or clique merging callback
%         (9)  cpxcb_FRACCUT.m    Gomory fractional cut callback
%         (10) cpxcb_DISJCUT.m    Disjunctive cut callback 
%         (11) cpxcb_FLOWMIR.m    Mixed integer rounding cut callback
%
% PriLev    Printing level in the CPLEX m-file and CPLEX C-interface.
%           = 0    Silent
%           = 1    Warnings and Errors
%           = 2    Summary information
%           = 3    More detailed information
%
%           > 10   Pause statements, and maximal printing (debug mode)
%
% Prob      A structure. If TOMLAB calls CPLEX, then Prob is the standard
%           TOMLAB problem structure, otherwise the user optionally may set:
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%           If any callback is defined (see description of callback) then
%           problem arrays are set as fields in Prob, and the Prob structure
%           is always passed to the callback routines as the last parameter.
%           The defined fields are Prob.QP.c, Prob.QP.F, 
%           Prob.x_L, Prob.x_U, Prob.A, Prob.b_L, Prob.b_U.
%           (if input is [], then Prob.P=1 is set)
%          
%           NOTE: The cplex MEX is using values >= 1E12 to define infinite
%           values in lower/upper bounds of variables and constraints
%           By setting Prob.BIG to a nonempty positive big value, this value
%           will be used by the MEX. 
%           DO NOT USE the MATLAB inf or NaN value in any arrays.
%           Convert any inf value to 1E12, -inf to -1E12 (or the value set
%           as Prob.BIG).
%
% IntVars   Defines which variables are integers, of general type I or binary B
%           Variable indices should be in the range [1,...,n].
%           IntVars is a single integer ==> Variable 1:IntVars are integer 
%           IntVars is a logical vector ==> x(find(IntVars > 0)) are integers 
%           IntVars is a vector of indices ==> x(IntVars) are integers 
%           (if [], then no integers of type I or B are defined)
%           CPLEX checks which variables has x_L=0 and x_U=1, i.e. binary.
%
% PI        Integer variables of type Partially Integer (PI), i.e. takes an
%           integer value up to a specified limit, and any value above that
%           limit.
%           PI must be a structure array where:
%           PI.var  Vector of variable indices in the range [1,...,n]
%           PI.lim  A vector of limit values, for each of the variables
%                   specified in PI.var, i.e. for variable i, 
%                   that is the PI variable with index j in PI.var:
%                   x(i) takes integer values in [x_L(i),PI.lim(j)] and
%                   continous values in [PI.lim(j),x_U(i)].
%
% SC        A vector with indices for the integer variables of type 
%           Semi-continuous (SC), i.e. that takes either the value 0 or a 
%           real value in the range [x_L(i),x_U(i)], assuming for some j,
%           i = SC(j), where i is an variable number in the range [1,n].
%
% SI        A vector with indices for the integer variables of type 
%           Semi-integer (SI), i.e. that takes either the value 0 or 
%           an integer value in the range [x_L(i),x_U(i)], assuming for some j,
%           i = SI(j), where i is an variable number in the range [1,n].
%
% sos1      A structure defining the Special Ordered Sets of Type One (sos1). 
%           Assume there are k sets of type sos1, then
%           sos1(1).var is a vector of indices for variables in sos1, set 1.
%           sos1(1).row is the row number for the reference row identifying
%                       the ordering information for the sos1 set, i.e.
%                       A(sos1(1).row,sos1(1).var) identifies this information
%           sos1(2).var is a vector of indices for variables in sos1, set 2.
%           sos1(2).row is the row number for the reference row of sos1 set 2.
%           ...
%           sos1(k).var is a vector of indices for variables in sos1, setk.
%           sos1(2).row is the row number for the reference row of sos1 set k.
%
% sos2      A structure defining the Special Ordered Sets of Type Two (sos2). 
%           Specified exactly as sos1 sets, see sos1 input variable description
%
% F         Square dense or sparse matrix. Empty if non-quadratic problem.
%
% logfile   Name of file to write CPLEX log to. If empty, no log is
%           written. 
%
% savefile  Name of file to write CPLEX problem just prior to calling the 
%           CPLEX solver. If empty, nothing is written. Also see the savemode 
%           parameter below.
%
% savemode  Integer flag indicating which format to use for the save file.
%           The following values are possible:
%
%        1: SAV  Binary SAV file
%        2: MPS  MPS format, original format
%        3: LP   CPLEX LP format, original format
%        4: RMP  MPS file with generic names
%        5: REW  MPS file with generic names
%        6: RLP  LP  file with generic names
%
%            The SAV format is a binary format suitable for submission to
%            ILOG help desk.
%
% qc       Structure array defining quadratic constraints ("qc").
%
%           Please note that CPLEX 9.0 only handles single-sided bounds 
%           on qc's. An arbitrary number of qc's is set using the Prob.QP.qc
%           structure array:
%
%          qc(1).Q   = sparse( <quadratic coefficient nxn matrix> );
%          qc(1).a   = full  ( <linear coefficient nx1 vector   > );
%          qc(1).r_U = <scalar upper bound> ;  
%
%           And similarly for qc(2), ... , qc(n_qc).
%
%           The standard interpretation is x'*Q*x + c'*x <= r_U, but it is 
%           possible to define an alternative sense x'*Q*x + c'*x >= r_L 
%           by setting qc(i).sense to a nonzero value and specifying a
%           lower bound in qc(i).r_L.
%
%           Observe that the Q matrix must be sparse, non-empty and positive 
%           semi-definite for all qc's. The linear coefficient vector qc(i).a 
%           may be omitted or set empty, in which case all zeros are assumed. 
%
%           Likewise, if a bound r_U or r_L is empty or not present, it
%           is assumed to be 0.0. Note that this is contrary to the usual
%           Tomlab standard, where an empty or omitted bound is assumed
%           to be +/- Inf. The reason is that a single-sided constraint with 
%           an infinite bound would have no meaning. 
%
% iisRequest    Flag indicating whether to compute an Irreducible
%               Infeasible Set (IIS) and return it to MATLAB. This
%               option can only be set for an LP problem being
%               solved using primal or dual simplex, or barrier with
%               crossover.
%               
%               = 0     Don't return IIS to MATLAB (default).
%               = 1     Compute IIS and return it to MATLAB if an LP problem 
%                       has been proven infeasible.
%
%               The IIS is returned through the output parameter 'iis'.
%
% iisFile       Name of file to write the IIS to. No file is written if
%               this input parameter is empty or if no IIS is available.
%
% saRequest     Structure telling whether and how you want CPLEX to perform a
%               sensitivity analysis (SA). You can complete an SA on the 
%               objective function, right hand side vector, lower and 
%               upper bounds. The saRequest structure contains four
%               sub structures:
%
%                   .obj, .rhs, .xl, .xu
%               
%               Each one of these contain the field:
%
%                   .index
%
%               .index contain indices to variables or constraints 
%               of which to return possible value ranges.
%
%               The .index array has to be sorted, ascending.
%
%               To get an SA of objective function on the four variables 120 
%               to 123 (included) and variable 19, the saRequest structure 
%               would look like this:
%
%                   saRequest.obj.index = [19 120 121 122 123];
%
%               The result is returned through the output parameter 'sa'.
%
% basis         Vector with CPLEX starting basis. 
%               If re-solving a similar problem several times, this can be
%               set to the 'basis' output argument of an earlier call
%               to cplex.m. The length of this vector must be equal to the sum
%               of the number of rows (m) and columns (m).
%
%               The first m elements contain row basis information, with the
%               following possible values for non-ranged rows:
%
%             0 associated slack/surplus/artificial variable nonbasic at value 0.0 
%             1 associated slack/surplus/artificial variable basic 
%
%               and for ranged rows (both upper and lower bounded)
% 
%             0 associated slack/surplus/artificial variable nonbasic at its lower bound  
%             1 associated slack/surplus/artificial variable basic 
%             2 associated slack/surplus/artificial variable nonbasic at its upper bound  
%
%
%               The last n elements, i.e. basis(m+1:m+n) contain column 
%               basis information:
%
%             0 variable at lower bound 
%             1 variable is basic 
%             2 variable at upper bound 
%             3 variable free and nonbasic 
%
% ------------------------------------------------------------------------------
%
% OUTPUT: 
%
% x         Solution vector with decision variable values (n x 1 vector)
% slack     Slack variables (m x 1 vector)
% v         Lagrangian multipliers (dual solution vector) (m x 1 vector)
% rc        Reduced costs. Lagrangian multipliers for simple bounds on x.
% f_k       Objective function value at optimum
%
% ninf      Number of infeasibilities
% sinf      Sum of infeasibilities
%
% Inform    Result of CPLEX run: (S=Simplex, B=Barrier) 
%   
%       1 (S,B) Optimal solution is available
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
%      32 Converged, dual feasible, primal infeasible
%      33 Converged, primal feasible, dual infeasible
%      34 Converged, primal and dual infeasible
%      35 Primal objective limit reached
%      36 Dual objective limit reached
%      37 Primal has unbounded optimal face
%      38 Non-optimal solution found, primal-dual feasible
%      39 Non-optimal solution found, primal infeasible
%      40 Non-optimal solution found, dual infeasible
%      41 Non-optimal solution found, primal-dual infeasible
%      42 Non-optimal solution found, numerical difficulties
%      43 Barrier found inconsistent constraints
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
%
% basis     basis status of constraints + variables, (m + n x 1 vector)
%
% lpiter    Number of LP iterations
%
% glnodes   Number of nodes visited
%
% iis       Structure of information about an IIS if requested. The fields:
%
%               iisStatus   Status information. Positive on
%                           success, negative on failure. Possible 
%                           negative values and their meaning:
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
% sa        Structure with information about the requested SA, if requested.
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
% basis     Vector containing basis at solution
%
%           The first m elements contain row basis information, with the
%           following possible values for non-ranged rows:
%
%         0 associated slack/surplus/artificial variable nonbasic at value 0.0 
%         1 associated slack/surplus/artificial variable basic 
%
%           and for ranged rows (both upper and lower bounded)
% 
%         0 associated slack/surplus/artificial variable nonbasic at its lower bound  
%         1 associated slack/surplus/artificial variable basic 
%         2 associated slack/surplus/artificial variable nonbasic at its upper bound  
%
%
%           The last n elements, i.e. basis(m+1:m+n) contain column 
%           basis information:
%
%         0 variable at lower bound 
%         1 variable is basic 
%         2 variable at upper bound 
%         3 variable free and nonbasic 
%
% -----------------------------------------------------------------------
%
% NOTE! After each call of cplex, output information is defined in a
% vector called cpxRetVec, to access this global vector do
%        global cpxRetVec
%
% See the TOMLAB /CPLEX User's Guide for a description of each element
%
% -----------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 9.0.2$
% Written Aug 5, 2002.    Last modified Feb 1, 2005.
%

function [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, iis, sa] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, ...
    PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2, F, ...
    logfile,savefile,savemode,qc, iisRequest, iisFile, saRequest, basis)

%#function cpx2cbinfo cpx2retvec
if nargin < 25
   basis = [];
if nargin < 24
 saRequest = [];
 if nargin < 23
  iisFile = '';
  if nargin < 22
   iisRequest = [];
   if nargin<21    
    qc = [];
    if nargin < 20
     savemode = [];
     if nargin < 19
      savefile = '';
      if nargin < 18
       logfile = '';
       if nargin < 17
        F = [];
        if nargin < 16
         sos2 = [];
         if nargin < 15
          sos1 = [];
          if nargin < 14
           SI = [];
           if nargin < 13
            SC = [];
            if nargin < 12
             PI = [];
             if nargin < 11
              IntVars = [];
              if nargin < 10
               Prob = [];
               if nargin < 9
                PriLev = [];
                if nargin < 8
                 callback = [];
                 if nargin < 7
                  cpxControl = [];
                  if nargin < 6
                   b_U = [];
                   if nargin < 5
                    error('cplex needs at least 5 arguments');
end, end, end, end, end, end, end, end, end, end, end, end, end, end , end, end
end, end, end, end, end

global cpxParameters cpxProblemAttributes

% Set empty to avoid old values from last run to still be used
cpxParameters        = [];
cpxProblemAttributes = [];

DEBUG=0;

n = max( [length(c),size(F,1),size(A,2),length(x_L),length(x_U)] );
if n==0, error('cplex: cannot determine problem dimension'); end

b_L = b_L(:);
b_U = b_U(:);
x_L = x_L(:);
x_U = x_U(:);

if isempty(x_L),    x_L = zeros(n,1); end
if isempty(b_U),    b_U = b_L; end
if isempty(PriLev), PriLev = 0; end
if isempty(callback) 
   callback=zeros(12,1); 
%   callback(11)=1;  % Default use the User Output Callback 
end

[m,nA] = size(A);

% Error checking 

if nA ~= n & m > 0
   fprintf('n = %d. Number of columns in A = %d\n',n,nA);
   error('cplex: Illegal length of A');
end

if length(b_L)~=m
   fprintf('m = %d. Length of b_L = %d\n',m,length(b_L));
   error('cplex: Illegal length of b_L');
end
if length(b_U)~=m
   fprintf('m = %d. Length of b_U = %d\n',m,length(b_U));
   error('cplex: Illegal length of b_U');
end
if length(x_L)~=n
   fprintf('n = %d. Length of x_L = %d\n',n,length(x_L));
   error('cplex: Illegal length of x_L');
end
if length(x_U)~=n
   fprintf('n = %d. Length of x_U = %d\n',n,length(x_U));
   error('cplex: Illegal length of x_U');
end

if ~isempty(IntVars) | ~isempty(sos1) | ~isempty(sos2) | ~isempty(PI) | ...
   ~isempty(SC) | ~isempty(SI)
   MIP=1;
else
   MIP=0;
end

% if ~isempty(qc) & MIP==1
%   error('cplex: MILP/MIQP problems cannot have quadratic constraints');
% end


% Vector for integers indicating type
iv = zeros(n,1);

if MIP
   callback(8)=0;  % Avoid simplex log callback for MIP

   if isempty(IntVars)
      % No binary variables B or integer variables of type I
   elseif length(IntVars)==1
      IntVars = min(IntVars,n);
      iv(1:IntVars) = 1;
   elseif any(IntVars==0) | all(IntVars==1)
      % Assume binary logical vector given
      iv(1:length(IntVars)) = IntVars;
   else
      if any(IntVars < 1 | IntVars > n)
         error('cplex: Illegal IntVars vector');
      end
      iv(IntVars)=1;
   end

   % Semi-continuous variables, o and real interval
   if ~isempty(SC)
      if any(SC < 1 | SC > n)
         SC
         error('cplex: Illegal index in SC vector');
      end
      iv(SC) = 3;
   end
   % Semi-integer variables, 0 and integer interval
   if ~isempty(SI)
      if any(SI < 1 | SI > n)
         SI
         error('cplex: Illegal index in SI vector');
      end
      iv(SI) = 4;
   end
   % Partially integer variables
   if ~isempty(PI)
      if isfield(PI,'var')
         iv(PI.var) = 5;
      else
         PI.var = [];
      end
   else
      PI.var = [];
   end

   glcolidx = find(iv);

   ngl=length(glcolidx);
   if ngl > 0
      % Identify binary variables among the general Integer variables
      ix = find(x_L(glcolidx)==0 & (x_U(glcolidx)==1) & (iv(glcolidx)==1));
      if ~isempty(ix)
         iv(glcolidx(ix)) = 2;
      end
   end

   %if ngl > 0
 
      %gltype = iv(glcolidx);
      %gltype(gltype==1)=0; % Set as general Integer variables

      %% Identify binary variables among the general Integer variables
      %ix = find(x_L(glcolidx)==0 & (x_U(glcolidx)==1) & (gltype==0));

      %gltype(ix)     = 1;

      %gllim          = x_U(glcolidx);

      %% Partially integer variables
      %ix = find(gltype==2);
      %if ~isempty(ix)
      %   gllim(ix)  = PI.lim;
      %end

      %if DEBUG
      %   xprinti(gltype,'gltype:');
      %   xprinte(gllim,'gllim:');
      %end
   %else
      %gltype=[];
      %gllim=[];
   %end

   if isempty(sos1) & isempty(sos2)
      settype   = [];
      setbeg    = [];
      setcolidx = [];
      setref    = [];
   else
      ns1       = length(sos1);
      ns2       = length(sos2);
      settype   = [ones(ns1,1);2*ones(ns2,1)];
      setbeg    = 1;
      setcolidx = [];
      setref    = [];
      for i=1:ns1
          if ~isfield(sos1(i),'var')
             fprintf('sos1 set %d. ',i);
             error('sos1 field var is missing');
          end
          if ~isfield(sos1(i),'row')
             fprintf('sos1 set %d. ',i);
             error('sos1 field row is missing');
          end
          ix        = sos1(i).var;
          if ~(all(ix >= 1 & ix <= n))
             fprintf('sos1 set %d. ',i);
             fprintf('\n');
             error('cplex: Illegal sos1 input variable vector');
          end
          row       = sos1(i).row;
          if ~(row >= 0 & row <= m)
             fprintf('sos1 set %d. ',i);
             fprintf('Illegal row number  %d.',row);
             fprintf('\n');
             error('cplex: Illegal sos1 row data');
          end
          k         = length(ix);
          setbeg    = [setbeg; setbeg(length(setbeg))+k];
          setcolidx = [setcolidx; ix(:)];
          if row==0
             key = full(c(ix));
             %key = full(c(ix)) + 0.05*rand(length(ix),1);
          else
             key = full(A(row,ix))';
             %key = full(A(row,ix))' + 0.05*rand(length(ix),1);
          end
          % Sort the ordering key, to get unique sequence of variables
          [zz,skey]=sort(key);
          setref = [setref; skey(:)];
      end
      if DEBUG
         xprinti(setbeg,'setbeg');
         xprinti(setref,'setref');
         xprinti(setcolidx,'setcolidx');
         pause
      end
      for i=1:ns2
          if ~isfield(sos2(i),'var')
             fprintf('sos2 set %d. ',i);
             error('sos2 field var is missing');
          end
          if ~isfield(sos2(i),'row')
             fprintf('sos2 set %d. ',i);
             error('sos2 field row is missing');
          end
          ix        = sos2(i).var;
          if ~(all(ix >= 1 & ix <= n))
             fprintf('sos2 set %d. ',i);
             fprintf('\n');
             error('cplex: Illegal sos2 input variable vector');
          end
          row       = sos2(i).row;
          if ~(row >= 0 & row <= m)
             fprintf('sos2 set %d. ',i);
             fprintf('Illegal row number  %d.',row);
             fprintf('\n');
             error('cplex: Illegal sos2 row data');
          end
          k         = length(ix);
          setbeg    = [setbeg; setbeg(length(setbeg))+k];
          setcolidx = [setcolidx; ix(:)];
          if row==0
             key = full(c(ix));
          else
             key = full(A(row,ix))';
          end
          % Sort the ordering key, to get unique sequence of variables
          [zz,skey]=sort(key);
          setref = [setref; skey(:)];
      end
   end
else
%   callback([1:7,9])=0;
   gltype    = []; glcolidx  = []; gllim     = []; settype   = [];
   setbeg    = []; setcolidx = []; setref    = [];
end

% Check that Prob.P is set.
if isempty(Prob)
   Prob = struct('P',1);
end
if ~isfield(Prob,'P')
   Prob.P=1;
end

%PriLev=2
%callback(11)=0
%callback(1:11)=1
%callback(9)=1

if any(callback)
   % Define fields in Prob for use in callback m-files
   
   disp('cplex.m: setup prob');
   Prob.A    = sparse(A);
   Prob.QP.c = c;
   Prob.QP.F = sparse(F);
   Prob.QP.qc = qc;
   Prob.x_L  = x_L;
   Prob.x_U  = x_U;
   Prob.b_L  = b_L;
   Prob.b_U  = b_U;

   %Prob.qgtype    = gltype;
   %Prob.mgcols    = glcolidx;
   %Prob.dlim      = gllim;
   Prob.iv        = iv;
   Prob.qstype    = settype;
   Prob.msstart   = setbeg;
   Prob.mscols    = setcolidx;
   Prob.dref      = setref;
end


% The cplex MEX has 1E12 as default, but is changed if Prob.BIG is noempty
BIG=DefPar(Prob,'BIG',1E12);

%PriLev = 15;

% iv vector indicating integer variable type:
% 0    Continuous
% 1    General integer
% 2    Binary
% 3    Semi-continuous variables, 0 and real interval
% 4    Semi-integer variables, 0 and integer interval
% 5    Partially integer variables (not in CPLEX, still used now)
% callback(5) = 1;
% callback
% PriLev = 0

[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, ...
    iis, sa] ...
   = cplexmex(c, sparse(F), sparse(A), callback , x_L, x_U, b_L, b_U,...
   BIG, Prob.P, PriLev, Prob, ...
   iv, settype, setbeg, setcolidx, setref, cpxControl, [], [], ...
   logfile, savefile, savemode, qc, iisRequest, iisFile, saRequest, basis);
%gltype, glcolidx, gllim, settype, setbeg, setcolidx, setref, cpxControl);

function x=DefPar(s,f,xdef)

if(nargin < 3)
   xdef = [];
end

if isfield(s,f)
   x = getfield(s,f);
   if isempty(x), x = xdef; end
else
   x = xdef;
end

% MODIFICATION LOG:
%
% 020805 fhe  Modified from xpress.m  to be used with CPLEX
% 020806 hkh  Revision for CPLEX use, change integer variable handling
% 020922 hkh  Revision of comments
% 030731 ango Correction of comments on Inform
% 030925 ango New dump mode (LP) added; logfile feature added. 
% 030930 ango Removed previous handling of problem dump
% 031125 ango CPLEX 9.0 adaptations. Error codes reviewed. 
% 040101 hkh  Wrong computation of n, must use x_L,x_U, not b_L and b_U
% 040803 med  Added pragmas for MATLAB Compiler
% 040818 fhe  Added iis and sa and some help text for iis and sa.
% 040825 ango Slight change for sa structure fields
% 041213 hkh  Set BIG hard coded as 1E12, if not Prob.BIG is set
% 041213 hkh  Add comments about BIG and inf
% 050128 med  Function DefPar added
% 050201 ango Basis support added
