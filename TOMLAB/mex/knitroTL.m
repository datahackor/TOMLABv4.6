% TOMLAB KNITRO 4.0 NLP Solver
%
% function Result = knitroTL(Prob)
%
% INPUT:  
%
% Prob   Problem structure in TOMLAB format
%
% -----------------------------------------------
% Fields used in input structure Prob:
%
% x_0         Initial vector.
%
% x_L, x_U    Bounds on variables. 
%
% A           Linear constraint matrix. 
% b_L, b_U    Bounds on linear constraints. 
%
% c_L, c_U    Bounds on nonlinear constraints. 
%
% LargeScale  If 1, indicates that the problem should be treated as sparse.
%             This makes knitroTL send a sparse version of Prob.A to the
%             solver. To avoid poor performance, sparse ConsPattern and
%             d2LPattern should be given. 
%             (Default 1) 
%
% ConsPattern Sparsity pattern of the gradient of the nonlinear
%             constraints. Should be sparse if given. 
%
% d2LPattern  Sparsity pattern of the Hessian of the Lagrangian function. 
%
% f_Low       A lower bound on the objective function value. KNITRO will stop 
%             if it finds an x for which f(x)<=f_Low.
%
% PriLevOpt   Print level in solver. TOMLAB /KNITRO supports solver
%             progress information being printed to the MATLAB command
%             window as well as to a file.
%
%             For output only to the command window: set PriLevOpt
%             to a positive value (1-6) and set:
%
%          >> Prob.KNITRO.PrintFile = '';  
%
%             If PrintFile is set to a valid filename, the same
%             information will be written to the file (if PriLevOpt
%             is nonzero).
%
%             For printing only to the PrintFile, set a name (or
%             knitro.txt will be used) and a _negative_ PriLevOpt
%             value. For example: 
%
%          >> Prob.KNITRO.PrintFile = 'myprob.txt';
%          >> Prob.PriLevOpt        = -2;
%
%             To run TOMLAB /KNITRO silently, set PriLevOpt=0; (default)
%
% -----------------------------------------------
% Fields used in Prob.KNITRO:
% -----------------------------------------------
%
% PrintFile   Name of file to print solver progress and status
%             information to. Please see PriLevOpt parameter
%             described above. 
%
% options     Structure with options values for KNITRO. The following
%             subfields may be set:
%
%  ALG        Choice of algorithm 
%	      
%             0 - Automatic selection (default)
%             1 - Barrier/Direct
%             2 - Barrier/CG 
%             3 - SLQP 
%	      
%  BARRULE    Barrier parameter strategy. For ALG=1, all strategies
%             are available. For ALG=2, only 0-2 can be used. See
%             the TOMLAB /KNITRO User's Guide for detailed information.
% 	      
%             0 - Automatic (default)
%             1 - Monotone
%             2 - Adaptive
%             3 - Probing
%             4 - Damp MPC
%             5 - Full MPC
%	      
%  FEASIBLE   Indicates whether or not to use the feasible version
%             of KNITRO.
%	      
%             0 - Iterates may be infeasible (default).
%	      
%             1 - Given an initial point which sufficiently satisfies
%                 all inequality constraints as defined by
% 	      
%                     cl + tol <= c(x) <= cu - tol        (*)
%	      
%                 for cl(i) ~= cu(i), the feasible version of KNITRO
%                 ensures that all subsequent iterates strictly
%                 satisfy the inequality constraints. However, the
%                 iterates may not be feasible with respect to the
%                 equality constraints.
%	      
%                 If the initial point is infeasible (or not
%                 sufficiently feasible according to (*)) with respect
%                 to the inequality constraints, then KNITRO will run
%                 the infeasible version until a point is obtained
%                 which sufficiently satisfies all the inequality
%                 constraints. At this point it will switch to
%                 feasible mode.
%	      
%                 The tolerance 'tol' in (*) is determined by the
%                 parameter FEASMODETOL described below. 
%	      
%             NOTE: The FEASIBLE option can be used only with ALG=2.
%	      
%  GRADOPT    Specifies how to calculate the gradients of the objective
%             and constraint functions. In addition to the available
%             numeric differentiation methods available in Tomlab,
%             KNITRO can internally estimate forward or centered
%             numerical derivatives. 
% 	       
%             The following settings are available:
% 	       
%             1 - KNITRO expects the user to provide exact first derivatives (default).  
% 	       
%                 However, note that this can imply Tomlab numerical derivatives,
%                 calculated by any of the available methods. See the Tomlab 
%                 Users Guide for more information on numerical derivatives.  
% 	       
%             2 - KNITRO will estimate forward finite differential
%                 derivatives.
% 	       
%             3 - KNITRO will estimate centered finite differential
%                 derivatives.
% 	       
%             4 - KNITRO expects exact first derivatives and performs gradient 
%                 checking by internally calculating forward finite differences. 
% 	       
%             5 - KNITRO expects exact first derivatives and performs gradient 
%                 checking by internally calculating centered finite
%                 differences. 
%	      
%  HESSOPT    Specifies how to calculate the Hessian of the Lagrangian function. 
% 	       
%             1 - KNITRO expects the user to provide the exact Hessian (default).
% 	       
%             2 - KNITRO will compute a (dense) quasi-Newton BFGS Hessian.
% 	       
%             3 - KNITRO will compute a (dense) quasi-Newton SR1 Hessian.
% 	       
%             4 - KNITRO will compute Hessian-vector products using finite differences.
% 	       
%             5 - KNITRO expects the user to provide a routine to
%                 compute the Hessian-vector products. If this option is
%                 selected, the calculation is handled by the Tomlab
%                 function "nlp_d2Lv.m".
% 	       
%             6 - KNITRO will compute a limited-memory quasi-Newton BFGS.
% 	       
%             NOTE: Typically if exact Hessians (or exact Hessian-
%             vector products) cannot be provided by the user but exact
%             gradients are provided and are not too expensive to
%             compute, option 4 above is recommended. The finite-
%             difference Hessian-vector option is comparable in terms
%             of robustness to the exact Hessian option (assuming exact
%             gradients are provided) and typically not too much slower
%             in terms of time if gradient evaluations are not the
%             dominant cost.
%             
%             However, if exact gradients cannot be provided
%             (i.e. finite-differences are used for the first
%             derivatives), or gradient evaluations are expensive, it
%             is recommended to use one of the quasi-Newton options, in
%             the event that the exact Hessian is not
%             available. Options 2 and 3 are only recommended for small
%             problems (n < 1000) since they require working with a
%             dense Hessian approximation. Option 6 should be used in
%             the large-scale case.
%             
%             NOTE: Options HESSOPT=4 and HESSOPT=5 are not available
%             when ALG=1.
%	      
% HONORBNDS   Determines whether or not to enforce satisfaction of
%             the simple bounds x_L <= x <= x_U throughout the
%             optimization. 
% 	       
%             0 - KNITRO does not enforce that the bounds on the
%                 variables are satisfied at intermediate iterates.
% 	       
%             1 - KNITRO enforces that the initial point and all
%                 subsequent solution estimates satisfy the bounds
%                 on the variables.
% 	       
%  INITPT     Determines whether or not an initial point strategy is
%             used. 
% 	       
%             0 - No initial point strategy is employed. 
% 	       
%             1 - Initial values for the variables are computed. 
% 	       
%  ISLP       These two options are used to tell KNITRO that the 
%  ISQP       problem being solved is a linear/quadratic
%             problem, enabling it to perform certain specializations. 
% 	       
%             Unless specifically set by the user, TOMLAB will set
%             appropriate ISLP and ISQP values.
% 	       
%  MAXCGIT    Determines the maximum allowable number of inner
%             conjugate gradient (CG) iterations per KNITRO minor
%             iteration.
% 	       
%             0 - Automatically determination based on problem
%                 size. (default)
% 	       
%             n - At most n CG iterations may be performed during
%                 one KNITRO minor iteration. n must be positive. 
% 	       
%  MAXIT      Specifies the maximum number of iterations before
%             termination. Unless set by the user, the value from
%             Prob.optParam.MaxIter is used. 
% 	       
%  SCALE      Performs a scaling of the objective and constraint
%             functions based on their values at the initial point. 
% 	       
%             0 - no scaling is performed
% 	       
%             1 - scaling may be used (default) 
%	      
%  SHIFTINIT  Determines whether or not the interior-point
%             algorithm shifts the initial point 
%	      
%             0 - KNITRO will not shift the given initial point to
%                 satisfy the variable bounds before starting the
%                 optimization. 
%	      
%             1 - KNITRO will shift the initial point (default)
%	      
%  SOC        Determines whether or not to use the second order
%             correction (SOC) options. A second order correction
%             may be beneficial for problems with highly nonlinear
%             constraints. 
%	      
%             0 - No second order correction steps are attempted
%	      
%             1 - Second order correction steps may be attempted on
%                 some iterations (default)
%	      
%             2 - Second order corrections are always attempted if
%                 the original setp is rejected and there are
%                 nonlinear constraints. 
%
% -- KNITRO real valued options: 
%
%  DELTA      Specifies the initial trust region radius scaling factor
%             used to determine the initial trust region size.
% 
%             Default value: 1.0
%
%  FEASMODETOL Specifies the tolerance in (*) by which an iterate
%              must be feasible with respect to the inequality
%              constraints before the feasible mode is
%              activated. See the FEASIBLE parameter. 
%
%              Default value: 1.0E-4
%
%  FEASTOL    Specifies the final relative stopping tolerance for the
%             feasibility error. Smaller values of feastol result in a
%             higher degree of accuracy in the solution with respect
%             to feasibility.
%
%             Default value: 1.0e-6
%
%  FEASTOLABS Specifies the final absolute stopping tolerance for the
%             feasibility error. Smaller values of feastol abs result
%             in a higher degree of accuracy in the solution with
%             respect to feasibility.
%
%             Default value: 0.0e0
%
%  OPTTOL     Specifies the final relative stopping tolerance for the
%             KKT (optimality) error.  Smaller values of opttol result
%             in a higher degree of accuracy in the solution with
%             respect to optimality.
%
%             Default value: 1.0e-6
%          
%  OPTTOLABS  Specifies the final absolute stopping tolerance for the
%             KKT (optimality) error. Smaller values of OPTTOLABS
%             result in a higher degree of accuracy in the solution
%             with respect to optimality.
%
%             Default value: 0.0e0
%
%  OBJRANGE   Determines the allowable range of values for the
%             objective function for determining unboundedness. If
%             the magnitude of the objective function is greater
%             than OBJRANGE and the iterate is feasible, then the
%             problem is determined to be unbounded. 
%
%             Default value: 1.0e20
%
%  MAXTIME    Specifies, in seconds, the maximum allowable time before
%             termination. Prob.MaxCPU is used if MAXTIME is not
%             explicitly set by the user.
%
%             Default value: 1e8 
%
%  PIVOT      Specifies the initial pivot threshold used in the
%             factorization routine. The value must be in the range
%             0.0 to 0.5, with higher values resulting in more
%             pivoting (more stable factorization). Values less than
%             0 will be set to 0 and and values larger than 0.5 will
%             be set to 0.5. If pivot is non-positive initially no
%             pivoting will be performed. Smaller values may improve
%             the speed but higher values are recommended for more
%             stability, for example if the problem appears to be
%             very ill-conditioned. 
% 
%             Default value: 1.0e-8
%
%  MU         Specifies the initial barrier parameter value.
%
%             Default value: 1.0e-1
%
% -----------------------------------------------------------------------
%
% Termination Message: At the end of the run a termination message is
% printed indicating whether or not the optimal solution was found and
% if not, why the code terminated. Below is a list of possible
% termination messages and a description of their meaning and
% corresponding info value.
%
%  0: EXIT: OPTIMAL SOLUTION FOUND.
%
%           Knitro found a locally optimal point which satisfies the
%           stopping criterion
%
% -1: EXIT: Iteration limit reached.
%
%           The iteration limit was reached before being able to
%           satisfy the required stopping criteria.
%
% -2: EXIT: Convergence to an infeasible point. Problem may be infeasible.
%
%           The algorithm has converged to a stationary point of
%           infeasibility. This happens when the problem is
%           infeasible, but may also occur on occasion for feasible
%           problems with nonlinear constraints. It is recommended to
%           try various initial points. If this occurs for a variety
%           of initial points, it is likely the problem is infeasible.
%     
% -3: EXIT: Problem appears to be unbounded.
%
%           The objective function appears to be decreasing without
%           bound, while satisfying the constraints.
%
% -4: EXIT: Current point cannot be improved.
%
%           No more progress can be made. If the current point is
%           feasible it is likely it may be optimal, however the
%           stopping tests cannot be satisfied (perhaps because of
%           degeneracy, ill-conditioning or bad scaling).
%    
% -5: EXIT: Current point cannot be improved. Point appears to be
%           optimal, but desired accuracy could not be achieved.
%
%           No more progress can be made, but the stopping tests are
%           close to being satisfied (within a factor of 100) and so
%           the current approximate solution is believed to be
%           optimal.
%
% -50 to -99: Termination values in this range imply some input error.
%
%
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
% H_k      Hessian of the Lagrangian function.
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
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
%
% ExitFlag Exit status.
%
% Inform   KNITRO information parameter.
%
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
% Solver   Name of the solver (knitro).
% SolverAlgorithm  Description of the solver.
%
%
% -------------------------------------------------------------------------
%
% Anders Goran, Tomlab Optimization Inc, E-mail: anders@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written May 5, 2003.  Last modified Dec 10, 2004.
%

function Result = knitroTL(Prob)

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print


if nargin < 1, error('knitroTL needs the Prob structure as input'); end

Prob.solvType = 3; % NLP (CON) solver

KNITRO   = DefPar(Prob,'KNITRO',[]);

options = DefPar(KNITRO,'options',[]);

if isfield(options,'gradopt')
   GRADOPT = options.gradopt;
elseif isfield(options,'GRADOPT')
   GRADOPT = options.GRADOPT;
else
   GRADOPT = 1;
end
if isempty(GRADOPT), GRADOPT = 1; end
if Prob.NumDiff == 6 & ~any(GRADOPT==[2 3]) 
   % Internal differentiation
   options.gradopt = 2;
   GRADOPT = 2;
end

% Hessian information, if required by choice in specs.hessopt
HL = triu(DefPar(Prob,'d2LPattern',[]));

if isfield(options,'hessopt')
   HESSOPT = options.hessopt;
elseif isfield(options,'HESSOPT')
   HESSOPT = options.HESSOPT;
else
   HESSOPT = [];
end
if isempty(HESSOPT) 
   if isempty(Prob.USER.g) | (isempty(Prob.USER.dc) & Prob.mNonLin > 0)
      % Avoid two levels of differentiation, by default
      if Prob.LargeScale == 1 
         if isempty(HL)
            options.hessopt = 6;   % Sparse BFGS
            HESSOPT = 6;
         else
            options.hessopt = 4;   % Sparse Hessian vector products
            HESSOPT = 4;
         end
      else                     
         options.hessopt = 2;   % Dense BFGS
         HESSOPT = 2;
      end
   else
      HESSOPT = 1; 
   end
end

if abs(Prob.NumDiff) == 6 & ~any(HESSOPT==[2 3 4 6]) 
   % Internal differentiation?
   if Prob.LargeScale == 1 
      if isempty(HL)
         options.hessopt = 6;   % Sparse BFGS
         HESSOPT = 6;
      else
         options.hessopt = 4;   % Sparse Hessian vector products
         HESSOPT = 4;
      end
   else                     
      options.hessopt = 2;   % Dense BFGS
      HESSOPT = 2;
   end
end

if any(HESSOPT == [1 5])
   Prob = iniSolve(Prob,3,2,2);
else
   Prob = iniSolve(Prob,3,1,1);
end

Result=ResultDef(Prob);

Result.Solver='KNITRO';
% Result.SolverAlgorithm='Interior Point NLP KNITRO 4.0';

PriLev=Prob.PriLevOpt;

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

BIG=1E20;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);

m  = m1+m2;

% Variable bounds separately
xl = bl(1:n);
xu = bu(1:n);

cl = bl(n+1:end);
cu = bu(n+1:end);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
  if nA~=n,  error('Linear constraints A MUST have n columns!'); end 
  if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end 
end 

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( bl(1:n),min(bu(1:n),x_0(:) ) ); 

Result.f_0 = nlp_f(x_0,Prob);
Result.x_0 = x_0;

% Lower bound on f(x)
f_Low = DefPar(Prob,'f_Low',-1E300);

LargeScale = DefPar(Prob,'LargeScale',1);

if LargeScale
   A = sparse(Prob.A);
   
   if isempty(Prob.ConsPattern)
      ConsPattern = sparse( ones(m2,n) );
   else
      ConsPattern = sparse( Prob.ConsPattern );
   end
   
   nz = nnz(A) + nnz(ConsPattern);
   
else
   A = full(Prob.A);
   ConsPattern = [];
   nz = n*m;
end

if m2 > 0
   % Determine the sparse problem structure
   if ~isempty(ConsPattern)
      [ix,iy]=find(ConsPattern);
      
      % Send linear index from multiple subscripts for nonzero pattern
      Prob.ConsIdx = sub2ind(size(ConsPattern),ix,iy);
   end
end

specs    = DefPar(KNITRO,'specs',[]); % Obsolete

if isfield(options,'ALG')
  alg = options.ALG;
elseif isfield(options,'alg')
  alg = options.alg;
else
  alg = 0;
end

switch(alg)
 case 1,
  S='Interior/Direct NLP KNITRO 4.0';
 case 2,
  S='Interior/CG NLP KNITRO 4.0';
 case 3,
  S='Active Set SLQP NLP KNITRO 4.0';
 otherwise
  S='Interior Point NLP KNITRO 4.0';
end
Result.SolverAlgorithm=S;

optParam = Prob.optParam;

% Max iters
if ~isfield(options,'MAXIT') & ~isfield(options,'maxit')
   options.MAXIT = DefPar(optParam,'MaxIter',2000);
end

% Max time
if ~isfield(options,'MAXTIME') & ~isfield(options,'maxtime')
  options.MAXTIME = DefPar(Prob,'MaxCPU',1e8);
end




% Problem type - LP or QP - but only if the user has not already set this.
if checkType('qp',Prob.probType)==1
   if ~isfield(options,'ISQP') & ~isfield(options,'isqp')
      options.ISQP = 1;
   end
end
if checkType('lp',Prob.probType)==1
   if ~isfield(options,'ISLP') & ~isfield(options,'islp')
      options.ISLP = 1;
   end
end

% KNITRO callback to m-file?
cb=DefPar(KNITRO,'callback','');
if ~isempty(cb)
   % SAL users - change this!
   switch(exist(cb))
      case {2,3,6}
         % m, mex, or p-file's are ok. Could check # of in/out too. 
      otherwise
         warning(['The callback function given (' callback ') is not an m-, mex-, or p-file. Callbacks will be disabled']);
         callback = '';
   end
end

% PrintFile
PrintFile = DefPar(KNITRO,'PrintFile',[]);

% Negative printlevel means we want printing to a file only. If no
% name is given, set knitro.txt as default name. 
if PriLev < 0 & isempty(PrintFile)
  PrintFile = 'knitro.txt';
end


% No Hessian d2L(x) information available? 
% Build a sparse structure, just a diagonal for LargeScale problems 
% and a full upper triangular 1's matrix for non-largescale.

if isempty(HL)
   if LargeScale
      HL = speye(n,n);
      % ioptions(8) = 6;
   else  
      HL = sparse(triu(ones(n,n))); 
   end
end

if ~issparse(HL), HL = sparse(HL); end
nHess = nnz(HL);

% Determine the sparse problem structure
if ~isempty(HL)
    % row and col indices. 
    [hessrow,hesscol]=find(HL);
    
    % Send linear index from multiple subscripts for nonzero pattern
    Prob.D2Lidx = sub2ind(size(HL),hessrow,hesscol);
end

[Inform,iter,x_k,f_k,g_k,c_k,v_k,cJac,H_k] = ...
   Tknitro(n,m1,m2,...
   xl,xu,x_0,...
   A,cl,cu,...
   ConsPattern,nz,...
   hesscol-1,hessrow-1,...
   f_Low,...
   PriLev,PrintFile,options,...
   cb,Prob);

% Tknitro version 4.0 returns H_k only for HESSOPT = 1;
if ~isempty(H_k)
   H_k = triu(H_k,1) + diag(diag(H_k)) + triu(H_k,1)';
end

Result.Iter = iter;
Result.x_k  = x_k;
Result.x_0  = x_0;
% Result.f_k  = nlp_f(x_k,Prob);
% Result.g_k  = nlp_g(x_k,Prob);
%Result.H_k  = nlp_H(x_k,Prob);
Result.f_k  = f_k;
Result.g_k  = g_k;
Result.H_k  = H_k;
Result.v_k  = [ v_k(m+1:end) ; v_k(1:m) ];
Result.cJac = cJac;

% Separate linear and nonlinear constraints
if m1>0, Result.Ax  = c_k(1:m1);   else Result.Ax = [];  end
if m2>0, Result.c_k = c_k(m1+1:m); else Result.c_k = []; end
    

% StateDef wants bl,bu to be [vars,nonlin,lin], unlike what is returned by
% defblbu above, therefore set additional input Order to 1 to reverse order.

% The barrier method will put x variables slightly outside bounds
FEASTOL = DefPar(options,'FEASTOL',1E-6);

Result = StateDef(Result, x_k(1:n), Result.Ax, Result.c_k, ...
                  FEASTOL, optParam.bTol, optParam.cTol, bl, bu, 1);

switch(Inform)
 case {0,-4,-5}
  ExitFlag = 0;
  
  switch(Inform)
   case 0
    ExitText = 'Optimal solution found';
   case -4
    ExitText = 'The current point cannot be improved. Optimality conditions not fulfilled';
   case -5
    ExitText = 'The current point cannot be improved. Optimality conditions close to fulfilled';
  end
  
 case -1 % iterlimit
  ExitFlag = 1;
  ExitText = 'Iteration limit reached';
 
 case -2 % infeasible
  ExitFlag = 4;
  ExitText = 'Problem may be infeasible';
 
 case -3 % unbounded
  ExitFlag = 2;
  ExitText = 'Problem seems to be unbounded';
 
 case -6 % time limit
  ExitFlag = 1;
  ExitText = 'Time limit reached';
  
 case -99
  ExitText = 'Insufficient memory to solve problem';
  ExitFlag = 10;
  
 case -98
  ExitText = 'Evaluation error (sqrt(neg) or div. by zero)';
  ExitFlag = 10;
  
 case -97
  ExitText = 'Unrecoverable error in LP solver';
  ExitFlag = 10;
  
 case -100
  ExitText = 'Termination requested in user routine';
  ExitFlag = 10; % ?
  
 otherwise
  ExitText = sprintf('Unknown status returned: %d',Inform);
  ExitFlag = 10;
end

% Values in the range -50 to -60 indicate input errorsless than -50 indicate 
if Inform >= -60 & Inform <= -50
   ExitFlag = 10;
   ExitText = 'Input error';
end


Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG
% 030505 ango Wrote file
% 030708 ango Hessian handling added
% 030729 ango Comments update
% 030829 ango Further comments updates
% 030903 ango ExitFlag/Inform handling
% 030917 ango f_Low implemented; LargeScale=1 is default.
% 031211 ango Changed the name of this file from TknitroTL
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040109 hkh  Skip call to ProbCheck, use returned f_k and g_k 
% 040203 ango LargeScale Hessian handling modified. 
% 040602 med  Help fix PriLev to PriLevOpt
% 041005 ango ISLP, ISQP checked and set for LP, QP problems. 
% 041110 ango Help updated for KNITRO 4, MaxCPU added
% 041111 med  Help fixes, MAXIT and MAXTIME fixed.
% 041128 hkh  Different calls to iniSolve dependent on HESSOPT
% 041129 med  optPar settings removed
% 041202 hkh  Removed inefficient call to Prob.optParam, use optParam instead
% 041202 hkh  Use new StateDef (and defblbu), avoid reversing bl/bu order
% 041202 hkh  Use FEASTOL instead of xTol, interior point methods are inexact
% 041210 hkh  Revise GRADOPT,HESSOPT use, avoid two level of differentiation
