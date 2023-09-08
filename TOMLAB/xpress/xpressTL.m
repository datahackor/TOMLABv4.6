% TOMLAB/ Xpress LP, QP, MILP and MIQP Solver
%
% xpressTL converts the problem from the Tomlab structure format and 
% calls xpress.m. On return converts the result to the Tomlab structure format.
% Also see the help for xpress.m
% 
% xpressTL.m solves the following mixed integer (linear or quadratic) 
% programming problem (LP, QP, MILP, MIQP):
%
%   minimize   0.5 * x'*F*x + c'x     subject to:	
%      x             x_L <=    x   <= x_U
%                    b_L <=   Ax   <= b_U
%   where 
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the Xpress-MP sparse matrix format.
%   c, x_L, x_U has dimension n
%   b_L, b_U has dimension m
%   F is a n x n symmetric matrix, sparse or dense. 
%   If F is empty, an LP or MILP problem is solved
%   Some or all x may be integer valued as specified by other input variables
%
% function Result = xpressTL(Prob)
%
% INPUT:  
% Prob   Problem structure in TOMLAB format
%
%
% Fields used in input structure Prob 
% Use lpAssign, mipAssign or qpAssign to define the Prob structure
%
%
% x_L, x_U  Lower and upper bounds on variables. 
% b_L, b_U  Lower and upper bounds on linear constraints. 
% A         Linear constraint matrix, dense or sparse m x n matrix
% QP.c      Linear coefficients in objective function, size n x 1
% QP.F      Quadratic matrix of size n x n 
% PriLevOpt Print level in xpressTL, the xpress m-file and xpressmp C-interface.
%           = 0  Silent
%           = 1  Warnings and Errors
%           = 2  Summary information
%           = 3  More detailed information
%           > 10 Pause statements, and maximal printing (debug mode)
%
% -----------------------------------------------
% Fields used in Prob.XPRESS: 
% -----------------------------------------------
%
% LogFile   File to write Xpress-MP log output to. Default is empty '' in which 
%           case nothing is written. Please note that Xpress-MP appends
%           it's output to the log file. 
%
% SaveFile  Filename for writing the problem prior to calling the Xpress-MP
%           solver. If empty, no file is written. The type of output is
%           determined by the SaveMode parameter. 
%           Xpress-MP will always add an extension to the filename given here.
%           The extension depends on the SaveMode chosen, see below.
%
% SaveMode  Character string with any combination of the following
%           character flags: 
%
%           p  full precision of numerical values;
%           o  one element per line
%           n  scaled
%           s  scrambled vector names
%           l  output in LP format
%
%           The extension added to the SaveFile name is .mat, 
%           unless the 'l' flag is used in which case the extension is .lp.
%
% iis        Flag indicating whether to compute an IIS and return it to
%            MATLAB. This option can only be set for an LP problem. If an
%            IIS is found, XPRESS automatically changes the problem to make
%            it feasible and reoptimizes it.
%               
%            = 0     Don't return IIS to MATLAB (default).
%            = 1     Compute IIS and return it to MATLAB if an LP problem 
%                    has been proven infeasible.
%
%            The IIS is returned through the output parameter 'iis'.
%
% iisFile    Flag indicating whether to write a file describing the IIS set
%            or not. If is set to 1, a file: LPprob.iis will be written.
%            Otherwise, no file is written.
%
% sa         Structure telling whether and how you want XPRESS to perform a
%            sensitivity analysis (SA). You can complete an SA on the 
%            objective function and right hand side vector. The saRequest
%            structure contains two sub structures:
%
%                .obj and .rhs
%               
%            They have one field each:
%
%                .index
%
%            In case of .obj.index, .index contains the indices of the columns 
%            whose objective function coefficients sensitivity ranges are 
%            required.
%
%            In case of .rhs.index, .index contains the indices of the rows 
%            whose RHS coefficients sensitivity ranges are required.
%
%            In both cases, the .index array has to be sorted, ascending.
%
%            To get an SA of objective function on the four variables 120 
%            to 123 (included) and variable 6 the saRequest structure would
%            look like this:
%
%                saRequest.obj.index = [6 120 121 122 123];
%
%            The result is returned through the output parameter 'sa'.
%
% -----------------------------------------------
% Fields used in Prob.MIP:
% -----------------------------------------------
% See the corresponding variables in xpress.m for an explanation
%
% MIP.IntVars
%           Defines which variables are integers, of general type I or binary B
%           Variable indices should be in the range [1,...,n].
%           IntVars is a single integer ==> Variable 1:IntVars are integer 
%           IntVars is a logical vector ==> x(find(IntVars > 0)) are integers 
%           IntVars is a vector of indices ==> x(IntVars) are integers 
%           (if [], then no integers of type I or B are defined)
%           xpress checks which variables has x_L=0 and x_U=1, i.e. binary.
%
% MIP.PI
%           Integer variables of type Partially Integer (PI), i.e. takes an
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
% MIP.SC
%           A vector with indices for the integer variables of type 
%           Semi-continuous (SC), i.e. that takes either the value 0 or a 
%           real value in the range [x_L(i),x_U(i)], assuming for some j,
%           i = SC(j), where i is an variable number in the range [1,n].
%
% MIP.SI
%           A vector with indices for the integer variables of type 
%           Semi-integer (SI), i.e. that takes either the value 0 or 
%           an integer value in the range [x_L(i),x_U(i)], assuming for some j,
%           i = SI(j), where i is an variable number in the range [1,n].
%
% MIP.sos1
%           A structure defining the Special Ordered Sets of Type One (sos1). 
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
% MIP.sos2 
%           A structure defining the Special Ordered Sets of Type Two (sos2). 
%           Specified exactly as sos1 sets, see sos1 input variable description
%
% MIP.xpcontrol
%           Structure, where the fields are set to the Xpress-MP control
%           parameters (Xpress-Optimizer Reference Manual Section 7) that the
%           user wants to specify values for. The prefix XPRS_ is not used. 
%           Example with 1 int, 1 double and 1 character valued: 
%
%  xpcontrol.LPITERLIMIT   = 50;   Setting max number of global iterations
%  xpcontrol.OPTIMALITYTOL = 1E-5; Changing reduced cost tolerance 
%  xpcontrol.OBJNAME ='ObjF1'; New name of objective function( <= 255 chars)
%
% callback  Logical vector defining which callbacks to use in Xpress-MP
%  If the i:th entry of logical vector callback is set, the corresponding 
%  callback is defined. See Section 5.3 in Xpress-Optimizer Reference Manual.
%  The callback calls the m-file specified below. The user may edit this file,
%  or make a new copy, which is put before in the Matlab path.
%  See xpress.m for a complete list of all callbacks available.
%
% -----------------------------------------------
% Fields used in Prob.optParam: (Structure with optimization parameters)
% -----------------------------------------------
% MaxIter   Limit of iterations  (if not MIP.xpcontrol.LPITERLIMIT is set)
%
% -----------------------------------------------
% OUTPUT: 
% -----------------------------------------------
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum
% x_k      Solution vector
% x_0      Initial  solution vector not known, set as empty
% g_k      Exact gradient computed at optimum, computed as c or c + Fx
% xState   State of variables.Free==0; On lower == 1; On upper == 2; Fixed == 3;
% bState   State of constraints. Free==0; Lower == 1; Upper == 2; Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector)
%          v_k = [rc;v]. rc n-vector of reduced costs. v holds m dual variables
% ExitFlag - exit status from xpress.m (similar to Tomlab)
% Inform   Xpress-MP information parameter. 
%           for LP:  xpProblemAttrib.LPSTATUS;
%           for MIP: xpProblemAttrib.MIPSTATUS;
%           0 = Optimal solution found
%           2 = Unbounded solution
%           4 = Infeasible problem
%           5 = Some error occured
% rc       Reduced costs. If ninf=0, last m == -v_k
% Iter     Number of iterations / nodes visited
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
%           FuncEv and ConstrEv set to Iter. GradEv=0.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in Tomlab QP standard 
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver  (Xpress-MP)
% SolverAlgorithm  Description of the solver
%
% -----------------------------------------------
% Output fields in Result.MIP:
% -----------------------------------------------
% MIP.slack     Slack variables (m x 1 vector)
% MIP.ninf      Number of infeasibilities
% MIP.sinf      Sum of infeasibilities
% MIP.lpiter    Number of LP iterations
% MIP.glnodes   Number of nodes visited
% MIP.basis     Basis status of constraints + variables, (m + n x 1 vector)
%               in the Xpress-MP format, fields xState and bState has the same
%               information in the Tomlab format.
%
% -----------------------------------------------
% Output fields in Result.XPRESS:
% -----------------------------------------------
% iis       Structure containing IIS information (niis x 1). niis is the
%           number of IISs found.
%           Fields:
%
%               iisStatus   Status flag. (Only set in the first element of
%                           the iis array.) Possible values:
%
%                           2  IIS was written to file LPprob.iis.
%                           1  IIS was obtained.
%                           -1  Problem was infeasible but no IIS found.
%                           -2  Problem was not infeasible.
%
%               iisMessage  Error message on error. (Only set in the first 
%                           element of the iis array.)
%               colind      The column indices of the IIS set.
%               rowind      The row indices of the IIS set.
%
% sa        Structure with information about the requested SA, if requested.
%           The fields:
%
%               obj         Ranges for the variables in the objective function.
%
%               rhs         Ranges for the right hand side values.
%
%           These fields are structures themselves. Both structures 
%           have identical field names:
%
%               status      Status of the SA operation. Possible values:
%
%                            1  Successful
%                            0  SA not requested.
%                           -1  Error: MIP problem was presolved.
%
%               lower       The lower range.
%
%               upper       The upper range.
%
%
% Also set into the Result.MIP output structure is:
%
% MIP.xpControlVariables  All Xpress-MP control variables
% MIP.xpProblemAttrib     All Xpress-MP problem attributes
%
% They are also available as global variables after the run:
% global xpProblemAttrib xpControlVariables
%
% -----------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 2004.0.0 $
% Written July 8, 1999.      Last modified Jan 17, 2005.

function Result = xpressTL(Prob)

%#function lp_f lp_g lp_H xp2control xp2problem

if nargin < 1 
    error('xpressTL needs the Prob structure as input');
end

global MAX_x MAX_c % Max number of variables and constraints to print

global xpProblemAttrib xpControlVariables

Prob.solvType = 11; % MIQP solver

Prob = iniSolve(Prob,11,0,0);

%DEBUG=0;

PriLev=Prob.PriLevOpt;
%if DEBUG
%   % Now for debugging
%   PriLev=1000
%end

Result=ResultDef(Prob);
Result.Solver='Xpress-MP';
Result.SolverAlgorithm='Xpress-MP LP/QP/MIP/MIQP Solver';

% Initial checks on the inputs

%
% Define lower and upper bound arrays for Xpress
%
% Inf are changed to BIG (=1E10), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%

% The xpress-mp MEX has 1E20 as default, but is changed if BIG is noempty
% Note that xpress-mp is changing all values >= BIG to XPRS_PLUSINFINITY
% which is hard coded to 1E20
% Similar all values <= -BIG are set to XPRS_MINUSINFINITY (hard coded 1E20)
Prob.BIG=DefPar(Prob,'BIG',1E20);

[bl, bu, n, m, m2] = defblbu(Prob, Prob.BIG);

nTot=n+m;

% Initial values (xpress does not use x_0)

Fzero = isempty(Prob.QP.F);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m, error('Linear constraints A MUST have m rows!'); end 
end 

Result.f_0=0;

% Check if any linear part
c = Prob.QP.c(:);

if isempty(c), c=zeros(n,1); end

IntVars = DefPar(Prob.MIP,'IntVars');
PI      = DefPar(Prob.MIP,'PI');
SC      = DefPar(Prob.MIP,'SC');
SI      = DefPar(Prob.MIP,'SI');
sos1    = DefPar(Prob.MIP,'sos1');
sos2    = DefPar(Prob.MIP,'sos2');
    
if isempty(IntVars) & isempty(sos1) & isempty(sos2) & isempty(SC) & isempty(SI) & isempty(PI)
    MIP=0;
else
    MIP=1;
end

xpcontrol = DefPar(Prob.MIP,'xpcontrol',[]);
callback  = DefPar(Prob.MIP,'callback',zeros(15,1));

if length(callback) < 15
   callback=[callback(:);zeros(15-length(callback),1)];
end

XPRESS = DefPar(Prob,'XPRESS',[]);
SaveFile = DefPar(XPRESS,'SaveFile','');
SaveMode = DefPar(XPRESS,'SaveMode','');
LogFile  = DefPar(XPRESS,'LogFile','');

% Timing
MaxCPU = DefPar(Prob,'MaxCPU',inf);
if isfinite(MaxCPU) & ~isfield(xpcontrol,'MAXTIME')
   xpcontrol.MAXTIME = MaxCPU;
end

if MIP
   PI      = DefPar(Prob.MIP,'PI',[]);
   SC      = DefPar(Prob.MIP,'SC');
   SI      = DefPar(Prob.MIP,'SI');
   sos1    = DefPar(Prob.MIP,'sos1');
   sos2    = DefPar(Prob.MIP,'sos2');

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
   %   callback(8)=1;
   %   callback(9)=0;
   %   callback(1:7)=0;
   %   callback(10)=1;
   %   callback(11)=1;
   %end
end

%if ~isfield(xpcontrol,'PRESOLVE')
%   % Default use no preSolve
%   xpcontrol.PRESOLVE=0;    % No presolve
%end
if ~isfield(Prob.optParam,'MaxIter')
   Prob.optParam.MaxIter=2000;
end
if ~isfield(xpcontrol,'LPITERLIMIT') & Prob.optParam.MaxIter~=2000
   xpcontrol.LPITERLIMIT=Prob.optParam.MaxIter;
end

%if DEBUG
%   xpcontrol.OPTIMALITYTOL=1E-6;
%
%   % Test of Illegal field
%   xpcontrol.URK=[];
%
%   % Test of short string < 256 characters
%   xpcontrol.MPSOBJNAME='Obj1';
%
%   %xpcontrol.DEFAULTALG=1;
%   xpcontrol
%
%   xpcontrol.MIPLOG=3;       % Global log callback at every node
%   xpcontrol.LPLOG=1;        % Simplex log callback at every node
%   xpcontrol.CUTSTRATEGY=2;  % Aggressive cut strategy
%   xpcontrol.CUTSTRATEGY=1;  % Conservative cut strategy
%   xpcontrol.CUTSTRATEGY=0;  % No cut strategy
%end

Prob.XPRESS = DefPar(Prob, 'XPRESS',[]);
iisRequest  = DefPar(Prob.XPRESS, 'iis', []);
iisFile     = DefPar(Prob.XPRESS, 'iisFile', '');
saRequest   = DefPar(Prob.XPRESS, 'sa', []);

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('xpress-mp',2,...
%         nnObj, 0, size(Prob.A,1), Prob);

[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, iis, sa] = ...
  xpress(c, Prob.A, bl(1:n), bu(1:n), bl(n+1:nTot), bu(n+1:nTot), xpcontrol,...
  callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2, sparse(Prob.QP.F),...
  LogFile, SaveFile, SaveMode, iisRequest, iisFile, saRequest);


%xprint(x,'x:');
%xprint(slack,'s:');
%xprint(v,'v:');
%xprint(rc,'rc:');
%xprint(basis(1:m),'cbasis:');
%xprint(basis(m+1:m+n),'xbasis:');

if isempty(xpProblemAttrib)
  % If something goes really wrong in the mex, this might happen. 
  lpstatus = 0;
  glstatus = 0;
else
  lpstatus=xpProblemAttrib.LPSTATUS; 
  glstatus=xpProblemAttrib.MIPSTATUS;
end
 
if MIP
   switch glstatus
     case 0
       Text = 'No problem has been loaded';
     case 1
       Text = 'LP/QP has not been optimized';
     case 2
       Text = 'LP/QP has been optimized';
     case 3
       Text = 'Global search incomplete - no integer solution found';
     case 4
       Text = 'Global search incomplete - integer solution found';
     case 5
       Text = 'Global search complete - No integer solution found ';
     case 6
       Text = 'Global search complete - integer solution found';
     otherwise
       Text = sprintf('Unknown glstatus %d returned',glstatus);
   end

else
   switch lpstatus
     case 1
         Text = 'Optimal solution found';
     case 2
         Text = 'Infeasible';
     case 3
         Text = 'Objective worse than cutoff';
     case 4
         Text = 'Unfinished';
     case 5
         Text = 'Unbounded';
     case 6
         Text = 'Cutoff in dual';
     case 7
         Text = 'Unsolved';
    otherwise
       Text = sprintf('Unknown lpstatus %d returned',lpstatus);
   end
end

% Special cases for some fatal errors, where there's no solution at
% all: 
if Inform < 0 
  switch(Inform)
    
   case -8,
    Text = ['Licensing error. Please check your license manager and' ...
	    ' xpress.lic file.' ];
   otherwise,
    Text = [ 'Unknown fatal error ' num2str(Inform)] ;
  end
end  

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / Xpress solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');

   if MIP
      fprintf('MIPSTATUS = %d. ',glstatus);

   else
      fprintf('LPSTATUS  = %d. ',lpstatus);
   end
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

Result.ExitText = Text;

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
      Result.g_k = Prob.QP.F*x + c;
      % z          = Prob.QP.F*x;
      % Result.g_k = z+c;
      % SICK. Xpress does only return linear part objective
      % Result.f_k = Result.f_k + 0.5*x'*z;
   else
      Result.g_k=c;
   end
elseif Fzero
   Result.g_k=[];
else
   Result.g_k = Prob.QP.F*x;
   % z          = Prob.QP.F*x;
   % Result.g_k = z;
   % SICK. Xpress does only return linear part objective
   % Result.f_k = Result.f_k + 0.5*x'*z;
end

Result.c_k=[];
Result.cJac=[];

% State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
if isempty(basis)
   Result.xState = [];
   Result.bState = [];
   Result.QP.B   = [];
else
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
   B       = xbasis;
   B(B==2) = -1;
   Result.QP.B=B;
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
Result.ExitFlag = Inform;
if MIP
   Result.Inform   = glstatus;
else
   Result.Inform   = lpstatus;
end

%xTol=Prob.optParam.xTol;


Result.MIP.ninf=ninf;
Result.MIP.sinf=sinf;
Result.MIP.slack=slack;
Result.MIP.lpiter=lpiter;
Result.MIP.glnodes=glnodes;
Result.MIP.basis=basis;
Result.MIP.xpControlVariables=xpControlVariables;
Result.MIP.xpProblemAttrib=xpProblemAttrib;

Result.XPRESS.iis = iis;
Result.XPRESS.sa = sa;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 991026 hkh  Making xpressTL.m separate from xpress.m
% 010710 hkh  Revision for v3.0
% 010713 hkh  Restrict BIG - infinity, to 1E10
% 010726 hkh  Add output in Result.MIP, sending all output from xpress.m back
% 011210 hkh  Revision for Rel13.
% 011213 hkh  Changed SC into SC and SI, semi-continuous and semi-integer
% 011226 hkh  Improved comments
% 020618 hkh  Improve speed for big Prob.QP.F matrices, avoid copying, testing
% 020702 hkh  For QP, basis empty, avoid using basis then
% 020810 hkh  Remove some wrong comments
% 030813 ango Better handling of Prob.MIP
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040525 ango Add SaveFile and SaveMode parameters
% 040528 hkh  Avoid lowering LPITERLIMIT, unless user has set optParam.MaxIter 
% 040528 hkh  Comment out all DEBUG parts
% 040803 med  Added pragmas for MATLAB Compiler
% 040830 fhe  Added support for iis and sa.
% 041213 hkh  Use BIG as Prob.BIG is nonempty, otherwise 1E20
% 041213 hkh  Report glnodes as Result.Iter for MIP problems
% 041213 hkh  Use lpiter for Result.FuncEv, GradEv, ConstrEv 
% 041213 hkh  Must use Prob.BIG, not BIG, in call to defblbu
