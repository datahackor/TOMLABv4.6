% TOMLAB /XA LP, MILP, and QP Solver
%
% xaTL converts the problem from the TOMLAB structure format and 
% calls XA.m. On return converts the result to the TOMLAB structure format.
% Also see the help for XA.m
% 
% XA solves the following types of problems:
%
% linear mixed integer (MILP)
% linear or quadratic continuous (LP,QP)
%
%   minimize   0.5 * x'*F*x + c'x    
%      x            
% 
%   subject to  x_L <=    x   <= x_U
%               b_L <=   Ax   <= b_U
%
%   where 
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the XA sparse matrix format.
%   c, x_L, x_U has dimension n x 1
%   b_L, b_U has dimension m x 1
%   F is a n x n symmetric matrix, sparse or dense. 
%   If F is empty, an LP or MILP problem is solved
%
% ---------------------------------------------------------------------------
%
% function Result = xaTL(Prob)
%
% INPUT:    
%
% Prob        Problem structure in TOMLAB format.
%             Use lpAssign, mipAssign or qpAssign to 
%             define the Prob structure.
%
%             Fields used in input structure Prob: 
%
%
%  x_L, x_U   Lower and upper bounds on variables, size n x 1
%  b_L, b_U   Lower and upper bounds on linear constraints, size m x 1 
%  A          Linear constraint matrix, dense or sparse m x n matrix
%   
%             NOTE - all bounds vectors - if [], +/- Inf is assumed
%   
%  QP.c       Linear objective function coefficients, size n x 1
%  QP.F       Quadratic matrix of size n x n 
%   
%  PriLevOpt  Print level in XATL, the XA m-file and XAmex C-interface.
%             = 0  Silent
%             = 1  Warnings and Errors
%             = 2  Summary information
%             = 3  More detailed information
%             > 10 Pause statements, and maximal printing (debug mode)
%
%
% Fields used in Prob.XA (Structure with XA specific parameters)
%
%  LogFile    Name of file to receive the XA iteration and results log. 
%             If empty or not present, no log is written. 
%
%  SaveFile   Name of file for saving the problem in MPS format. The
%             extension ".mps" is added automatically. 
%             If empty or not present, the save feature is disabled. 
%
%  callback   0/1 vector of length 10. Defines the active callbacks during
%             the solving process. 
%
%             The callback routines and their corresponding indices in the 
%             callback vector are: 
%
%              1 - xacb_begin.m   - before solving
%              2 - xacb_infeas.m  - infeasible iteration
%              3 - xacb_feas.m    - feasible iteration
%              4 - xacb_node.m    - integer node generated
%              5 - xacb_intsol.m  - integer solution found
%              6 - xacb_branch.m  - user selects b & b variable
%              7 - xacb_barrier.m - barrier iteration  
%              8 - xacb_resolve.m - problem solve -> fast modify resolve
%              9 - N/A            - currently not used
%             10 - xacb_end.m     - after solving 
%
%             The user can either modify the existing xacb_*.m
%             files (use i.e. "which xacb_feas") to find them,
%             or copies could be made and placed before the original
%             files in the MATLAB path. 
%  
%             The calling syntax for all XA callbacks is
%  
%                 function xacb_*( xacbInfo, Prob )
%  
%             and is described in more detail in xacb.m (help xacb.m)
%
%  iis        A flag telling whether to obtain an IIS when a
%             problem is determined infeasible. The flag can be set
%             to one of the following values:
%
%              0 - Do not search for an IIS (default).
%              1 - Implementation of irreducible inconsistent systems
%                  (IIS) of constraints. Algorithm: IIS.
%              2 - Method of locationing a minimal number of constraints
%                  such that if all are removed the model is
%                  feasible. Algorithm: Block.
%
%             The IIS is returned to the field: Result.XA.iis
%             Information about the IIS is automatically written to
%             the LogFile if a LogFile is defined.
%
% Fields used in Prob.MIP
%
%  xaControl  Struct array with fields for XA parameters. For example: 
%
%             xaControl.Barrier = 'Yes';   (also 1/0 for Yes/No is accepted)
%             xaControl.ITERATIONS = 1000; 
% 
%  IntVars    Defines which variables are integers, of general type I or binary B
%             Variable indices should be in the range [1,...,n].
%             IntVars is a single integer ==> Variable 1:IntVars are integer 
%             IntVars is a 0-1 vector ==> x(find(IntVars > 0)) are integers 
%             IntVars is a vector of indices ==> x(IntVars) are integers 
%             (if [], then no integers variables are defined)
%             XA checks which variables has x_L=0 and x_U=1, i.e. binary.
%
%  VarWeight  Weight for each variable in the variable selection phase.
%             A lower value gives higher priority. XA uses integer values,
%             ideally 1 or higher. 
%   
%  SC         A vector with indices for the Type 1 Semi-Continous variables,
%             i.e. that takes either the value 0 or a value in the range 
%             [ x_L(i) , x_U(i) ].
%
%  SC2        A vector with indices for the Type 2 Semi-Continous variables,
%             i.e. that takes either the value of the corresponding upper bound,
%             or a value in the range [ x_L(i) , 0 ].
%
% Fields used in Prob.optParam: (Structure with optimization parameters)
% 
%   MaxIter   Limit of iterations (if not XA.xaControl.ITERATION is set)
% ------------------------------------------------------------------------------
%
% OUTPUTS: 
%
% Result   Structure with results (see ResultDef.m):
%
%  f_k        Function value at optimum
%  x_k        Solution vector
%  x_0        Initial  solution vector not known, set as empty
%  g_k        Exact gradient computed at optimum, computed as c or c + Fx
%
%  xState     State of variables.   Free==0; On lower == 1; On upper == 2; Fixed == 3;
%  bState     State of constraints. Free==0; On lower == 1; On upper == 2; Equality == 3;
%
%  v_k        Lagrangian multipliers (for bounds + dual solution vector)
%             v_k = [rc;v]. rc n-vector of reduced costs. v holds m dual variables
%
%  rc         Reduced costs. If ninf=0, last m == -v_k
%
%  ExitFlag   Exit status, TOMLAB standard
%
%  Iter       Number of iterations / nodes visited
%
%  FuncEv     Number of function evaluations. Set to Iter.
%
%  GradEv     Number of gradient evaluations. Set to Iter if
%             QP/MIQP, otherwise 0. 
%             FuncEv and ConstrEv set to Iter. GradEv=0.
%
%  ConstrEv   Number of constraint evaluations. Set to 0.
%
%  QP.B       Basis vector in TOMLAB QP standard.
%
%  Solver           Name of the solver  (XA)
%  SolverAlgorithm  Description of the solver
%
% Fields used in Result.XA:
%
%  iis        Structure containing IIS information. 
%             Fields:
% 
%                iisStatus   Status flag. Possible values:
% 
%                             1 - IIS was obtained.
%                             0 - IIS was not requested.
%                            -1 - Problem is infeasible but no IIS found.
%                            -2 - Problem is not infeasible.
% 
%                iisMessage  Status message.
%                rowind      The row indices of the IIS set.
% 
% -----------------------------------------------------------------------
%
% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 14.0 $
% Written Oct 16, 2003.      Last modified Jan 17, 2005.

function Result = xaTL(Prob)

if nargin < 1 
    error('xaTL needs the Prob structure as input');
end

global MAX_x MAX_c % Max number of variables and constraints to print

% These fields now undefined
% -----------------------------------------------
% Output fields in Result.MIP:
% -----------------------------------------------
%   MIP.slack     Slack variables (m x 1 vector)
%   MIP.ninf      Number of infeasibilities
%   MIP.sinf      Sum of infeasibilities
%   MIP.lpiter    Number of LP iterations
%   MIP.glnodes   Number of nodes visited
%   MIP.basis     basis status of constraints + variables, (m + n x 1 vector)
%                 in the XA format, fields xState and bState has the same
%                 information in the Tomlab format.
%

MIP = DefPar(Prob,'MIP',[]);
IntVars   = DefPar(MIP,'IntVars',[]);

if isempty(IntVars)
   Prob.solvType = 2; % QP solver
else
   Prob.solvType = 8; % MILP solver
end

Prob = iniSolve(Prob,Prob.solvType,0,0);

PriLev=Prob.PriLevOpt;

Result=ResultDef(Prob);
Result.Solver='XA';
Result.SolverAlgorithm='XA LP/MILP/QP solver';

% Initial checks on the inputs

%
% Define lower and upper bound arrays for XA
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%

BIG=2E23;

[bl, bu, n, m, m2] = defblbu(Prob, BIG, 1);

% Check if any linear part
c = DefPar(Prob.QP,'c',zeros(n,1));

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   %   fprintf('mA = %i, m = %i\n', mA, m);
   if mA~=m, error('Linear constraints A MUST have m rows!'); end 
end 

Result.f_0=0;

VarWeight = DefPar(MIP,'VarWeight',[]);
SC1       = DefPar(MIP,'SC',[]);
SC2       = DefPar(MIP,'SC2',[]);

if(isempty(IntVars) & isempty(VarWeight) & isempty(SC1) & isempty(SC2))
  mip = 0;
else
  mip = 1;
end

xaControl = DefPar(MIP,'xaControl',[]);
xaclp     = DefPar(MIP,'xaCtrl',{});

if ~isfield(Prob.optParam,'MaxIter')
   Prob.optParam.MaxIter=2000;
end
%if ~isfield(xaControl,'ITERATION') & Prob.optParam.MaxIter~=2000
if ~isfield(xaControl,'ITERATION')
   xaControl.ITERATION=Prob.optParam.MaxIter;
end

MaxCPU = DefPar(Prob, 'MaxCPU', 2e9);
if ~isfield(xaControl,'TIMELIMIT')
    xaControl.TIMELIMIT = MaxCPU;
end

% Unnecessary. tom2xaclp does this.
% if isfield(xaControl,'ITERATION')
%    % Iteration iteration limit
%    xaclp{end+1} = sprintf('Set Iteration %d',xaControl.ITERATION);
% end

XA = DefPar(Prob,'XA');
LogFile    = DefPar(XA,'LogFile','');
SaveFile   = DefPar(XA,'SaveFile','');
callback   = DefPar(XA,'callback',zeros(1,10));
iisRequest = DefPar(XA,'iis', 0);

xaclp = tom2xaclp(xaControl,xaclp);

if PriLev>4
   xaControl
   celldisp(xaclp)
end

[x_k,f_k,Inform,modsts,solsts,act,dact,status,iis] = xa(c,Prob.A,bl(1:n),bu(1:n), ...
   bl(n+1:n+m),bu(n+1:n+m),xaclp,callback,PriLev,LogFile,SaveFile,Prob,...
   IntVars,VarWeight,SC1,SC2,Prob.QP.F,iisRequest);

Result.x_k = x_k;
Result.f_k = f_k;
Result.x_0=[];
Result.v_k=act+dact;

%Result.g_k=[];

Result.Inform   = Inform;

Result.XA.rc     = Inform;
Result.XA.modsts = modsts;
Result.XA.solsts = solsts;
Result.XA.act    = act;
Result.XA.dact   = dact;
Result.XA.status = status;
Result.XA.iis    = iis;

switch(modsts)
 case 1  , if(mip), ModText = 'Optimal Integer Solution'; else ModText='Optimal Solution'; end
 case 2  , ModText = 'Integer Solution (not proven the optimal integer solution)';
 case 3  , ModText = 'Unbounded solution';
 case 4  , ModText = 'Infeasible solution';
 case 5  , ModText = 'Callback function indicates Infeasible solution';
 case 6  , ModText = 'Intermediate infeasible solution';
 case 7  , ModText = 'Intermediate nonoptimal solution';
 case 9  , ModText = 'Intermediate Non-integer solution';
 case 10 , ModText = 'Integer Infeasible';
 case 13 , ModText = 'More memory required to load/solve model. Increase memory request in XAINIT call';
 case 32 , ModText = 'Integer branch and bound process currently active, model has not completed solving';
 case 99 , ModText = 'Currently solving model, model has not completed solving';
  
 otherwise, ModText = 'Unknown XA Model Status';
end

% Again.. for ExitFlag
switch(modsts) 
 case {1,2}, ExitFlag = 0;  % OK
 case 3    , ExitFlag = 2;  % Unbounded
 case {4,5,6,7,9,10}, ExitFlag = 4; % Infeasible
 case 13   , ExitFlag = 11; % Memory errors
 otherwise, ExitFlag = -1;
end
  
% Note: Integer problems return a Model Status code of 1 if the optimal
% integer solution is found. If XA has not proven that its integer
% solution is optimal, then a Model Status code of 2 is returned.

switch(solsts)
 case 1  , SolText = 'Normal Completion';
 case 2  , SolText = 'Iteration Interrupt';
 case 3  , SolText = 'Resource Interrupt, like time limit exceeded';
 case 4  , SolText = 'Terminated by User, probably a Ctrl+Z';
 case 8  , SolText = 'Node Table Overflow';
 case 10 , SolText = 'Solver Failure';
  
 otherwise, SolText = sprintf('Unknown XA Solver Status %d',solsts);
end


% Some solsts values should change the value of ExitFlag
switch(solsts)
  case {2,3}, ExitFlag = 1; 
end

Result.ExitFlag = ExitFlag;
ExitText = [ SolText ' : ' ModText ];
Result.ExitText = ExitText;

switch(Inform)
   case 2  , RCText='XAINIT XACLP XAMPSI Successful.';
   case 4  , RCText='XARCCI/M/R XADBFI/M XALOAD Successful.';
   case 6  , RCText='XASOLV Successful.';
   case 8  , RCText='XAUNDO Successful.';
   case 10 , RCText='XADONE Successful.';
   case 101, RCText='Allocation amount in XAINIT must be greater than or equal to 0. Make sure argument is a positive long_int integer.';
   case 102, RCText='Allocation amount in XAINIT must be greater than or equal to 0. Make sure argument is a positive long_int integer.';
   case 103, RCText='Memory allocation error, no memory available.';
   case 104, RCText='Memory allocation error, requested memory not available.';
   case 110, RCText='XACLP called out of sequence.';
   case 111, RCText='XADONE called out of sequence.';
   case 112, RCText='XASOLV called out of sequence.';
   case 113, RCText='XAUNDO called out of sequence.';
   case 114, RCText='XAACT,XADUAL,XAACTC,XADUALC,XABIA,XACBIA XAVAR, XAANAL, XAANALA called out of sequence.';
   case 115, RCText='XALOAD called out of sequence.';
   case 116, RCText='XAMPSI, XADBFI, XARCCI called out of sequence.';
   case 117, RCText='XAOK called out of sequence.';
   case 118, RCText='XARPRT called out of sequence.';
   case 122, RCText='BXA.DLL not in path. Either move BXA.DLL to a directory in your path statement, or add the directory BXA.DLL is located in to your path statement in your autoexec.bat. You will have to restart Windows and possibly your computer to correct your PATH command.';
   case 124, RCText='Wrong version of BXA.DLL or your version has been damaged.';
   case 126, RCText='XA.... routine called after XADONE.';
   case 128, RCText='Restart Windows, you are having programming difficulties and it has caught up with you.';
   case 132, RCText='In call to XART1, the XAInfo structure (XAInfoSize value) size is incorrect. Check to see that dataarea.XAInfoSize = sizeof(XAInfo).';
   case 134, RCText='In call to XART2, the XAMsgLine structure (XAMsgSize value) size is incorrect. Check to see that dataarea.XAMsgSize = sizeof(XAMsgLine).';
   case 136, RCText='In call to XAMODEL, the XAModel structure (XAModel Size value) size is incorrect. Check to see that dataarea.XAModel Size = sizeof(XAModel).';
   case 204, RCText='XACLP is passed a command line parameter which expects No or Yes as a value.';
   case 205, RCText='XACLP is passed a misspelled command line parameter.';
   case 207, RCText='XACLP is passed an unreasonable command line parameter value.';
   case 208, RCText='XACLP is passed a PAGESIZE, TMARGIN, or BMARGIN command which is incompatible with each others value.';
   case 209, RCText='XACLP is passed a LINESIZE or LMARGIN command which is incompatible with each others value.';
   case 210, RCText='XACLP is passed a command line which ends prematurely.';
   case 211, RCText='XACLP is unable to process the OUTPUT command line parameter value.';
   case 214, RCText='XACLP is passed a FIELDSIZE or DECIMALS command which is incompatible with each others value.';
   case 300, RCText='Column number out of range, or misspelled row or column name.';
   case 301, RCText='Too many or zero rows in problem, increase maxrow.';
   case 302, RCText='Too many or zero columns in problem, increase maxcol.';
   case 303, RCText='Free memory error, XA could not release it''s memory. Probably XA data storage area has been clobbered.';
   case 304, RCText='Column lower bound > upper bound.';
   case 305, RCText='Row lower bound > upper bound.';
   case 306, RCText='Problem too large for available memory, increase use parameter in your call to XAINIT.';
   case 307, RCText='Too many nonzeros, increase maxnonzero.';
   case 308, RCText='Row number out of range.';
   case 309, RCText='Duplicate column and/or row name.';
   case 310, RCText='Arrays overlap. Run XAOK to identify overlapping arrays.';
   case 312, RCText='Bad colptr values. Run XAOK to locate problem.';
   case 314, RCText='Unreasonable rownos values. Run XA to locate problem.';
   case 316, RCText='Duplicate rownos value for the same column.';
   case 318, RCText='Unreasonable technology coefficient in a array.';
   case 320, RCText='Unreasonable value in status array.';
   case 322, RCText='Unreasonable value in priority array, check XASETA call.';
   case 324, RCText='Unreasonable value in increment array, check XASETA call.';
   case 326, RCText='Semi-continuous type 1 or type 2 column missing lower bound.';
   case 328, RCText='Column can not have both semi-continuous type 1 and type 2 at the same time.';
   case 330, RCText='Unrecognized line in MPS file.';
   case 400, RCText='MPS formatted file is missing NAME line. The NAME must start in column 1 and be the first line in the file.';
   case 402, RCText='Incomplete MPS file. File probably truncated.';
   case 404, RCText='Row relationship error. Rerun with LISTINPUT YES.';
   case 406, RCText='Expecting a number. Rerun with LISTINPUT YES.';
   case 408, RCText='Bound relationship error. Rerun with LISTINPUT YES.';
   case 410, RCText='Split Column declaration. All rows a column intersects must be grouped together. Rerun with LISTINPUT YES. Or if your column names are case sensitive issue an XACLP with "Set CaseSensitive Yes".';
   case 412, RCText='In Column section, a column has duplicate row entries. Rerun with LISTINPUT YES, or if your row names are case sensitive issue an XACLP with "Set CaseSensitive Yes".';
   case 420, RCText='Column has a Power sequence with a negative lower bound, either correct lower bound or add XA_LBBASE status setting.';
   case 600, RCText='Misspelled or missing DBF table filename. Check for missing or incorrect path.';
   case 602, RCText='Error reading DBF table.';
   case 604, RCText='Special row and column name appear together.';
   case 606, RCText='Row increment setting; increment settings only apply to integer columns.';
   case 608, RCText='Row priority setting, priority settings only apply to integer columns.';
   case 610, RCText='Row field name not found in DBF table.';
   case 612, RCText='Column field name not found in DBF table.';
   case 614, RCText='Coefficient field name not found in DBF table.';
   case 616, RCText='No data available in DBF table.';
   case 700, RCText='No data previously loaded for solving.';
   case 702, RCText='Row name size not equal to previously loaded name size.';
   case 704, RCText='Column name size not equal to previously loaded size.';
   case 706, RCText='Specified name size does not match row or column name size.';
   case 900, RCText='Unknown error.';
   case 901, RCText='The first argument in the call to this function must be a pointer to the XA Environment Structure XAOSL. If the pointer address is correct then the XAOSL Structure data is corupted/destroyed. You must preserse the XAOSL Structure data from the time you call XAINIT until you call XADONE.';
   case 930, RCText='XA is solving an integer programming problem. A node is generated splitting the feasible region of each integer variables. The default maximum number of nodes is sqrt( number of Integers variables ) + number of 0/1 and semi-continuous columns + 4000 This default settings is too small. Use the XACLP command line parameter SET MAXNODES # to increase.';
   case 996, RCText='Missing or invalid c:\xa.lic file.';
   case 997, RCText='Email c:\xa.lic for activation.';
   case 998, RCText='Internal system administration error.';
   case 999, RCText='Your code has clobbered XA''s internal memory. Check and make sure you do not have any uninitialized pointers, arrays are mallocÂ’åä with the proper lengths, ...; the problem occurs in your code that was executed between the last two XA.... function calls.';
      
      %  --------------------------------------------
      %  Warning Message all begin with 2xxxx.
   case 20001, RCText='Ranging a free/null row with MPS file.';
   case 20030, RCText='Row name referenced in COLUMNS section is undefined.';
   case 20040, RCText='Column has no technological coefficients, XA is fixing with zero primal activity.';
   case 20045, RCText='Column name referenced in BOUNDS section is undefined.';
   case 20050, RCText='Unable to save "advance basis" solution. Disk is probably full or you do not have write access to it''s disk. Relocate this file on another disk with more space ore remove unnecessary files from disk.';
   case 20060, RCText='A free column only appears in a FREE row. Use the XACLP command line parameter SET MPSXCOMPATIBLE YES to FIX this column to zero. If you don''t XA may find your problem is unbounded.';
   case 20080, RCText='Duplicate row names in ROWS section.';
   case 20090, RCText='Unrecognized line in MPS formatted file.';
      
   otherwise
      if(rc>=902 & rc<= 919)
         RCText='919 (xx=02 thru 19) Argument xx has a unreasonable value. This can also occur when you leave off an argument. Check calling sequence or rerun with SET DEBUG YES to help identify arguments.';
      end
      
end

Result.XA.RCText = RCText;


% Result printing
if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / XA solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   
   fprintf('%s',ExitText);
   fprintf('\n');
   
   fprintf('\nObjective function at x (obj) %25.16f\n\n',f_k);
   %if mip
   %   fprintf('LP iterations%7d. ',lpiter);
   %else
   %   fprintf('Nodes visited%7d. ',glnodes);
   %end
   %fprintf('\n');
   %fprintf('Number of infeasibilities%7d. ',ninf);
   %fprintf('Sum of infeasibilities %e. ',sinf);
   %fprintf('\n');
end


if PriLev > 1
   if isempty(MAX_x)
      MAX_x=length(x_k);
   end
   fprintf('Optimal x = \n');
   xprinte(x_k(1:min(length(x_k),MAX_x)),'x:  ');
end


if PriLev > 3
   if isempty(MAX_c)
      MAX_c=20;
   end
   
   fprintf('Primal activities for all columns and rows\n');
   xprint(act,'act:  ',' %14.9f',5);

   fprintf('Dual activities for all columns and rows  = \n');
   xprint(dact,'dact: ',' %14.9f',5);
end

%if ~isempty(c)
%   if ~Fzero
%      Result.g_k=Prob.QP.F*x+c;
%   else
%      Result.g_k=c;
%   end
%elseif Fzero
%   Result.g_k=[];
%else
%   Result.g_k=Prob.QP.F*x;
%end
if mA > 0
   Result.Ax = Prob.A*x_k;
   Result = StateDef(Result, x_k, Result.Ax, [], ...
                     Prob.optParam.xTol, Prob.optParam.bTol, [], bl, bu, 1);
else
   Result.Ax = [];
   Result = StateDef(Result, x_k, [], [], ...
                     Prob.optParam.xTol, Prob.optParam.bTol, [], bl, bu, 1);
end

Result.c_k=[];
Result.cJac=[];

Result.Iter     = 1;

Result.Inform   = Inform;

%Result.MIP.ninf=ninf;
%Result.MIP.sinf=sinf;
%Result.MIP.slack=slack;
%Result.MIP.lpiter=lpiter;
%Result.MIP.glnodes=glnodes;
%Result.MIP.basis=basis;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 031016 ango Wrote file, based on cplexTL.m
% 031203 ango Corrected comments about xaControl and xaCtrl.
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040118 hkh  Additional clean-up
% 040528 hkh  Check on field ITERATIONS, use MaxIter if not set
% 040923 ango Omit unnecessary code for ITERATIONS - Move Prob.XA.xaControl
%             to Prob.MIP.xaControl
% 041011 ango Changes to exit texts.
% 041012 frhe Added IIS.
% 041013 med  Help updated.
% 041102 frhe Default TIMELIMIT 2e9 seconds. MaxCPU is not ignored anymore.
% 041103 frhe Default ITERATION now works. This way of setting defaults
%             made me have to change to case sensitive parameters in
%             tom2xaclp.
% 041109 frhe BIG is 1e20 instead of 1e10
% 041119 frhe BIG is 2e23 instead of 1e20
% 041126 ango XA 14 released
% 041202 hkh  Added call to StateDef, revised call to defblbu
% 050117 med  mlint review