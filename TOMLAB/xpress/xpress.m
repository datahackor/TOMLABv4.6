% TOMLAB /Xpress LP, QP, MILP and MIQP Solver
%
% -----------------------------------------------------------------------
%
%   xpress.m solves the following 
%   mixed integer (linear or quadratic) programming problem (MILP, MIQP):
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
%
% ------------------------------------------------------------------------
% function [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
%          glnodes, iis, sa] = xpress(c, A, x_L, x_U, b_L, b_U, ...
%          xpcontrol, callback, PriLev, Prob, ...
%          IntVars, PI, SC, SI, sos1, sos2, F, ...
%          LogFile, SaveFile, SaveMode, ...
%          iisRequest, iisFile, saRequest);
%
% INPUT:  
% c         Linear objective function cost coeffs, n x 1.
% A         Linear constraint matrix, dense or sparse m x n matrix.
%           xpress.m converts the matrix to a sparse format.
% x_L       Lower bounds on x. (if [], then x_L=0 assumed)
% x_U       Upper bounds on x
% b_L       Lower bounds on linear constraints
%
%    The following parameters are optional: 
%
% b_U       Upper bounds on linear constraints (if [], then b_U=b_L assumed)
%
% xpcontrol Structure, where the fields are set to the Xpress-MP control
%           parameters that the user wants to specify values for.
%           See the TOMLAB /Xpress user's guide for more information.
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
%
% callback(1)  m-file: xpcb_usn   User Select Node Callback
%         (2)  m-file: xpcb_upn   User Preprocess Node Callback
%         (3)  m-file: xpcb_uon   User Optimal Node Callback
%         (4)  m-file: xpcb_uin   User Infeasible Node Callback
%         (5)  m-file: xpcb_uis   User Integer Solution Callback
%         (6)  m-file: xpcb_ucn   User Node Cut-off Callback
%         (7)  m-file: xpcb_ucb   User Choose Branching Variable Callback
%         (8)  m-file: xpcb_il    Simplex Log Callback
%         (9)  m-file: xpcb_gl    Global Log Callback
%         (10) m-file: xpcb_bl    Barrier Log Callback
%         (11) m-file: xpcb_uop   User Output Callback
%         (12) m-file: xpcb_cmi   User Defined Cut Manager Init Routine
%         (13) m-file: xpcb_cms   User Defined Cut Manager Termination Routine
%         (14) m-file: xpcb_cm    User Defined Cut Manager Routine
%         (15) m-file: xpcb_tcm   User Defined Top Cut Manager Routine
%
% if callback is empty, by default callback(11) == 1.
% Then Xpress-MP error and warnings messages are printed in the Matlab
% command window. If xpcontrol.OUTPUTLOG == 1, then all Xpress-MP solver 
% information is printed in the Matlab command window. 
%
% NOTE! Currently only the first 11 callbacks are defined
%
% PriLev    Printing level in the xpress m-file and xpressmp C-interface.
%           = 0  Silent
%           = 1  Warnings and Errors
%           = 2  Summary information
%           = 3  More detailed information
%           > 10 Pause statements, and maximal printing (debug mode)
%
% Prob      Problem structure. If TOMLAB calls xpress through tomRun/tomSolve/xpressTL,
%           then Prob is the standard Tomlab problem structure, otherwise the 
%           user optionally may set: 
%           Prob.P = ProblemNumber, where ProblemNumber
%           is some integer (if input is [], then Prob.P=1 is set).
%           Prob.BIG = Big number representing infinity in Xpress/MP 
%           No numerical element should have abs value > BIG, see NOTE below
%
%           Any other desired information can also be passed to the
%           callback routines.
%
%           If any callback is enabled (see description of callback) then
%           problem arrays are set as fields in Prob, and the Prob structure
%           is always passed to the callback routines.
%
%           The defined fields are Prob.QP.c, Prob.QP.F, Prob.x_L,
%           Prob.x_U, Prob.A, Prob.b_L, Prob.b_U.  (if input is [],
%           then Prob.P=1 is set)
%
%           If Prob.MIP.KNAPSACK = 1 and callback(9) == 1, then the simple 
%           heuristic in xpcb_gl is used. If callback(9) is set, and 
%           Prob.MIP.KNAPSACK, or Prob.MIP is undefined, xpress is setting 
%           Prob.MIP.KNAPSACK = 0, to avoid the call to the heuristic. 
%          
%           NOTE: The xpress-mp MEX is using values >= Prob.BIG to define 
%           infinite values in lower/upper bounds of variables and constraints
%           By setting Prob.BIG to a nonempty positive big value, this value
%           will be used by the MEX.  By default, Prob.BIG is 1E20
%
%           DO NOT USE the MATLAB inf or NaN value in any arrays.
%           Convert any inf value to Prob.BIG (1E20), -inf to -Prob.BIG (-1E20) 
%
% IntVars   Defines which variables are integers, of general type I or binary B
%           Variable indices should be in the range [1,...,n].
%           IntVars is a single integer ==> Variable 1:IntVars are integer 
%           IntVars is a logical vector ==> x(find(IntVars > 0)) are integers 
%           IntVars is a vector of indices ==> x(IntVars) are integers 
%           (if [], then no integers of type I or B are defined)
%           xpress checks which variables has x_L=0 and x_U=1, i.e. binary.
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
% iisRequest Flag indicating whether to compute an IIS and return it to
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
% saRequest  Structure telling whether and how you want XPRESS to perform a
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
% -------------------------------------------------------------------------
% 
% OUTPUT: 
%
% x         Solution vector with decision variable values (n x 1 vector)
% slack     Slack variables (m x 1 vector)
% v         Lagrangian multipliers (dual solution vector) (m x 1 vector)
% rc        Reduced costs. Lagrangian multipliers for simple bounds on x.
% f_k       Objective function value at optimum
% ninf      Number of infeasibilities
% sinf      Sum of infeasibilities
% Inform    Result of Xpress-MP run
%           0 = Optimal solution found
%           2 = Unbounded solution
%           4 = Infeasible problem
%           5 = Some error occured
%           See the Xpress-MP problem attributes LPSTATUS (for LP) and
%           MIPSTATUS (for MIP) for more exact information. They are available 
%           in the global variable xpProblemAttrib, see below.
% basis     basis status of constraints + variables, (m + n x 1 vector)
% lpiter    Number of LP iterations
% glnodes   Number of nodes visited
%
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
% NOTE! After each call of xpress, defined are also the two global structures
%  global xpProblemAttrib
%  global xpControlVariables
%
%  The fields in these structures have all the values of the Xpress-MP
%  control variables and problem attributes, e.g.
%     xpControlVariables.DEFAULTALG    is the algorithm used
%     xpProblemAttrib.LPSTATUS         for the LP  status return
%     xpProblemAttrib.MIPSTATUS        for the MIP status return
%
% -----------------------------------------------------------------------
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release 2004.0.0 $
% Written July 8, 1999.   Last modified Jan 17, 2005.
%

function [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
         glnodes, iis, sa] = xpress(c, A, x_L, x_U, b_L, b_U, ...
         xpcontrol, callback, PriLev, Prob, IntVars, PI, SC, SI, ...
         sos1, sos2, F, LogFile, SaveFile, SaveMode, iisRequest,...
         iisFile, saRequest)
      
%#function xp2control xp2problem      

if nargin < 23
    saRequest = [];
if nargin < 22
   iisFile = '';
if nargin < 21
   iisRequest = [];
if nargin < 20
   SaveMode = '';
if nargin < 19
   SaveFile = '';
   if nargin < 18
      LogFile = '';
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
                                       xpcontrol = [];
                                       if nargin < 6
                                          b_U = [];
                                          if nargin < 5
                                             error('xpressmp needs at least 5 arguments');
end, end, end,
end, end, end, 
end, end, end, 
end, end, end,
end, end, end,
end, end, end,
end


global xpProblemAttrib xpControlVariables

% Set empty to avoid old values from last run to still be used
xpProblemAttrib=[];
xpControlVariables=[];

DEBUG=0;

% Safeguard for pure quadratic problem
if isempty(c) & ~isempty(F)
   n = size(F,1);
   c = zeros(n,1);
else
   c = c(:);
end

n = length(c);

if isempty(x_L),    x_L = zeros(n,1); end
if isempty(b_U),    b_U = b_L; end
if isempty(PriLev), PriLev = 0; end
if isempty(callback), callback=zeros(15,1); end

if PriLev>0
   if(callback(11)==0), callback(11)=1; end
   xpcontrol.OUTPUTLOG=1;
end

x_L = x_L(:);
x_U = x_U(:);
b_L = b_L(:);
b_U = b_U(:);

[m,nA] = size(A);

% Error checking 

if nA ~= n & m > 0
   fprintf('n = %d. Number of columns in A = %d\n',n,nA);
   error('xpress: Illegal length of A');
end

if length(b_L)~=m
   fprintf('m = %d. Length of b_L = %d\n',m,length(b_L));
   error('xpress: Illegal length of b_L');
end
if length(b_U)~=m
   fprintf('m = %d. Length of b_U = %d\n',m,length(b_U));
   error('xpress: Illegal length of b_U');
end
if length(x_L)~=n
   fprintf('n = %d. Length of x_L = %d\n',n,length(x_L));
   error('xpress: Illegal length of x_L');
end
if length(x_U)~=n
   fprintf('n = %d. Length of x_U = %d\n',n,length(x_U));
   error('xpress: Illegal length of x_U');
end

if ~isempty(IntVars) | ~isempty(sos1) | ~isempty(sos2) | ~isempty(PI) | ...
   ~isempty(SC) | ~isempty(SI)
   MIP=1;
else
   MIP=0;
end

if MIP
   callback(8)=0;  % Avoid simplex log callback for MIP

   % Logical vector for integers
   iv = zeros(n,1);

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
         error('xpress: Illegal IntVars vector');
      end
      iv(IntVars)=1;
   end

   % Semi-continuous variables, o and real interval
   if ~isempty(SC)
      if any(SC < 1 | SC > n)
         SC
         error('xpress: Illegal index in SC vector');
      end
      iv(SC) = 3;
   end
   % Semi-integer variables, 0 and integer interval
   if ~isempty(SI)
      if any(SI < 1 | SI > n)
         SI
         error('xpress: Illegal index in SI vector');
      end
      iv(SI) = 4;
   end
   % Partially integer variables
   if ~isempty(PI)
      if isfield(PI,'var')
         iv(PI.var) = 2;
      else
         PI.var = [];
      end
   else
      PI.var = [];
   end

   glcolidx = find(iv);

   ngl=length(glcolidx);
   if ngl > 0

      gltype = iv(glcolidx);

      gltype(gltype==1)=0; % Set as general Integer variables

      % Identify binary variables among the general Integer variables
      ix = find(x_L(glcolidx)==0 & (x_U(glcolidx)==1) & (gltype==0));

      gltype(ix)     = 1;

      gllim          = x_U(glcolidx);

      % Partially integer variables
      ix = find(gltype==2);
      if ~isempty(ix)
         gllim(ix)  = PI.lim;
      end

      if DEBUG
         xprinti(gltype,'gltype:');
         xprinte(gllim,'gllim:');
      end
   else
      gltype=[];
      gllim=[];
   end

   if isempty(sos1) & isempty(sos2)
      settype   = [];
      setbeg    = [];
      setcolidx = [];
      setref    = [];
   else
      ns1       = length(sos1);
      ns2       = length(sos2);
      nset      = ns1+ns2;
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
             error('xpress: Illegal sos1 input variable vector');
          end
          row       = sos1(i).row;
          if ~(row >= 0 & row <= m)
             fprintf('sos1 set %d. ',i);
             fprintf('Illegal row number  %d.',row);
             fprintf('\n');
             error('xpress: Illegal sos1 row data');
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
             error('xpress: Illegal sos2 input variable vector');
          end
          row       = sos2(i).row;
          if ~(row >= 0 & row <= m)
             fprintf('sos2 set %d. ',i);
             fprintf('Illegal row number  %d.',row);
             fprintf('\n');
             error('xpress: Illegal sos2 row data');
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
   callback([1:7,9])=0;
   gltype    = []; glcolidx  = []; gllim     = []; settype   = [];
   setbeg    = []; setcolidx = []; setref    = [];
end

% Check that Prob.P is set.
if isempty(Prob)
   Prob.P=1;
end
if ~isfield(Prob,'P')
   Prob.P=1;
end
if ~isfield(Prob,'MIP')
   Prob.MIP.KNAPSACK=0;
elseif ~isfield(Prob.MIP,'KNAPSACK')
   Prob.MIP.KNAPSACK=0;
end

%PriLev=2
%callback(11)=0
%callback(1:11)=1
%callback(9)=1

if any(callback)
   % Define fields in Prob for use in callback m-files
   if ~isfield(Prob,'QP')
      Prob.QP.c = c;
      Prob.x_L  = x_L;
      Prob.x_U  = x_U;
      Prob.b_L  = b_L;
      Prob.b_U  = b_U;
      if any(callback(1:7))
         % Only define matrices in struct if active callbacks are used
         Prob.QP.F = sparse(F);
         Prob.A    = sparse(A);
      end
   end

   Prob.qgtype    = gltype;
   Prob.mgcols    = glcolidx;
   Prob.dlim      = gllim;
   Prob.qstype    = settype;
   Prob.msstart   = setbeg;
   Prob.mscols    = setcolidx;
   Prob.dref      = setref;
end


% The xpress-mp MEX is using BIG as the value for testing the input for inf
% Note that xpress-mp is changing all values >= BIG to XPRS_PLUSINFINITY
% which is hard coded to 1E20
% Similar all values <= -BIG are set to XPRS_MINUSINFINITY (hard coded 1E20)

BIG=DefPar(Prob,'BIG',1E20);

% Call Xpress-MP using the xpressmp dll MEX interface!
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, iis, sa] ...
    = xpressmp(c, sparse(F), sparse(A), callback , x_L, x_U, b_L, b_U,...
      BIG, Prob.P, PriLev, Prob, ...
      gltype, glcolidx, gllim, settype, setbeg, setcolidx, setref, xpcontrol, ...
      SaveFile,SaveMode,LogFile, iisRequest, iisFile, saRequest);

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
% 980708 hkh  Written
% 981026 hkh  Clean up
% 010709 fre  Changed some things to get it work for r12 of xpress
% 010710 hkh  Remove Xpress path
% 010711 fre  Modified the call to xpressmp.
% 010713 hkh  Accept empty A
% 011211 hkh  Revision for R13
% 011213 hkh  Changed SC into SC and SI, semi-continuous and semi-integer
% 011227 hkh  Change SOS to sos
% 011229 hkh  Bug in Prob fields. Change field names to Release 13 names
% 020108 hkh  Add comments on callback(11) and xpcontrol.OUTPUTLOG
% 020628 hkh  Avoid setting the problem into Prob, if already defined. 
% 020628 hkh  Also avoid setting the problem, if no active callbacks used.
% 021202 ango Changed names of callbacks to lower case, for mcc.
% 040525 ango Add SaveFile and SaveMode parameters, change BIG to 1E10
% 040601 ango Add Logfile parameter
% 040604 ango Change order of last three parameters, to resemble cplex.m.
% 040608 med  Help fixes
% 040803 med  Added pragmas for MATLAB Compiler
% 040829 fhe  Added support for iis and sa.
% 040909 ango Set callback(11) if PriLev>0. 
% 041213 hkh  Set BIG hard coded as 1E20, if not Prob.BIG is set
% 041213 hkh  Add comments about BIG and inf
% 050117 med  mlint review
% 050128 med  DefPar added