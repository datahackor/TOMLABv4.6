% TOMLAB /XA LP, MILP and QP Solver
%
% -----------------------------------------------------------------------
%
%   xa.m is an interface routine to the XA solver which solves
%   continuous or mixed-integer linear programming problems (LP,MILP)
%   and also continuous quadratic programming problems (QP).
%
% ---  The LP and MILP problem types are defined as
%
%   minimize         c'x    subject to:	
%      x     x_L <=    x   <= x_U
%            b_L <=   Ax   <= b_U
%   where 
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the XA sparse matrix format
%   c, x_L, x_U has dimension n x 1 
%   b_L, b_U    has dimension m x 1
%
%   For MILP problems, a subset of the decision variables, x_i for
%   i in I are restricted to integer values only. 
%
% --- The QP problem type is defined as
%
%   minimize   0.5 * x'*F*x + c'x     subject to:	
%      x              x_L <=    x   <= x_U
%                     b_L <=   Ax   <= b_U
%   where 
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the XA sparse matrix format.
%   c, x_L, x_U has dimension n x 1
%   b_L, b_U has dimension m x 1
%   F is a n x n symmetric matrix, sparse or dense. 
%
%
%
% ------------------------------------------------------------------------
% function [x_k,f_k,Inform,modsts,solsts,act,dact,status,iis] = ...
%    xa(c,A,x_L,x_U,b_L,b_U,xaControl,callback,...
%       PriLev,LogFile,SaveFile,Prob,...
%       IntVars,VarWeight,SC,SC2,F,iisRequest);
%
% INPUT:  
% c         Linear objective function cost coeffs, n x 1.
% A         Linear constraint matrix, dense or sparse m x n matrix.
%           xa.m converts the matrix to a sparse format.
% x_L       Lower bounds on x. (if [], then x_L=0 assumed)
% x_U       Upper bounds on x
% b_L       Lower bounds on linear constraints
%
%    The following parameters are optional: 
%
% b_U       Upper bounds on linear constraints (if [], then b_U=b_L assumed)
%
% xaControl Struct array with fields for XA parameters. For example: 
%
%             xaControl.Barrier = 'Yes';   (also 1/0 for Yes/No is accepted)
%             xaControl.ITERATIONS = 1000; 
% 
%           tells XA to use the Barrier algorithm and send the maximum number
%           of iterations to 1000. 
%
% callback  0/1 vector of length 10. Defines the active callbacks during
%           the solving process. 
%
%           The callback routines and their corresponding indices in the 
%           callback vector are: 
%
%            1 - xacb_begin.m   - before solving
%            2 - xacb_infeas.m  - infeasible iteration
%            3 - xacb_feas.m    - feasible iteration
%            4 - xacb_node.m    - integer node generated
%            5 - xacb_intsol.m  - integer solution found
%            6 - xacb_branch.m  - user selects b & b variable
%            7 - xacb_barrier.m - barrier iteration  
%            8 - xacb_resolve.m - problem solve -> fast modify resolve
%            9 -                - currently not used
%           10 - xacb_end.m     - after solving 
%
%           The user can either modify the existing xacb_*.m
%           functions (use i.e. "which xacb_feas") to find them,
%           OR copies can be made and placed before the original
%           files in the MATLAB path. 
%
%           The calling syntax for all XA callbacks is
%
%               function xacb_(*)( xacbInfo, Prob )
%
%           and is described in more detail in xacb.m (help xacb.m)
%
% PriLev    Printing level in the XA m-file and XA C-interface.
%           = 0    Silent
%           = 1    Warnings and Errors
%           = 2    Summary information
%           = 3    More detailed information
%
%           > 10   Pause statements, and maximal printing (debug mode)
%
% LogFile   Name of file to write XA log to. If empty, no log is
%           written. 
%
% SaveFile  Name of file to write MPS representation of the problem to. The
%           name should be given without extension - .mps is added
%           automatically. If empty, the problem is not saved. 
%
% Prob      A structure. If TOMLAB calls XA, then Prob is the standard
%           TOMLAB problem structure, otherwise the user optionally may set:
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%
% IntVars   Defines which variables are integers, of general type I or binary B
%           Variable indices should be in the range [1,...,n].
%           IntVars is a single integer ==> Variable 1:IntVars are integer 
%           IntVars is a logical vector ==> x(find(IntVars > 0)) are integers 
%           IntVars is a vector of indices ==> x(IntVars) are integers 
%           (if [], then no integers of type I or B are defined)
%           XA checks which variables has x_L=0 and x_U=1, i.e. binary.
%
% VarWeight Vector of branching priorites for MILP problems. Should be a
%           n-vector, ideally of integers >=1. A lower value means higher
%           priority in the variable selection phase. 
%
% SC        A vector with indices for the Type 1 Semi-Continous variables,
%           i.e. that takes either the value 0 or a value in the range 
%           [ x_L(i) , x_U(i) ].
%
% SC2       A vector with indices for the Type 2 Semi-Continous variables,
%           i.e. that takes either the value of the corresponding upper bound,
%           or a value in the range [ x_L(i) , 0 ].
%
% F         Square dense or sparse matrix. Empty if non-quadratic problem.
%
% iisRequest A flag telling whether to obtain an IIS when a
%            problem is determined infeasible. The flag can be set
%            to one of the following values:
%
%             0 - Do not search for an IIS (default).
%             1 - Implementation of irreducible inconsistent systems
%                 (IIS) of constraints. Algorithm: IIS.
%             2 - Method of locationing a minimal number of constraints
%                 such that if all are removed the model is
%                 feasible. Algorithm: Block.
%
%            The IIS is returned through the output parameter: iis
%            Information about the IIS is automatically written to
%            the LogFile if a LogFile is defined.
%
% ------------------------------------------------------------------------------
%
% OUTPUT: 
%
% x_k       Solution vector with decision variable values (n x 1 vector)
% f_k       Objective function value at optimum
% Inform    Result of XA run
% modsts    Model status
% solsts    Solver status
% act       Primal activities for all columns+rows
% dact      Dual activities for all columns+rows
% status    Status of all rows+cols at optimum
% iis       Structure containing IIS information. 
%           Fields:
% 
%               iisStatus   Status flag. Possible values:
% 
%                            1 - IIS was obtained.
%                            0 - IIS was not requested.
%                           -1 - Problem is infeasible but no IIS found.
%                           -2 - Problem is not infeasible.
% 
%               iisMessage  Status message.
%               rowind      The row indices of the IIS set.
%
% -----------------------------------------------------------------------
%
% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 14.0$
% Written Oct 13, 2003.   Last modified Jan 17, 2005.
%

function [x_k,f_k,Inform,modsts,solsts,act,dact,status,iis] = ...
   xa(c,A,x_L,x_U,b_L,b_U,xaControl,callback,PriLev,LogFile, ...
      SaveFile,Prob,IntVars,VarWeight,SC,SC2,F, ...
      iisRequest)

if nargin < 18
  iisRequest = 0;
  if nargin < 17
    F = [];
    if nargin < 16
      SC2 = [];
      if nargin < 15
        SC = [];
        if nargin < 14
          VarWeight = [];
          if nargin < 13
            IntVars = [];
            if nargin < 12
              Prob = [];
              if nargin < 11
                PriLev = [];
                if nargin < 10
                  SaveFile = '';
                  if nargin < 9
                    LogFile = '';
                    if nargin < 8
                      callback = [];
                      if nargin < 7
                        xaControl = {};
                        if nargin < 6
                          b_U = [];
                          if nargin < 5
                            error('xa needs at least 5 arguments');
                          end,end,end,end,end,end,end,end,end, ...
                              end,end,end,end,end




% Determine # of variables
n = max([length(c),size(F,1),size(A,2)]);

if n==0
   error('xa: failed to determine problem dimension - c, F and A are all empty')
end

% Size checking 
if ~isempty(F)
   if(size(F,1)~=n)
      error(sprintf('xa.m: F has %d rows, expected %d',size(F,1),n));
   end
   if(size(F,2)~=n)
      error(sprintf('xa.m: F has %d columns, expected %d',size(F,2),n));
   end
end

[m,nA]=size(A);
if ~isempty(A)
   if nA~=n
      error(sprintf('xa.m: A has %d columns, expected %d',nA,n));
   end
end

if isempty(x_L),    x_L = zeros(n,1); end
if isempty(b_U),    b_U = b_L; end

if length(x_L)~=n
   error(sprintf('xa: Length of b_L = %d, expected %d\n',length(b_L),n));
end
if length(x_L)~=n
   error(sprintf('xa: Length of b_U = %d, expected %d\n',length(b_U),n));
end
if length(x_L)~=n
   error(sprintf('xa: Length of x_L = %d, expected %d\n',length(x_L),n));
end
if length(x_U)~=n
   error(sprintf('xa: Length of x_U = %d, expected %d\n',length(x_U),n));
end

if m==0 | isempty(c) 
   % Add a dummy constraint
   A = spalloc(1,n,1);
   A(1,1) = 1;
   b_L    = x_L(1);
   b_U    = x_U(1);
end

if isempty(PriLev), PriLev = 0; end

if isempty(IntVars) & isempty(SC) & isempty(SC2) & isempty(VarWeight)
   MIP=0;
else
   MIP=1;
end


if MIP
   % Vector for integers indicating type
   iv   = zeros(n,1);
   semi = zeros(n,1);
   
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
         error('xa: Illegal IntVars vector');
      end
      iv(IntVars)=1;
   end
   
   if isempty(VarWeight)
      prio = [];
   else
      prio = VarWeight(:);
      
      if length(prio) < n
         prio(end+1:n) = max(prio)+1;
      end
   end
   
   % Semi-continuous variables, between lower and upper bounds, 
   % or fixed at 0. 
   if ~isempty(SC)
      if any(SC < 1 | SC > n)
         SC
         error('xa: Illegal index in SC vector');
      end
      semi(SC) = 1;
   end
   
   % Semi-integer variables, between 0 and lower bound, 
   % or fixed at upper bound
   if ~isempty(SC2)    
      if any(SC2 < 1 | SC2 > n)
         SC2
         error('xa: Illegal index in SC2 vector');
      end
      semi(SC2) = 2;
   end
else
   iv   = [];
   prio = [];
   semi = [];
end

if( any(iv) & ~isempty(F) )
   error('xa: XA handles only continuous QP problems.');
end

% Check that Prob.P is set.
if isempty(Prob)
   Prob = struct('P',1);
end
if ~isfield(Prob,'P')
   Prob.P=1;
end


if length(callback) < 10
   callback(end+1:10) = 0;
end

if any(callback)
   % Define fields in Prob for use in callback m-files
   
   %    disp('xa.m: setup prob');
   %    Prob.A    = sparse(A);
   %    Prob.QP.c = c;
   %    Prob.QP.F = sparse(F);
   %    Prob.x_L  = x_L;
   %    Prob.x_U  = x_U;
   %    Prob.b_L  = b_L;
   %    Prob.b_U  = b_U;
   % 
   %    %Prob.qgtype    = gltype;
   %    %Prob.mgcols    = glcolidx;
   %    %Prob.dlim      = gllim;
   %    Prob.iv        = iv;
   %    Prob.qstype    = settype;
   %    Prob.msstart   = setbeg;
   %    Prob.mscols    = setcolidx;
   %    Prob.dref      = setref;
end

if ~isempty(LogFile)
   if ~ischar(LogFile)
      error('xa.m: input argument LogFile is not a string');
   end
end

if ~isempty(SaveFile)
   if ~ischar(SaveFile)
      error('xa.m: input argument SaveFile is not a string');
   end
end

% Compensation for XA definition of quadratic programs
% 041119 frhe Compensation not needed in latest XA14 release
%if ~isempty(F)
%    F = F*2;
%    for i =1:size(F,1)
%        F(i,i) = F(i,i) * 0.5;
%    end 
%    c = 2*c; 
%end

% If xaControl is a struct, xa.m has probably been called directly
% and we need to provide a cell array to the MEX:
if isstruct(xaControl)
   xaControl = tom2xaclp(xaControl,{});
end

[Inform,modsts,solsts,act,dact,status,iis]=xamex(c,sparse(F),sparse(A),b_L,b_U,x_L,x_U,...
   iv,semi,prio,PriLev,LogFile,SaveFile,xaControl,callback,Prob,iisRequest);

x_k = act(1:n);
f_k = act(n+1);

% 041119 frhe Compensation not needed in latest XA14 release
%if ~isempty(F)
%   f_k = 0.5*f_k;
%end
   
% if isempty(F)
%    f_k = c'*x_k;
% else
%    f_k = x_k(:)'*F*x_k(:) + c'*x_k;
% end


% MODIFICATION LOG:
%
% 031013 ango Wrote file
% 031027 ango Added callback handling and comments
% 031106 ango Comments updated
% 031121 ango Scaling of c for QP problems
% 031202 ango Fixed an erroneous ||
% 041012 frhe Added IIS.
% 041013 med  Help updated.
% 041109 frhe Added compensation for quadratic programming.
% 041119 frhe Commented compensation code. Not needed in latest release.
% 041126 ango XA 14 released
% 050117 med  mlint review