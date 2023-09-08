% clsAssign implements the TOMLAB (TQ) format for mixed-integer unconstrained
% and constrained nonlinear least squares problems.
% It is also suitable to define vector valued function problems with
% a corresponding Jacobian matrix as derivative.
% For example minimax problems are solved with infSolve after using clsAssign
% to define the problem.
% L1 fitting problems are solved with L1Solve after using clsAssign
% to define the problem.
%
% clsAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% -----------------------------------------------------
%
% Mixed-integer constrained nonlinear least squares problem:
% 
%
%        min    f(x) = 0.5 * (r' * r),  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% Linear    equality equations: Set b_L==b_U
% Nonlinear equality equations: Set c_L==c_U
% Fixed     variables:          Set x_L==x_U. Both x_L and x_U must be finite.
% x(IntVars) are integer values, IntVars is an index set, a subset of [1:n].
%
% -----------------------------------------------------
%
% Syntax of clsAssign:
%
% function Prob = clsAssign(r, J, JacPattern, x_L, x_U, Name, x_0, ...
%                           y, t, weightType, weightY, SepAlg, fLowBnd, ...
%                           A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
%                           x_min, x_max, f_opt, x_opt, ...
%                           IntVars, VarWeight, fIP, xIP);
%
% INPUT (Call with at least eight parameters)
%
% r           Name of function that computes the residual 
% J           Name of function that computes the Jacobian m x n - matrix 
% JacPattern  m x n zero-one sparse or dense matrix, where 0 values indicate 
%             zeros in the Jacobian and ones indicate values that might 
%             be non-zero. If empty indicates estimation of all elements 
%             JacPattern is used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, JacPattern==[]
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as  a nx1 Inf vector.
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector 
%
% Note:       The number n of the unknown variables x are taken as
%             max(length(x_L),length(x_U),length(x_0))
%             You must specifiy at least one of these with correct length,
%             then the others are given default values
%
% y           m x 1 vector with observations y(t) to be fitted
% Note 1:     In the nonlinear least squares computations, the length of y
%             is used, and if weightType == 1, the y values are used to
%             weight the residual. However, the user should still subtract
%             the observations y from the model in the residual computations
% Note 2:     If y(t) is not used, then provide y = zeros(m,1), because Tomlab
%             needs the length m to setup memory correct for some solvers
%
% -----------------------------------------------------------------------------
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
% -----------------------------------------------------------------------------
%
% t           m x 1 vector with time values
% weightType  Type of weighting
%             0 = No weights
%             1 = weight with observation vector y
%             2 = weight with user given weight vector/matrix in Prob.LS.weightY
%                 Must either have length m x 1 or m x m (could be sparse)
%             3 = A user given function defined in Prob.LS.weightY will be
%                 called to compute the weights, defined as for weightType=2
%                 The function, say "DefWht" should be defined as
%                 function wLS = DefWht(x,r,Prob)
%                 function wLS = DefWht(x,r) or function wLS = DefWht(x)
% weightY     Vector of weights (or function defining the weights)
% SepAlg      Flag if to use separable nonlinear least squares
%
% fLowBnd     A lower bound on the function value at optimum. Default 0
%             A good estimate is not critical. Use [] if not known at all.
%
% L I N E A R   C O N S T R A I N T S
% A           Matrix A in linear constraints b_L<=A*x<=b_U. Dense or sparse.
% b_L         Lower bound vector in linear constraints, b_L<=A*x<=b_U. 
% b_U         Upper bound vector in linear constraints, b_L<=A*x<=b_U. 
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints 
% dc          Name of function that computes the constraint Jacobian mN x n 
% c_L         Lower bound vector in nonlinear constraints, c_L<=c(x)<=c_U. 
% c_U         Upper bound vector in nonlinear constraints, c_L<=c(x)<=c_U. 
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate 
%             zeros in the constraint Jacobian and ones indicate values that 
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
%
% A D D I T I O N A L   P A R A M E T E R S
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.
%
% A D D I T I O N A L MIP P A R A M E T E R S
%
% If solving a mixed-integer nonlinear least squares problems
%
% IntVars     The set of integer variables. Can be given in one of three ways:
%
%             1) a scalar N<=n, in which case variables x(1)-x(N) are 
%                restricted to integer values.
%
%             2) a vector of indices, e.g. [1 2 5]
%
%             3) a 0-1 vector of length n=length(x) where nonzero elements    
%                indicate integer variables
%
% VarWeight   Priorities for each variable in the variable selection phase
%             A lower value gives higher priority. 
%
% fIP         An upper bound on the IP value wanted. Makes it possible to
%             cut branches and avoid node computations.
% xIP         The x-values giving the fIP value.
%
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Oct 12, 2000.    Last modified Jan 17, 2005.

function Prob = clsAssign(r, J, JacPattern, x_L, x_U, Name, x_0, ...
                          y, t, weightType, weightY, SepAlg, fLowBnd, ...
                          A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
                          x_min, x_max, f_opt, x_opt, ...
                          IntVars, VarWeight, fIP, xIP)

if nargin < 29
   xIP=[];
   if nargin < 28
      fIP=[];
      if nargin < 27 
         VarWeight=[];
         if nargin < 26
            IntVars=[];
if nargin < 25
   x_opt=[];
   if nargin < 24
      f_opt=[];
      if nargin < 23
         x_max=[];
         if nargin < 22
            x_min=[];
            if nargin < 21
               c_U=[];
               if nargin < 20
                  c_L=[];
                  if nargin < 19
                     ConsPattern=[];
                     if nargin < 18
                        dc=[];
                        if nargin < 17
                           c=[];
                           if nargin < 16
                              b_U=[];
                              if nargin < 15
                                 b_L=[];
                                 if nargin < 14
                                    A=[];
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

if nargin < 13
   fLowBnd=[];
   if nargin < 12
      SepAlg=[];
      if nargin < 11
         weightY=[];
         if nargin < 10
            weightType=[];
            if nargin < 9
               t=[];

end, end, end, end, end, 

n = max([length(x_L),length(x_U),length(x_0)]);

Prob       = ProbDef(1);
Prob.P     = 1;
Prob.N     = n;
Prob.Name  = deblank(Name);

if isempty(x_min)
   Prob.x_min = -1*ones(n,1);
else
   Prob.x_min = x_min;
end
if isempty(x_max)
   Prob.x_max = 1*ones(n,1);
else
   Prob.x_max = x_max;
end

Prob.f_opt = f_opt;
Prob.x_opt = x_opt;
Prob.f_Low = 0;

Prob.LS.y  = y(:);
if isempty(y)
   disp('clsAssign: WARNING - empty y, e.g. solver NLSSOL will not work');
end
Prob.LS.t  = t(:);

if isempty(weightType) 
   Prob.LS.weightType=0;
else
   Prob.LS.weightType=weightType;
end
if isempty(SepAlg) 
   Prob.LS.SepAlg=0;
else
   Prob.LS.SepAlg=SepAlg;
end

Prob.LS.weightY = weightY;
Prob.JacPattern   = JacPattern;
Prob.ConsPattern  = ConsPattern;

if isempty(fLowBnd) 
   Prob.LineParam.fLowBnd=0;
else
   Prob.LineParam.fLowBnd=max(0,fLowBnd); 
end

if isempty(x_0)
   Prob.x_0 = zeros(n,1);
elseif length(x_0) ~= n
   fprintf('Length of x_0 %d, should be %d\n',length(x_0),n);
   error('Illegal length of x_0');
else
   Prob.x_0 = x_0(:);
end

if isempty(x_L) 
   Prob.x_L = -Inf*ones(n,1);
elseif length(x_L) ~= n
   fprintf('Length of x_L %d, should be %d\n',length(x_L),n);
   error('Illegal length of x_L');
else
   Prob.x_L = x_L(:);   
end

if isempty(x_U) 
   Prob.x_U = Inf*ones(n,1);
elseif length(x_U) ~= n
   fprintf('Length of x_U %d, should be %d\n',length(x_U),n);
   error('Illegal length of x_U');    
else
   Prob.x_U = x_U(:);
end

[mA,mN] = size(A);
Prob.mLin = mA;

if mA > 0
   if mN ~= n
      fprintf('Number of variables %d\n',n);
      fprintf('Number of columns in linear constraint matrix A %d\n',mN);
      fprintf('These lengths should be the same, check input!!!\n');
      error('Illegal number of columns in A');
   end
   Prob.A = A;
   if isempty(b_L)
      Prob.b_L=-Inf*ones(mA,1);
   elseif mA == length(b_L)
      Prob.b_L=b_L(:);
   else
      fprintf('Length of b_L %d, Rows in A are %d\n',length(b_L),mA);
      error('Illegal length of b_L');
   end
   if isempty(b_U) 
      Prob.b_U=Inf*ones(mA,1);
   elseif mA == length(b_U)
      Prob.b_U=b_U(:);
   else
      fprintf('Length of b_U %d, Rows in A are %d\n',length(b_U),mA);
      error('Illegal length of b_U');
   end
end

mN = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if ~isempty(IntVars)
   Prob.probType = checkType('minlp');
elseif mA+mN > 0
   Prob.probType = checkType('cls');
else
   Prob.probType = checkType('ls');
end

Prob.probFile=0;

if mN > 0
   if isempty(c_L) 
      Prob.c_L=-Inf*ones(mN,1);
   else
      Prob.c_L=c_L(:);
   end
   if isempty(c_U) 
      Prob.c_U=Inf*ones(mN,1);
   else
      Prob.c_U=c_U(:);
   end
end
if mN == 0 & ~isempty(c)
   fprintf('WARNING in clsAssign!!! ')
   fprintf('Constraint c is given. But no lower or upper bounds.\n')
end

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

Prob.MIP = struct('IntVars',IntVars, 'VarWeight',VarWeight, ...
                  'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP, ...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

%Prob.PartSep.pSepFunc= 1;
%Prob.PartSep.index   = 0;

Prob = tomFiles(Prob, 'ls_f', 'ls_g', 'ls_H', c, dc, [], r, J);

% MODIFICATION LOG
%
% 001012  hkh  Written, based on probAssign
% 011204  hkh  Set Prob.f_Low = 0
% 020409  hkh  Improve comments for vector valued functions
% 020411  hkh  Add comments about L1
% 030210  hkh  Errors in comments about c_L and c_U
% 031201  hkh  Check if c defined and no bounds given, issue warning
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040330  hkh  Demand input y, 8 parameters. Display warning if y empty
% 040414  hkh  Expand to handle MIP parameters, mixed-integer NLLS
% 040607  med  Help updates
% 040728  med  tomFiles used instead
% 041201  hkh  Added check number of columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  mlint review