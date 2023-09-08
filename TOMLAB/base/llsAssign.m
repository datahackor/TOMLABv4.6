% llsAssign implements the TOMLAB (TQ) format for the linear least
% squares problem.
%
% llsAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% -----------------------------------------------------
%
% LLS minimization problem:
%
%
%        min   f(x) =  0.5 * ||Cx - d||.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of llsAssign:
%
% function Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
%                           t, weightType, weightY, ...
%                           A, b_L, b_U, ...
%                           x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least five parameters)
%
% C           Matrix m x n in objective ||Cx -y(t)||
% y           Vector m x 1 with observations in objective ||Cx -y(t)||
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as  a nx1 Inf vector.
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
%
% t           m x 1 vector with time values (if empty assumed to be 1:m)
% weightType  Type of weighting
% weightY     Vector of weights
%
% L I N E A R   C O N S T R A I N T S
% A           Matrix A in linear constraints b_L<=A*x<=b_U. Dense or sparse.
% b_L         Lower bound vector in linear constraints b_L<=A*x<=b_U.
% b_U         Upper bound vector in linear constraints b_L<=A*x<=b_U.
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
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Nov 8, 2000.    Last modified Jan 17, 2005.

function Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
                          t, weightType, weightY, ...
                          A, b_L, b_U, x_min, x_max, f_opt, x_opt)

if nargin < 16
   x_opt=[];
   if nargin < 15
      f_opt=[];
      if nargin < 14
         x_max=[];
         if nargin < 13
            x_min=[];
            if nargin < 12
               b_U=[];
               if nargin < 11
                  b_L=[];
                  if nargin < 10
                     A=[];
                     if nargin < 9
                        weightY=[];
                        if nargin < 8
                           weightType=[];
                           if nargin < 7
                              t=[];
                              if nargin < 6
                                 x_0=[];
end, end, end, end, end, end, end, end, end, end, end

n = size(C,2);
m = size(C,1);

Prob       = ProbDef(1);

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
Prob.LineParam.fLowBnd=0;

Prob.LS.C  = C;
if length(y) ~= m
   fprintf('Length of y %d, Columns in C are %d\n',length(y),m);
   error('Illegal length of y');
end
Prob.LS.y  = y(:);
if isempty(t)
   Prob.LS.t  = [1:m]';
else
   Prob.LS.t  = t(:);
end

if isempty(weightType)
   Prob.LS.weightType=0;
else
   Prob.LS.weightType=weightType;
end

Prob.LS.weightY = weightY;
Prob.LS.SepAlg=0;

Prob.LS.damp=0;
Prob.LS.L=[];

if issparse(C)
   Prob.JacPattern = spones(C);
end

[mA,mN]  = size(A);
Prob.mLin    = mA;
Prob.mNonLin = 0;

if mA > 0
   Prob.A = A;
   if mN ~= n
      fprintf('Number of variables %d\n',n);
      fprintf('Number of columns in linear constraint matrix A %d\n',mN);
      fprintf('These lengths should be the same, check input!!!\n');
      error('Illegal number of columns in A');
   end
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
   if issparse(A)
      Prob.ConsPattern  = spones(A);
   end
end

Prob.probType = checkType('lls');

Prob.probFile=0;

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

global MAX_x MAX_r % Max number of variables/constraints/resids to print

if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

Prob = tomFiles(Prob, 'ls_f', 'ls_g', 'lls_H', [], [], [], 'lls_r', 'lls_J');

% MODIFICATION LOG
%
% 001012  hkh  Written, based on probAssign
% 011204  hkh  Set Prob.f_Low = 0, Prob.LineParam.fLowBnd=0, Prob.LS.SepAlg=0;
% 011204  hkh  Set JacPattern and ConsPattern if sparse matrices
% 021216  hkh  Add fields damp=0 and L empty in Prob.LS
% 030918  ango Default parameters counted wrong, fixed.
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040414  hkh  Use spones to set sparsity patterns in JacPattern, ConsPattern
% 040429  hkh  Safeguard y in LS.y to be a column vector
% 040607  med  Help updates, extra checks removed
% 040728  med  tomFiles used instead
% 041201  hkh  Added check on number of columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  MAX_c removed as global, mlint check done
