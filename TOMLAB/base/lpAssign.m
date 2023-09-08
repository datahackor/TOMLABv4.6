% lpAssign is a direct way of setting up a Linear Programming (LP) problem
% in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = lpAssign(...)
%
% It is then possible to solve the LP problem using the TOMLAB LP solver
% lpSolve with the call:  Result = lpSolve(Prob);
%
% or any quadratic or general constrained solver:
%
% Result = tomRun('qpSolve', Prob, 1);
% Result = tomRun('conSolve', Prob, 1);
% Result = tomRun('minos', Prob, 1);
%
% It is also possible to run the linprog.m interface, similar to 
% linprog in Optimization Toolbox.
%
% See the file tomlab\examples\testlinprog.m for an example
%
% lpAssign may also create an Init File in the TOMLAB Init File format,
% see the input argument setupFile.
% -----------------------------------------------------
%
% LP minimization problem:
% 
%
%        min    c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
%
% Syntax of lpAssign:
%
% function Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
%                 setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (One parameter c must always be given)
%
% c            The vector c in c'x in the objective function
% A:           The linear constraint matrix 
% b_L:         The lower bounds for the linear constraints
% b_U:         The upper bounds for the linear constraints
% x_L:         Lower bounds on x
% x_U:         Upper bounds on x
%
%              b_L, b_U, x_L, x_U must either be empty or of full length
%
% x_0:         Starting point x (may be empty)
% Name         The name of the problem (string)
% setupFile    The (unique) name as a TOMLAB Init file. If nonempty lpAssign
%              will create a executable m-file with this name and the given
%              problem defined as the first problem in this file.
%              See lp_prob.m, the TOMLAB predefined LP Init File.
%              If empty, no Init File is created.
% nProblem     Number of problems to predefine in the setupFile
%              Not used if setupFile is empty.
%              If empty assumed to be one. Then text are included in the
%              setupFile on how to create additional problems.
% fLowBnd      A lower bound on the function value at optimum. Only used if
%              running nonlinear TOMLAB solvers with line search.
% x_min        Lower bounds on each x-variable, used for plotting
% x_max        Upper bounds on each x-variable, used for plotting
% f_opt        Optimal function value(s), if known (Stationary points)
% x_opt        The x-values corresponding to the given f_opt, if known.
%              If only one f_opt, give x_opt as a 1 by n vector
%              If several f_opt values, give x_opt as a length(f_opt) x n matrix
%              If adding one extra column n+1 in x_opt, 
%              0 indicates min, 1 saddle (nonlinear problems), 2 indicates max.
%              x_opt and f_opt is used in printouts and plots.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.6.0$
% Written July 2, 1999. Last modified Jan 17, 2005.

function Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
                setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 15
   x_opt=[];
   if nargin < 14
      f_opt=[];
      if nargin < 13
         x_max=[];
         if nargin < 12
            x_min=[];
            if nargin < 11
               fLowBnd=[];
               if nargin < 10
                  nProblem=[];
                  if nargin < 9
                     setupFile=[];
                     if nargin < 8
                        Name=[];
                        if nargin < 7
                           x_0=[];
                           if nargin < 6
                              x_U=[];
                              if nargin < 5
                                 x_L=[];
                                 if nargin < 4
                                    b_U=[];
                                    if nargin < 3
                                       b_L=[];
                                       if nargin < 2
                                          A=[];
                                             if nargin < 1
                           error('lpAssign requires at least one parameter c'); 
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

if isempty(nProblem)
   nProblem=1;
elseif nProblem < 1
   nProblem=1;
end

probType = checkType('lp');

if nargout > 0

Prob=ProbDef(1);

Prob.QP.c = c(:);
Prob.P    = 1;
n         = max(size(A,2),length(c));
Prob.N = n;

Prob.probType = probType;
Prob.probFile = 0;

if ~isempty(fLowBnd), Prob.LineParam.fLowBnd=fLowBnd; end 

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

if ~isempty(A) 
   Prob.A=A;
   [m,mN]=size(A);
   if mN ~= n
      fprintf('Number of variables %d\n',n);
      fprintf('Number of columns in linear constraint matrix A %d\n',mN);
      fprintf('These lengths should be the same, check input!!!\n');
      error('Illegal number of columns in A')
   end
   Prob.mLin = m;
   if isempty(b_L) 
      Prob.b_L=-Inf*ones(m,1);
   elseif m==length(b_L)
      Prob.b_L=b_L(:);
   else
      fprintf('Length of b_L %d, Rows in A are %d\n',length(b_L),m);
      error('Illegal length of b_L'); 
   end
   if isempty(b_U) 
      Prob.b_U=Inf*ones(m,1);
   elseif m==length(b_U)
      Prob.b_U=b_U(:);
   else
      fprintf('Length of b_U %d, Rows in A are %d\n',length(b_U),m);
      error('Illegal length of b_U');
   end
else
   Prob.mLin = 0;
end

Prob.mNonLin = 0;

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');

end

[FName, Name] = PrintAssign(probType, Name, setupFile, nProblem);

F=[];

if ~isempty(FName)
   save(FName,'setupFile','Name','F','c','A','b_L','b_U','x_0','x_L','x_U', ...
        'x_min','x_max','x_opt','f_opt');
else
   Prob.Name = deblank(Name);
end

% MODIFICATION LOG
%
% 990914  hkh  Modify for new automatic setup file generation
% 010528  hkh  Add definition of problem name if no IF file created
% 030524  hkh  Use tomFiles to define the Prob.USER structure
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040607  med  problemName to Name, help fixed for Init File
% 040728  med  tomFiles used instead
% 041123  hkh  Change call to tomRun in help
% 041201  hkh  Add check on columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0