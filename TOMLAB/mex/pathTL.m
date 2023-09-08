% TOMLAB PATH Solver
%
% function Result = pathTL(Prob) or tomRun('path', Prob)
%
% Solves LP, QP, LCP and MCP problems.
% See lcpAssign, mcpAssign for a problem description.
%
% INPUT:  
%
% Prob   Problem structure in TOMLAB format
%
% -----------------------------------------------
% Fields used in input structure Prob (call Prob=*Assign; to define Prob)
% -----------------------------------------------
%
% x_0       Starting values.
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% PriLevOpt Print level in solver.
%
% -----------------------------------------------
% Fields used in Prob.PATH:
% -----------------------------------------------
% M         For LCP problems. n x n matrix.
%
% q         Vector of linear functions. n x 1.
%
% mu        Starting value for multipliers.
%
% PrintFile Name of PATH Print file. Default: '' (silent)
%
% OptFile   Name of options file. Default: '' (no options file is read).
%           See the user's guide for available settings and syntax.
%
% -----------------------------------------------------------------------
%
% OUTPUT: 
%
% Result   Structure with results (see ResultDef.m):
%
% x_k      Solution vector.
% x_0      Initial solution vector.
% v_k      Lagrangian multipliers.
%
% f_k      Function value at optimum. (for LP and QP).
% r_k      Function values at optimum. (for LCP and MCP).
% J_k      Jacobian at optimum. (for LCP and MCP).
%
% xState   State of variables. Free == 0; On lower == 1; On upper == 2; 
%          Fixed == 3;
%
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2; 
%          Equality == 3;
%
% ExitFlag Exit status.
% ExitText Exit text from solver.
%
% Inform   PATH information parameter.
%
% Iter      Number of iterations
% MinorIter Number of minor iterations
% 
% FuncEv   Number of function evaluations. (MCP only)
% GradEv   Number of gradient evaluations. (MCP only)
% Solver   Name of the solver (PATH LCP or PATH MCP).
% SolverAlgorithm  Description of the solver.
%
%
% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2004 by Tomlab Optimization Inc., $Release 4.6.0 $
% Written Apr 20, 2004.  Last modified Dec 22, 2004.
%

function Result = pathTL(Prob)

if nargin < 1, error('pathTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

BIG=1E20;
probType = Prob.probType;

if probType == 18
    Prob.solvType = 18; % MCP solver
end

if any(probType == [17 8 2])
    Prob.solvType = 17; % LCP solver
end

Result=ResultDef(Prob);     % Define the Result struct

% Conversion of constraints.
% A vector t identifies if Ax<=b, Ax>=b, Ax=b

[bl, bu, n, nnLin, nnCon] = defblbu(Prob, BIG, 1);

if ~isempty(Prob.A)
   if isempty(Prob.b_L)
      Prob.AixUpp = find(~isinf(Prob.b_U));
      Prob.AixLow = [];
      Prob.AixEQ  = [];
   elseif isempty(Prob.b_U)
      Prob.AixLow = find(~isinf(Prob.b_L));
      Prob.AixUpp = [];
      Prob.AixEQ  = [];
   else
      ixEQ        = Prob.b_L==Prob.b_U & ~isinf(Prob.b_L);
      Prob.AixLow = find(~isinf(Prob.b_L) & ~ixEQ); % Ax>=b
      Prob.AixUpp = find(~isinf(Prob.b_U) & ~ixEQ); % Ax<=b
      Prob.AixEQ  = find(ixEQ);                     % Ax==b
   end
else
   Prob.AixEQ  = [];
   Prob.AixLow = [];
   Prob.AixUpp = [];
end

if ~isempty(Prob.A)
    % Redefine A, b and define the t vector
    %  -1: less than or equal
    %   0: equation
    %   1: greater than or equal
    
    A = Prob.A(Prob.AixLow,:);
    A = [A; Prob.A(Prob.AixUpp,:)];
    A = [A; Prob.A(Prob.AixEQ,:)];
    if ~issparse(A)
        A = sparse(A);
    end
    
    t(1:length(Prob.AixLow)) = 1;
    t(length(Prob.AixLow)+1:length(Prob.AixLow)+length(Prob.AixUpp)) = -1;
    t(length(Prob.AixLow)+length(Prob.AixUpp)+1:length(Prob.AixLow)+length(Prob.AixUpp)+length(Prob.AixEQ)) = 0;
    
    b = Prob.b_L(Prob.AixLow);
    b = [b; Prob.b_U(Prob.AixUpp)];
    b = [b; Prob.b_U(Prob.AixEQ)];
    
    m = size(A,1);
else
    m = 0;
    A = [];
end

if any(probType == [17 8 2]) % LCP, LP, QP
    Prob=iniSolve(Prob,17,0,0);    % LCP, LP or QP problem.
    Result.Solver='PATH LCP'; % Solver name
else
    Prob=iniSolve(Prob,18,1,0);    % MCP problem.
    Result.Solver='PATH MCP'; % Solver name
end

Result.SolverAlgorithm='Mixed Complimentary Solver PATH 4.6';
PriLev=Prob.PriLevOpt;      % Printing level in solver

% Initial checks on the inputs

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end


x_L = bl(1:n);
x_U = bu(1:n);

% Extract x_0 from Prob structure
x_0 = Prob.x_0;
if isempty(x_0)
    x_0 = zeros(n,1);
end
x_0 = min(max(x_0,x_L), x_U);

if probType == 8     % If LP problem
   q = Prob.QP.c;
   M = sparse(length(q),length(q));
   mu = zeros(m,1);
end
if probType == 2     % If QP problem
   q = Prob.QP.c;
   M = sparse(Prob.QP.F);
   mu = zeros(m,1);
end
if probType == 17    % If LCP problem
   M = sparse(Prob.PATH.M);
   q = Prob.PATH.q;
   mu = full(Prob.PATH.mu(:));  % extract and initialize mu
end
if probType == 18    % If MCP problem
   mu = full(Prob.PATH.mu(:));  % extract and initialize mu
end

if ~isempty(A)
   x_L_p = -BIG*ones(m,1);
   x_U_p =  BIG*ones(m,1);
   
   x_L_p(find(t > 0)) = 0;
   x_U_p(find(t < 0)) = 0;
   mu = min(max(mu,x_L_p),x_U_p);
   
   if probType ~= 18
       M = [M -A'; A sparse(m,m)];   
       q = [q; -b];
   end
   
   x_0 = [x_0; mu];
   x_L = [x_L; x_L_p];
   x_U = [x_U; x_U_p];
end

idx = find(x_L > x_U);
if length(idx) > 0
   error('Bounds infeasible.');
end

row = n + m;

PATH=DefPar(Prob,'PATH');
% Printfile and Optfile
PrintFile = DefPar(PATH,'PrintFile','');
OptFile   = DefPar(PATH,'OptFile','');

switch probType
   case {17, 2, 8} % LCP, QP, LP problem.
      nnzJ = nnz(M);
      [Inform, ttime, f_k, basis, x_k, iw] = lcppath(...
         row, nnzJ, x_0, x_L, x_U, M, q, PrintFile, PriLev, OptFile);
      
      v_k=[];
      if ( m > 0 )
         v_k = x_k(n+1:n+m);
         x_k = x_k(1:n);
      end
      
      % MCP only
      f = []; J = [];
      
      % iw(1) : residual             Value of residual at final point.
      % iw(2) : major_iterations     Major iterations taken.
      % iw(3) : minor_iterations     Minor iterations taken.
      % iw(4) : crash_iterations     Crash iterations taken.
      % iw(5) : function_evaluations Function evaluations performed.
      % iw(6) : jacobian_evaluations Jacobian evaluations performed.
      % iw(7) : gradient_steps       Gradient steps taken.
      % iw(8) : restarts             Restarts used.  

      Result.Iter      = iw(2);
      Result.MinorIter = iw(3);
      Result.PATH.iw = iw;
      
     case 18
        
        [f,J,domerr] = feval('mcp_funjac',x_0+1e-5*ones(size(x_0))+1e-5*abs(x_0),1,Prob);
        
        if (domerr > 0)
           [f,J,domerr] = feval('mcp_funjac',x_0,1,Prob);
        end
        
        if (domerr > 0)
           error(['mcp_funjac not defined at starting point']);
        end
        
        if ~issparse(J)
           error(['mcp_funjac must return a sparse Jacobian']);
        end
        
        nnzJ = nzmax(J);
        ele = nnzJ + 2*nzmax(A);
        
        if m > 0            
           Prob.PATH.fJ = 'mcp_funjac';
           Prob.PATH.A = A;
           Prob.PATH.b = b;
           Prob.PATH.m = m;
            
           [Inform, ttime, f, J, x_k] = mcppath(...
              row, ele, x_0, x_L, x_U, 'mcp_vifunjac', Prob, PrintFile, PriLev, OptFile);            
        else
           [Inform, ttime, f, J, x_k] = mcppath(...
              row, ele, x_0, x_L, x_U, 'mcp_funjac', Prob, PrintFile, PriLev, OptFile);
        end

        v_k = [];
        
        if m > 0
           v_k = x_k(n+1:n+m);
           x_k = x_k(1:n);
           
           J = J(1:n,1:n);
           f = f(1:n) + A'*v_k;
        end  
        
     otherwise
        error('Invalid problem type: use lp, qp, lcp or mcpAssign');
end

Result.Prob = Prob;
Result.Solver = 'PATH';
Result.x_k = x_k;
Result.x_0 = Prob.x_0;
Result.v_k = v_k;

optParam = Prob.optParam;

if m > 0
   Result = StateDef(Result, x_k(1:n), Prob.A*x_k(1:n), [], ...
                  optParam.xTol, optParam.bTol, [], bl, bu, 1);
else
   Result = StateDef(Result, x_k(1:n), [],[], optParam.xTol, [],[], bl, bu, 1);
end

if any(probType == [17 8 2])
   Result.SolverAlgorithm = 'PATH LCP Interface';
   if probType == 8
      Result.f_k = Prob.QP.c(:)'*x_k;
   end
   if probType == 2
      Result.f_k = x_k(:)'*Prob.QP.F*x_k(:)/2+Prob.QP.c(:)'*x_k;
   end
else
   Result.r_k = f;
   Result.J_k = J;
   
   global n_r n_J
   Result.FuncEv = n_r;
   Result.GradEv = n_J;
   Result.SolverAlgorithm = 'PATH MCP Interface';
end

ExitFlag = 99;

% Status code:
switch(Inform)
   case {1,2}
      ExitFlag=0;    % OK
   case {3,4,5}
      ExitFlag=1;    % Iter / time limit
   case {7,8,9}
      ExitFlag=4;    % Infeasibility / domain errors
   case 6
      ExitFlag=11;   % User request stop
   case 10
      ExitFlag=10;   % Internal error
   otherwise
      ExitFlag=-1;   % Should not happen.
end

% Once more, for the exit text
switch(Inform)
   case 1,
      ExitText='Optimal solution found';
   case 2,
      ExitText='No progress';
   case 3,
      ExitText='Major iteration limit reached';
   case 4,
      ExitText='Minor iteration limit reached';
   case 5,
      ExitText='Time limit reached';
   case 6,
      ExitText='User requested interrupt';
   case 7,
      ExitText='Bound error';
   case 8,
      ExitText='No starting point found';
   case 9,
      ExitText='Problem is infeasible';
   case 10,
      ExitText='Internal error';
   otherwise
      ExitText='Unknown status value';
end

Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 040418 med  Created PATH interface.
% 040419 med  Added LP and QP support.
% 040421 ango Minor fix for unconstrained problems.
% 040421 ango Iterations reported. Adapted to modified lcppath.
% 040503 ango Complete status flag handling added.
% 041202 hkh  Revise calls to defblbu and StateDef
