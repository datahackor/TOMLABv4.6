% TOMLAB BARQP Sparse QP Interface
%
% Use tomRun to call BARQP: 
%
%   Result = tomRun('barqp',Prob);
%
% -----------------------------------------------------------------------------
%
% BARQP implements a sparse barrier method to solve the linear or quadratic  
% minimization problem defined as
%
%       min  0.5 * x' * F * x + c' * x
%        x
%       s/t   x_L <=   x  <= x_U
%             b_L <= A x  <= b_U
%
% where x,x_L,x_U are dense vectors of length n;
%       A is a sparse or dense matrix with dimension m*n and
%       b_L,b_U are dense vectors of length m;
%
% -----------------------------------------------------------------------------
%
% INPUT:  
%
% Prob        Problem structure in TOMLAB format. 
%             Use lpAssign or qpAssign to create the Prob structure.
% 
%             Fields used:
%
% x_L, x_U    Bounds on variables. 
%
% A           Linear constraint matrix, sparse (recommended) or dense 
% b_L, b_U    Bounds on linear constraints. 
%
% PriLev      Print level in the MEX interface.
%
% Prob.BOS    Structure with solver specific information. Fields used: 
%
%  options    Structure with options for the BARQP solver. The user sets
%             fields with names corresponding to the options he/she wishes 
%             to change. For example:
%
%             Prob.BOS.options.CONTOL = 1E-7
%             Prob.BOS.options.ALGOPT = 'FM'
%
%             The following keywords are recognized:
%
%             CONTOL, OBJTOL, PGDTOL, MAXNFE,
%             NITMIN, IT1MAX, ALFLWR, ALFUPR,
%             LYNPLT, LYNPNT, LYNVAR, NITMAX
%             IOFLAG, IOFLIN, IOFMFR, IOFPAT,
%             MAXLYN, TOLFIL, TOLKTC, TOLPVT,
%             ALGOPT, KTOPTN, IPOSTO, LYNFNC,
%             NEWTON, IRELAX, BIGCON, FEATOL,
%             PMULWR, PTHTOL, RHOLWR, IMAXMU,
%             MXQPIT, IOFSHR, MUCALC
%
%  PrintFile  Name of file to write general optimization output
%             to. The amount of information to print is controlled by the
%             following options: IOFLAG, IOFLIN, IOFSHR, IOFSRC, IOFMFR,
%             IOFPAT. If no PrintFile is given, but any of the
%             options mentioned above are nonzero, a default
%             PrintFile 'bosout.txt' will be created.
%
%  morereal   Number of elements to add to the double 'hold' vector. The
%             default size of 'hold' is 2,000,000 elements, which
%             sometimes is not enough for large problems.
%
%  moreint    Number of elements to add to the integer 'ihold' vector. The
%             minimum size of 'ihold' is 2,000,000 elements, which
%             sometimes is not enough for large problems. 
% 
%             The MEX interface can and does reallocate these vectors by
%             itself, but if reallocation occurs frequently, it is better to 
%             give a suitable value here. 
%
% -----------------------------------------------------------------------------
%
% OUTPUT:
%
% Result     Structure with optimization results. The following fields are
%            set:
%
%  x_k       Optimal point, if one has been found       
%  f_k       Objective function value at x_k
%
%  x_0       Starting point
%  f_0       Objective function value at x_0
%
%  g_k       Gradient of f(x) at final point x_k
%  H_k       Hessian of f(x) at final point x_k
%
%  Ax        Value of linear constraints A*x at x_k
%  c_k       Value of nonlinear constraints c(x) at x_k
%  cJac      Gradient of nonlinear 
%
%  ExitFlag  Flag telling if convergence or failure 
%  Inform    BARQP exit status
%
%  BOS       Structure with solver specific result information.
%
% -----------------------------------------------------------------------------
%
% Anders Goran, Tomlab Optimization Inc, E-mail: anders@tomlab.biz
% Copyright (c) 2004-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Nov 23, 2004. Last modified Jan 25, 2005.
%

function Result = barqpTL(Prob);

if nargin < 1
   error(['barqpTL needs two inputs, the Prob structure']);
end

Result = bosqpTL(Prob, 3);

% MODIFICATION LOG
% 
% 050125 frhe Created modification log to this new file.