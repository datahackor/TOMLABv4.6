%		GetSolver.m
%
% function Solver = GetSolver(Type, LargeScale, Convex)
% 
% GetSolver returns the TOMLAB default solver for different problems
% All standard Types of optimization and also type: FP and DLP, see below
%
% The general minimization problem of QP type is:
% 
%
%        min   0.5 * x' * F * x + c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% INPUT PARAMETERS
% Type:      String with type of problem:
%  'qp'      Quadratic programming, F is nonempty (QP)
%  'lp'      Linear programming, F is empty, c is nonempty, and c~=0 (LP)
%  'fp'      Feasible point (phase 1) linear programming, F and c empty or 0
%  'dlp'     Dual linear programming. A standard LP problem with a 
%            dual feasible initial point available.
% Also Type can be any of the following standard types
%  'uc'      Unconstrained optimization
%  'con'     Nonlinear programming (constrained optimization) (NLP)
%  'ls'      Nonlinear least squares (NLLS)
%  'lls'     Linear least squares (LLS)
%  'cls'     Constrained nonlinear least squares
%  'mip'     Mixed-integer (linear) programming (MIP or MILP)
%  'glb'     Global optimization (box-bounded)
%  'glc'     Global optimization (box-bounded, integer and constrained)
%  'miqp'    Mixed-integer quadratic programming (MIQP)
%  'minlp'   Mixed-integer nonlinear programming (MINLP)
%  'sdp'     Semidefinite programming (SDP) - Linear SDP with LMI constraints
%  'bmi'     Linear SDP with BMI constraints (BMI)
%  'exp'     Exponential sum fitting
%            
% LargeScale If the flag Prob.LargeScale > 0 is set, LargeScale is set true
%
% Convex     If Convex > 0, problem is convex, i.e. F is positive semidefinite
%            for QP problems. Only used for type QP.
%
% OUTPUT PARAMETERS
% Solver     String with name of default solver
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Nov 6, 2000.    Last modified Dec 21, 2004.
%

function Solver = GetSolver(Type, LargeScale, Convex)

if nargin < 3
   Convex = 1;
   if nargin < 2
      LargeScale = 0;
   end
end

[TomV,os,TV] = tomlabVersion;

switch lower(Type)
  case 'qp'
     % Default TOMLAB QP solver
     % QPOPT
     if Convex & LargeScale
        if TV(9)
           Solver = 'CPLEX';
        elseif TV(4)
           Solver = 'sqopt';
        elseif TV(2)
           Solver = 'qp-minos';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(8)
           Solver = 'xpress-mp';
        elseif TV(3)
           Solver = 'qpopt';
        else
           Solver = 'qld';
        end
     elseif Convex  & ~LargeScale
        if TV(9)
           Solver = 'CPLEX';
        elseif TV(3)
           Solver = 'lssol';
        elseif TV(2)
           Solver = 'qpopt';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(4)
           Solver = 'sqopt';
        elseif TV(8)
           Solver = 'xpress-mp';
        else
           Solver = 'qld';
        end
     elseif ~Convex &  LargeScale
        if TV(4)
           Solver = 'snopt';
        elseif TV(2)
           Solver = 'qp-minos';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(3)
           Solver = 'lssol';
        else
           Solver = 'qpSolve';
        end
     else
        if TV(2)
           Solver = 'qpopt';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(7)
           Solver = 'bqpd';
        else
           Solver = 'qpSolve';
        end
     end
  case 'lp'
    if TV(2)
       if LargeScale
          %Solver = 'LP-MINOS';
          Solver = 'MINOS';
       else
          Solver = 'MINOS';
          % Solver = 'lpopt'; % Still MINOS is more reliable
       end
    elseif TV(9)
       Solver = 'CPLEX';
    elseif TV(8)
       Solver = 'xpress-mp';
    else
       Solver = 'qld';
    end
  case 'fp'
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(2)
       %Solver = 'LP-MINOS';
       Solver = 'MINOS';
    elseif TV(8)
       Solver = 'xpress-mp';
    else
       Solver = 'qld';
    end
  case 'dlp'
    if TV(2)
       if LargeScale
          %Solver = 'LP-MINOS';
          Solver = 'MINOS';
       else
          Solver = 'MINOS';
          % Solver = 'lpopt'; % Still MINOS is more reliable
       end
    elseif TV(9)
       Solver = 'CPLEX';
    elseif TV(8)
       Solver = 'xpress-mp';
    else
       Solver = 'qld';
    end
  case 'con'
    if LargeScale
        if TV(4)
           Solver = 'snopt';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'minos';
        else
           Solver = 'conSolve';
        end
    else
        if TV(3)
           Solver = 'npsol';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(2)
           Solver = 'minos';
        else
           Solver = 'conSolve';
        end
    end
  case 'uc'
    if LargeScale
        if TV(4)
           Solver = 'snopt';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        else
           Solver = 'ucSolve';
        end
    else
        if TV(3)
           Solver = 'npsol';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        else
           Solver = 'ucSolve';
        end
    end

  case {'ls','cls','exp'}
    if LargeScale
        if TV(4)
           Solver = 'slsSolve';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'slsSolve';
        elseif TV(7)
           Solver = 'filterSQP';
        else
           Solver = 'clsSolve';
        end
    else
        if TV(3)
           Solver = 'nlssol';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        else
           Solver = 'clsSolve';
        end
    end
  case {'lls'}
    if LargeScale
        if TV(4)
           Solver = 'sqopt';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(3)
           Solver = 'clsSolve'; % clsSolve instead of MINOS
           % Solver = 'minos';
        elseif TV(2)
           Solver = 'clsSolve'; % clsSolve instead of MINOS
           % Solver = 'minos';
        else
           Solver = 'clsSolve';
        end
    else
        if TV(3)
           Solver = 'lssol';
        elseif TV(1)
           Solver = 'lsei'; % Will be picked if not LSSOL present
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(4)
           Solver = 'sqopt';
        elseif TV(2)
           Solver = 'minos';
        else
           Solver = 'clsSolve';
        end
    end
  case {'mip'}
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(8)
       Solver = 'xpress-mp';
    else
       Solver = 'mipSolve';
    end
  case {'glb'}
    Solver = 'glbFast';
  case {'glc'}
    Solver = 'glcCluster';
  case {'miqp'}
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(8)
       Solver = 'xpress-mp';
    else
       Solver = 'miqpbb';
    end
  case {'minlp'}
    if TV(7)
       Solver = 'minlpbb';
    else
       Solver = 'glcCluster';
    end
  case {'sdp'}
    Solver = 'PENSDP';
  case {'bmi'}
    Solver = 'PENBMI';
  otherwise
    disp(Type)
    error('Illegal type of optimization problem')
end

% MODIFICATION LOG
% 010726 hkh Differentiate further between sparse and dense con / uc problems
% 020105 hkh Use glbFast and glcCluster as default
% 020701 hkh Use more general license handling, revise all selections
% 020701 hkh Add miqp, miqq, minlp and sdp types
% 030117 hkh Change type miqq to bmi
% 030211 hkh Change QP solver selection
% 030309 hkh Change LP,DLP,CON. First select SOL, otherwise /MINLP
% 041221 hkh Complete revision, prefer cplex to xpress-mp
