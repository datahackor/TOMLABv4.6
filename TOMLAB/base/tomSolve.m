% function Result = tomSolve(Solver,Prob)
% 
% tomSolve solves a sub problem without any check on the Prob structure 
%
% It is intended for use by other solvers when solving a subproblem like a
% QP, LP, Dual LP or FP (feasible point) problem
%
% Global variables are saved before the solver call, and restored afterwards
%
% It could also be used for recursive solutions of user problems
% or control loops when speed is a demand.
%
% INPUT PARAMETERS
% Solver    Name of the solver
% Prob      Input structure, feeded to the solver
%
% OUTPUT PARAMETERS
% Result    Output result structure, feeded back to the caller.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Nov 6, 2000.    Last modified Jan 12, 2005.
%

function Result = tomSolve(Solver,Prob)

if nargin < 2
   error('SolveQP needs input structure Prob');
end

%global solvType
%solvTypeSave=solvType;

% Set solvType to the subproblem type, when solving this subproblem
% solvType=Prob.solvType;

global GlobalLevel

if isempty(GlobalLevel)
   GlobalLevel=1;
else
   GlobalLevel=GlobalLevel+1;
end
Level=GlobalLevel;
globalSave(Level);

switch lower(Solver)
   case {'qpopt'}
      Result = qpoptTL(Prob);
   case {'lpopt'}
      Result = lpoptTL(Prob);
   case {'lp-minos'}
      Result = minoslpTL(Prob);
   case {'qld'}
      Result = qldTL(Prob);
   case {'qp-minos'}
      Result = minosqpTL(Prob);
   case {'minos'}
      Result = minosTL(Prob);
   case {'qpsolve'}
      Result = qpSolve(Prob);
   case {'lpsolve'}
      Result = lpSolve(Prob);
   case {'dualsolve'}
      Result = DualSolve(Prob);
   case {'lsei'}
      Result = lseiTL(Prob);
   case {'tlsqr'}
      Result = TlsqrTL(Prob);
   case {'lssol'}
      Result = lssolTL(Prob);
   case {'sqopt'}
      Result = sqoptTL(Prob);
   case {'npsol'}
      Result = npsolTL(Prob);
   case {'nlssol'}
      Result = nlssolTL(Prob);
   case {'snopt'}
      Result = snoptTL(Prob);
   case {'ucsolve'}
      Result = ucSolve(Prob);
   case {'consolve'}
      Result = conSolve(Prob);
   case {'nlpsolve'}
      Result = nlpSolve(Prob);
   case {'clssolve'}
      Result = clsSolve(Prob);
   case {'glbsolve'}
      Result = glbSolve(Prob);
   case {'glbfast'}
      Result = glbFast(Prob);
   case {'glcsolve'}
      Result = glcSolve(Prob);
   case {'glcfast'}
      Result = glcFast(Prob);
   case {'glccluster'}
      Result = glcCluster(Prob);
   case {'rbfSolve'}
      Result = rbfSolve(Prob);
   case {'mipsolve'}
      Result = mipSolve(Prob);
   case {'strustr'}
      Result = sTrustr(Prob);
   case {'bqpd'}
      Result = bqpdTL(Prob);
   case {'filtersqp'}
      Result = filterSQPTL(Prob);
   case {'pensdp'}
      Result = pensdpTL(Prob);
   case {'penbmi'}
      Result = penbmiTL(Prob);
   case {'nlpql','nlpqlp'}
      Result = nlpqlTL(Prob);

   case {'pdco'}
      Result = pdcoTL(Prob);
   case {'pdsco'}
      Result = pdscoTL(Prob);
      
   case {'oqnlp'}
      Result = oqnlpTL(Prob);
   case {'msnlp'}
      Result = msnlpTL(Prob);   
   case {'msnlp'}
      Result = lsgrg2TL(Prob);
      
   case {'knitro'}
      Result = knitroTL(Prob);
   case {'conopt'}
      Result = conoptTL(Prob);

   case {'miqpbb'}
      Result = miqpBBTL(Prob);
   case {'minlpbb'}
      Result = minlpBBTL(Prob);
   case {'cutPlane'}
      Result = cutplane(Prob);
   case {'xpress-mp'}
      Result = xpressTL(Prob);
   case {'cplex'}
      Result = cplexTL(Prob);
   case {'xa'}
      Result = xaTL(Prob);

   case {'fmincon'}
      Result = opt20Run('fmincon',Prob);
   case {'fminsearch'}
      Result = opt20Run('fminsearch',Prob);
   case {'fminunc'}
      Result = opt20Run('fminunc',Prob);
   case {'lsqcurvefit'}
      Result = opt20Run('lsqcurvefit',Prob);
%  case {'lsqlin'}
%     Result = opt20Run('lsqlin',Prob);
   case {'lsqnonlin'}
      Result = opt20Run('lsqnonlin',Prob);
   case {'lsqnonneg'}
      Result = opt20Run('lsqnonneg',Prob);
   case {'constr'}
      Result = opt15Run('constr',Prob);
   case {'fmins'}
      Result = opt15Run('fmins',Prob);
   case {'fminu'}
      Result = opt15Run('fminu',Prob);
   case {'leastsq'}
      Result = opt15Run('leastsq',Prob);
   case {'ego'}
      Result = ego(Prob);
   %case {'infSolve'}
   case {'quadprog'}
      Result = opt20Run('quadprog',Prob);
   case {'linprog'}
      Result = opt20Run('linprog',Prob);
   case {'qp'}
      Result = opt15Run('qp',Prob);
   case {'lp'}
      Result = opt15Run('lp',Prob);
   otherwise
      Result = tomRun(Solver,Prob);
end

globalGet(Level);

GlobalLevel=GlobalLevel-1;

%solvType=solvTypeSave;

% MODIFICATION LOG:
%
% 001106  hkh  Written
% 010717  hkh  Added glbFast and Xpress-MP
% 010815  hkh  Added glcFast
% 011104  hkh  Added glcCluster
% 020526  hkh  Adding Dundee QP solver bqpd
% 020621  hkh  Adding Dundee MIQP solver MIQPbb
% 020630  hkh  Adding Dundee solvers MINLPbb, filterSQP; and PENSDP
% 020702  hkh  Adding CPLEX
% 030116  hkh  Change to Tlsqr, add PENBMI
% 030123  hkh  Add pdco, pdsco
% 030129  hkh  Remove empty varargin in call to opt15Run and opt20Run
% 040101  hkh  Add nlpql
% 040103  hkh  Direct call to all TL files
% 040312  ango Spelling of qpoptTL, lpoptTL fixed
% 041216  med  OQNLP and MSNLP added
% 050112  med  Added LSGRG2, KNITRO and CONOPT, XA