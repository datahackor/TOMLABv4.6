% CPLEX MEX-interface internal callback routine
%
% Makes the cpxRetVec vector global
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 9.0.0$
% Written Aug. 8, 2002 Last modified Sept 22, 2002
%

function cpx2retvec(a)

global cpxRetVec
cpxRetVec = a;
