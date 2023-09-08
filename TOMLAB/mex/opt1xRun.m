% opt1xRun: Code to run optimization solvers in MathWorks Optimization TB 1.5
%
% function Result = opt1xRun(Solver, Prob, varargin)
%
% Solver is currently one of:
%    CONSTR
%    FMINU
%    LP
%    QP
%    LEASTSQ
%
%    FMINS is outside opttb 1.x
%
% The test on the name in Solver is not case-sensitive
%
% This is a dummy routine when no Optimization Toolbox routines are present
%
% Prob   is the TOMLAB problem structure
%
% Result is the TOMLAB result structure
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written Sep 16, 1999.  Last modified April 14, 2004.
%

function Result = opt1xRun(Solver, Prob, varargin)

fprintf('\n\n\nopt1xRun: TOMLAB Dummy. ');
fprintf('Optimization Toolbox routines not available');
fprintf(' or not in PATH\n\n\n');

fprintf('If you have Optimization Toolbox 1.x see startup.m\n');
fprintf('The directory tomlab\\optim1.x must be placed before ');
fprintf('tomlab\\mex in the path\n\n\n');

Result = [];