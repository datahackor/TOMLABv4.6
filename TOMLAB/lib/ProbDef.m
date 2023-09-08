%
% Initialization of structure Prob
%
% function [Prob] = ProbDef(New);
%
% INPUT:
%
% New     Flag if a new type of problem is defined.
%         If New ~=0, initialize probType as empty
%         Also set counters: 
%              n_f n_g n_H n_c n_dc n_d2c n_J n_r n_d2r 
%         as double(0)
%
% OUTPUT:
%  Prob   Structure
%               
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written May 18, 1998.   Last modified Dec 13, 2004.
%

function [Prob] = ProbDef(New);

if nargin < 1
   New=[];
end

if isempty(New), New=0; end

if New
   %global solvType 
   %solvType = []; 
   global  probType
   probType = [];
   global n_f n_g n_H n_c n_dc n_d2c n_J n_r n_d2r 
   % Init counters to 0
   n_f=double(0); n_g=double(0); n_H=double(0); n_c=double(0); n_dc=double(0); 
   n_d2c=double(0); n_r=double(0); n_J=double(0); n_d2r=double(0);
end

LP = LineParamDef;
USER = struct('f',[], 'g',[], 'H',[], 'c',[], 'dc',[], 'd2c',[], ...
              'r',[], 'J',[], 'd2r',[],'fc',[],'gdc',[]);

SOL = struct( 'SpecsFile', [], 'PrintFile', [], 'SummFile', [], ...
    'xs', [], 'hs', [], 'nS', double(0), 'hElastic', [], ...
    'iState', [], 'cLamda', [], 'R', [], 'optPar', -999*ones(1,63), ...
    'optParN', 63); 


Prob = struct('Tomlab','v4.5','A', [], ...
  'ADObj', double(0), 'ADCons', double(0), 'BIG', [], 'b_L', [], 'b_U', [],...
  'c_L', [], 'c_U', [], 'CheckNaN', 0, 'cName', [], 'cols', [], 'ConIx', [],...
  'ConsDiff',double(0), 'ConsIdx', [], 'ConsPattern', [], ...
  'd2LPattern',[],'f_Low', -1E300, ...
  'f_opt', [], 'GradTolg', [], 'GradTolH', [], 'GradTolJ', [],...
  'HessIx',[],'HessPattern', [], 'JacIx', [],'JacPattern', [], ...
  'LargeScale', [], 'MaxCPU', inf ,'MENU', 0, ...
  'Mode', 2, 'nState',double(1), 'N', [], 'mLin',[],'mNonLin',[],'Name','',...
  'NumDiff', double(0), 'P', double(1),...
  'plotLine', 0, 'PriLev', double(0), 'PriLevOpt', double(0),  ...
  'probFile', [], 'probType', [], 'rows',[], 'simType', double(0), ...
  'smallA', double(0), ...
  'SolverDLP', [], 'SolverFP', [], 'SolverLP', [], 'SolverQP', [], ...
  'uP', [], 'uPName', [], 'WarmStart',0,'Warning',1, 'x_0', [], 'x_L', [],...
  'x_U', [], 'x_min', [], 'x_max', [], 'x_opt', [], 'xName', [], ...
  'QP', [], 'LS', [], 'MIP', [], 'GO', [], 'CGO', [], 'ExpFit', [], ...
  'NTS', [], 'LineParam', LP, 'optParam', [], 'PartSep', [], ...
  'SOL', SOL, 'Solver',[],'USER',USER,'DUNDEE',[],'PENOPT',[]);

% 'p_f', 'nlp_f', 'p_g', 'nlp_g', 'p_H', 'nlp_H', 'p_c', 'nlp_c', ...
% 'p_dc', 'nlp_dc', 'p_d2c', 'nlp_d2c', 'p_r', 'nlp_r', 'p_J', 'nlp_J', ...
% 'p_d2r', 'nlp_d2r', 'ExpFit', [], 'QP', [], 'MIP', [], 'LS', [], ...

Prob.LineParam = LineParamDef;

return

% Prob.MENU = 0;       % Flag used to tell if menu, GUI, or driver is calling
% Prob.uP   = [];      % User parameters
% Prob.uPName = [];    % Problem name which uP were defined for
%
%
% Prob.ExpFit  = [];    % Fitting of exponential sums
% Prob.QP      = [];    % Quadratic and Linear Programming
% Prob.MIP     = [];    % Mixed-Integer Programming
% Prob.LS      = [];    % Least squares problems
% Prob.NTS     = [];    % Nonlinear time series
% Prob.Solver  = [];    % Solver field. Solver.Alg; Solver.Name;Solver.Method
% Prob.PartSep = [];    % Partially separable functions

% Define the SOL parameters
% xs,hs,nS is used for warm start of SNOPT and SQOPT
% iState is used for warm start of QPOPT, NPSOL
% cLamda,R is used for warm start of NPSOL
% optParN: Number of parameters in optPar 
% SpecsFile: File to read SPECS parameters from
% PrintFile: Log Print File
% SummFile: Summary File (standard output, but not in Matlab)
% optPar: Vector with optimization parameters. 
% The meaning of each index is described in the m-file interface routines:
% minos.m (minosTL.m), snopt.m sqopt.m qpopt.m npsol.m etc.


% The User or Init routine should define the up to nine routine names
% The routines shall compute f, g, H, c, dc, d2c, r, J and d2r.
% c and dc are needed for constrained problems (and maybe d2c)
% r, J and d2r are needed for nonlinear least squares problems


% MODIFICATION LOG:
%
% 980825  hkh  Changed nts to NTS for Nonlinear Time Series.
% 980909  mbk  Structure field GLOBAL added.
% 980920  hkh  Change Prob.f_min to Prob.optParam.f_Low. Delete Prob.optPar
% 980921  mbk  optparam.f_Low changed to optParam.f_Low.
% 981006  hkh  Added name uPName, defined together with uP to avoid conflicts
% 981010  hkh  Changed function call logic, using TOMLAB gateway routines
% 981018  hkh  New field NLLS. t and Yt moved to NLLS field.
% 981020  hkh  New field PartSep, for partially separable functions
% 981023  hkh  Added Prob.p_d2r and Prob.USER.p_d2r
% 981026  hkh  Add field NumDiff for numerical differentiation
%              Changed back f_Low to top level of struct
% 981105  hkh  Delete field f_0
% 981111  hkh  Add definition of Solver field in Prob
% 981119  hkh  Add input variable New. 
% 990213  hkh  prob.x_min=[]; prob.x_max=[]; should be capital P 
% 000710  hkh  Add SOL field parameters
% 000928  hkh  Remove GLOBAL field
% 001014  hkh  Added field ConsDiff
% 001019  hkh  Always use fixed names for the TOMLAB gateway routines
% 001105  hkh  Make field NLLS into LS, add linear least squares
% 010903  hkh  Now 63 parameters in optPar
% 011110  hkh  Add two new fields, GO and RBF, used for global optimization
% 020103  hkh  Change field RBF to CGO
% 020409  hkh  Add field Mode = request for functions and gradients
% 020409  hkh  Add field nState,= 1 if 1st time call for functions & gradients
% 020512  hkh  Add field DUNDEE for Fletcher and Leyffer Dundee solvers
% 020630  hkh  Add field PENSDP for Tomlab /PENSDP
% 030107  hkh  Change field PENSDP to PENOPT, for /PENSDP and /PENBMI
% 030117  hkh  Add field CheckNaN, if NaN in derivatives should be checked
% 030129  hkh  Add field Tomlab with version number
% 030522  hkh  Add field fc and gdc in USER
% 030524  hkh  Add field simType, default 0, Use double(1), double(0)
% 031129  hkh  Add fields ADObj, ADCons, use double for nState and counters
% 031215  ango Change '4.0' to '4.2'
% 040102  hkh  Add fields mLin mNonLin for lengths of constraints
% 040106  hkh  Add fields ConIx, ConsIdx, and rearrange in alphabetic order
% 040106  hkh  Add fields JacIx, HessIx
% 040412  hkh  Add field Prob.Warning, put x_L and x_U after each other
% 040425  hkh  Add field Prob.smallA and Prob.MaxCPU
% 040608  hkh  Improve comments
% 040928  ango Change '4.2' to '4.4'
% 041201  med  Change '4.4' to '4.5'
% 041213  hkh  Add field BIG, default []
