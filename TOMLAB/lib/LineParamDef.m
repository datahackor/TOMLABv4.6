%
% Initialization of structure LineParam
%
% LineParam contains the line search parameters
%
% function LineParam = LineParamDef;
%
% INPUT:
%
% OUTPUT:
%   LineParam  Structure
%               
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.4.0$
% Written Sept 24, 1998. Last modified Oct 18, 2004.
%

function LineParam = LineParamDef;

LineParam = struct( 'LineAlg',1, 'sigma',[], 'InitStepLength',1, ...
      'MaxIter',15, 'fLowBnd',-realmax, 'rho',0.01, ...
      'tau1',9, 'tau2',0.1, 'tau3',0.5, 'eps1',1E-7, 'eps2',100*eps );

return

% Some comments about the optimization parameters

% (x) refers to the OPT TBX 1.x options vector

% (7) LineAlg: Line search algorithm. 0 = quadratic, 1 = cubic, 3 = curved 
% (if available)

% (18) Initial step length. (Default 1 or less). 

% (21) LineSearch.sigma: Line search accuracy; 0 < sigma < 1
% sigma=0.9 inaccurate line search. sigma 0.1 accurate line search



% MODIFICATION LOG:
%
% 000923  hkh  Written. Extracted from optParamDef.
% 041018  hkh  Set default sigma as [], not all solvers have 0.9 as default

