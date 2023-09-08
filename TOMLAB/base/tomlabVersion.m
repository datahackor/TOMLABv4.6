% function [tomV,OS] = tomlabVersion
%
% tomV Tomlab Version
% OS   Operating system
% 
% Returns the current Tomlab version
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., $Release: 4.0.7$
% Written Oct 1, 2000.     Last modified Jul 28, 2003.

function [tomV,OS,TV] = tomlabVersion

% Tomlab v3.1: Version any of MINI, v3.0, SOL, NPSOL, SNOPT, CGO
%
% Tomlab v3.2 and v4.0: 
% Version any of Base Module, MINOS, NPSOL, SNOPT, SOL, CGO, PENSDP, MINLP,
%                Xpress, CPLEX, PENBMI
%
% Values =1 in TV vector if license for module is OK. Indices described below:
%
% TV   Tomlab module
% 1.   Base Module 
% 2.   MINOS 
% 3.   NPSOL 
% 4.   SNOPT 
% 3+4. SOL 
% 5.   CGO 
% 6.   PENSDP 
% 7.   MINLP
% 8.   Xpress 
% 9.   CPLEX
% 10.  PENBMI
% 11.  KNITRO
% 12.  CONOPT
% 13.  AMPL
% 14.  OQNLP
% 15.  XA
% 16.  NLPQL
% 17.  LGO
% 18.  MAD
% 19.  DIDO

% This is the command in v3.0:
%[x1,x2,x3,x4,x5,x6,tomV,OS]=tomlab(1);

[x1,x2,x3,x4,x5,x6,tomV,OS]=tomlablic(1);

NoOfProducts = 40;

if isstr(tomV)
   TV    = [1; zeros(NoOfProducts,1)];
   TV(1) = 1;
   if findstr('v3',tomV)
      TV(2) = 1;
   elseif findstr('SO',tomV)
      TV(2:4) = 1;
   elseif findstr('SN',tomV)
      TV(2) = 1;
      TV(4) = 1;
   elseif findstr('NP',tomV)
      TV(2) = 1;
      TV(3) = 1;
   end
   if findstr('CGO',tomV)
      TV(5) = 1;
   end
else
   TV   = [tomV(2:end);zeros(NoOfProducts-length(tomV),1)];
   tomV = tomlablic(2);
end

