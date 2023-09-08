%
% Set values in structure LineParam
%
% The fields in LineParam store line search parameters
%
% function LineParam = LineParamSet(LineParam);
%
% INPUT:
%  LineParam   Structure
%
% OUTPUT:
%  LineParam   Structure
%               
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Sept 24, 2000. Last modified Sept 24, 2000.
%

function LineParam = LineParamSet(LineParam);

lP=LineParamDef;

if ~isstruct(LineParam)
   LineParam = lP;
else
   if ~isfield(LineParam,'LineAlg')
      LineParam.LineAlg=lP.LineAlg;
   end
   if ~isfield(LineParam,'sigma')
      LineParam.sigma=lP.sigma;
   end
   if ~isfield(LineParam,'InitStepLength')
      LineParam.InitStepLength=lP.InitStepLength;
   end
   if ~isfield(LineParam,'MaxIter')
      LineParam.MaxIter=lP.MaxIter;
   end
   if ~isfield(LineParam,'fLowBnd')
      LineParam.fLowBnd=lP.fLowBnd;
   end
   if ~isfield(LineParam,'rho')
      LineParam.rho=lP.rho;
   end
   if ~isfield(LineParam,'tau1')
      LineParam.tau1=lP.tau1;
   end
   if ~isfield(LineParam,'tau2')
      LineParam.tau2=lP.tau2;
   end
   if ~isfield(LineParam,'tau3')
      LineParam.tau3=lP.tau3;
   end
   if ~isfield(LineParam,'eps1')
      LineParam.eps1=lP.eps1;
   end
   if ~isfield(LineParam,'eps2')
      LineParam.eps2=lP.eps2;
   end
end

% MODIFICATION LOG
%
% 000924  hkh  Written from OptParamSet.
