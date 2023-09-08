% Check if MAD TB is installed correctly or not
%
% function madOK = checkMAD(PriLev);
% 
% if PriLev == 1 print message about how to install AD properly
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.2.0$
% Written Dec 3, 2003.  Last modified Dec 3, 2003.

function madOK = checkMAD(PriLev);

if nargin < 1
   PriLev = [];
end

if isempty(PriLev), PriLev = 1; end

% Check MAD Automatic differentiation


global MADMinSparseN
if isempty(MADMinSparseN)
   % if exist('fmad') & exist('getvalue')
   if PriLev
      fprintf('\n================================================\n');
      fprintf('The MAD TB is not correctly installed.\n');
      fprintf('No automatic differentiation possible.\n');
      fprintf('Install the MAD Toolbox or ');
      fprintf('if installed, change tomlab/startup.m to get ');
      fprintf('the paths correct.\n');
      fprintf('And call madinitglobals to set any cleared globals.\n\n');
   end
   madOK = 0;
else
   madOK = 1;
end

% MODIFICATION LOG
%
% 031203  hkh  Written
% 040901  med  getvalue lower case