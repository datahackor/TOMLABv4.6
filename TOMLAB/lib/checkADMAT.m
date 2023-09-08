% Check if ADMAT/ADMIT-1 TB is installed correctly or not
%
% function admatOK = checkADMAT(PriLev);
% 
% if PriLev == 1 print message about how to install AD properly
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.2.0$
% Written Dec 3, 2003.  Last modified Dec 6, 2003.

function admatOK = checkADMAT(PriLev);

if nargin < 1
   PriLev = [];
end

if isempty(PriLev), PriLev = 1; end

% Check ADMAT Automatic differentiation
global ADhess
if isempty(ADhess)
   % if ~exist('getHPI') | ~exist('evalH')
   if PriLev
      fprintf('\n================================================\n');
      fprintf('The ADMAT/ADMIT-1 TB is not correctly installed.\n');
      fprintf('No automatic differentiation is possible.\n');
      fprintf('Get the ADMIT-1 Toolbox from:\n');
      fprintf('http://www.cs.cornell.edu/home/verma/AD/research.html\n');
      fprintf('or if installed, change tomlab/startup.m to get ');
      fprintf('the paths correct.\n');
      fprintf('And call admatinitglobals to set any cleared globals.\n\n');
   end
   admatOK = 0;
else
   admatOK = 1;
end

% MODIFICATION LOG
%
% 031203  hkh  Written
% 031206  hkh  admatinitglobals to be called if globals missing
% 040510  med  Changed name for function
