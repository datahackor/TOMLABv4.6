%
% function n = xnargin(S)
%
% Number of function input arguments
% Matlab 5.2 and later code
%
% A bug in Matlab 5.1 on PC windows system makes it necessary to run
% different code on Matlab 5.1 and earlier versions.
%
% On 5.1 we must use lower on the argument to nargin.
%
% This routine does not check if S exist, xxnargin does.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.3.1$
% Written Nov 26, 1998.    Last modified Jan 31, 2003.

function n = xnargin(S)

if isa(S, 'function_handle')
 S = func2str(S);
end

% n = abs(nargin(S));

% --- STANDALONE (mcc) USERS, SEE HERE. ---
%
% The following try-catch code can sometimes be useful instead of the 
% original n = abs(nargin(S)) above:

try
   n = abs(nargin(S));
catch
   err = sprintf('xnargin(%s) failed, please check argument %s\n',S,S);
   error(err)
end

% --- END STANDALONE CODE ---

% MODIFICATION LOG
%
% 030131 hkh  Add function_handle treatment
% 040728 ango Alternative version useful for mcc users added, commented. 
% 041112 frhe This version is a copy of xnargin.m, but with the mcc
%             change applied.