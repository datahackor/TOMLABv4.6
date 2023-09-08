%
% function n = xnargout(S)
%
% Number of function output arguments
% Matlab 5.2 and later code
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.5.0$
% Written Aug 28, 2004.   Last modified Nov 17, 2004.

function n = xnargin(S)

if isa(S, 'function_handle')
 S = func2str(S);
end

% n = abs(nargout(S));

try
   n = abs(nargout(S));
catch
   err = sprintf('xnargout(''%s'') failed',S);
   error(err);
end

% MODIFICATION LOG
%
% 040828 ango Wrote file 
% 041117 ango try-catch construct is default.
