%
% function x=DefPar(s,f,xdef)
%
% If 'f' is a field in structure 's' and s.f is not empty, return it.
% If 'f' is not a field in s, or s.f is empty, return xdef.
%
%
function x=DefPar(s,f,xdef)

if(nargin < 3)
   xdef = [];
end

if isfield(s,f)
   x = getfield(s,f);
   if isempty(x), x = xdef; end
else
   x = xdef;
end
