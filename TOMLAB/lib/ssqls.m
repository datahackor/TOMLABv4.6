% Solution of sparse linear least squares problems using subspace minimization
%
% x=ssqls(A,b,epsRank) solves the sparse linear least squares problem
%                    min ||Ax-b||
%                     x          2
% using QR factorization and application of Q to the right-hand side b.
%
% epsRank is used to determine the pseudo rank pRank
% ------------------------------------------------------------------------
%
% INPUT PARAMETERS
% A        The matrix A
% b        The right hand side
% epsRank  Rank tolerance. Default 1E-12
% 
% OUTPUT PARAMETERS
%
% x        Solution
% pRank    Pseudo rank
% maxR     Maximal diagonal element in R
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Apr 16, 2002. Last modified Apr 16, 2002.
%

% ------------------------------------------------------------------------
function [x,pRank,maxR] = ssqls(A,b,epsRank)

if nargin < 3
   epsRank=1E-12;
end

[m n]= size(A);

% Minimum-degree ordering of A.

q = colmmd(A);

[R,p,c]=sqr2(A(:,q),b);

% Resulting column permutation.

Pc = q(p);

if n > 1 & m > 1
   maxR  = max(abs(diag(R)));
   pRank = full(nnz(abs(diag(R)) >= maxR*epsRank));
else
   maxR  = abs(R(1,1));
   pRank = 1;
end

%if ~isempty(R), maxR=abs(R(1,1)); else, maxR=0; end
%if n > 1 & m > 1
%   pRank = nnz(abs(diag(R)) >= maxR*epsRank);
%else
%   pRank = 1;
%end

% Compute the least squares solution.

%z1 = tomsol(3,full(R(1:pRank,1:pRank)),c(1:pRank));
%z2 = tomsol(3,R(1:pRank,1:pRank),c(1:pRank));


%diffZ=norm(z-z1)
%diffZ2=norm(z-z2)

x=R\c;

%MaxRefine = 0;
%if MaxRefine == 0
%   x=R\c;
%else
%   % One step of iterative refinement. Not coded yet
%
%   x = zeros(n,1);
%
%   %pRank
%   %n
%
%   x(1:pRank) = ItRef(A(:,q), b, z, R(1:pRank,1:pRank), pRank, 1);
%
%end

% Permute the computed solutions

x(Pc) = x;

%if maxR == 0
%   pRank=0;
%end

function [x] = ItRef(A, b, x, R, pRank, MaxIter)

for i=1:MaxIter
    r   = b-A(:,1:pRank)*x(1:pRank);
 rNorm = norm(r)
    %ATb = A'*r;
    %y   = R'\ATb;
    y   = tomsol(2,R',A(:,1:pRank)'*r);
    %dx  = R\y;
    dx  = tomsol(3,R,y);
    x(1:pRank)=x(1:pRank)+dx;
end

 r   = b-A(:,1:pRank)*x(1:pRank);
 NewNorm = norm(r)

