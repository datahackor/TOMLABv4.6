% daceInit.m
%
% daceInit finds an initial grid of "space-filling" design.
%
% function X = daceInit(k, x_L, x_U, M);
%
% k is in principle the dimension. It determines how many points to use
%   k    # of points
%   1    21
%   2    21
%   3    33
%   4    41
%   5    51
%   6    65
%  >6    65
%
% x_L       Lower bounds for each element in x.
% x_U       Upper bounds for each element in x.
%
% M         Number of points to generate, if to overrule the k value
%           Default empty, i.e. M is not used, only k 
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Oct 29, 1998.   Last modified March 9, 2005.

function X = daceInit(k, x_L, x_U, M)

if nargin < 4, M=0; end

m = length(x_L);

if M > 0, k = 0; end

if k==2  & m==2
   n = 21;           % number of sampled points
   %u = 0:0.05:1;
   v = [0.5 0.15 0.75 1 0.35 0.6 0.1 0.85 0.3 0.55 0.05 ...
        0.95 0.7 0.25 0.45 0 0.9 0.65 0.2 0.4 0.8];
   X = zeros(2,n);
   y = zeros(n,1);
   X(1,:) = x_L(1) + ( x_U(1) - x_L(1) )*[0:0.05:1];
   X(2,:) = x_L(1) + ( x_U(1) - x_L(1) )*v;
else
   if M > 0
      n = M;
   elseif k <= 2
      n = 21;
   elseif k == 3
      n = 33;
   elseif k == 4
      n = 41;
   elseif k == 5
      n = 51;
   else
      n = 65;
   end
   
   U = linspace(0,1,n);

   %rand('state',2);
   %rand('state',122);
   %rand('state',1313); % Hart6_2.log, 65 initiala punkter
   %rand('state',13); % Hart6_3.log, 81 initiala punkter
   %rand('state',112); % Hart6_4.log, 65 initiala punkter

   for i = 2:m
       v = linspace(0,1,n);
       dummy = rand(1,n);
       [a b]=sort(dummy);
       v = v(b);
       U = [U;v];
   end
   
   X = zeros(m,n);
   for i = 1:n
       X(:,i) = x_L + ( x_U - x_L).*U(:,i);
   end
end

if 0%k == 2  % Plot initial points
   plot(X(1,:),X(2,:),'*r');
   pause
end

% MODIFICATION LOG:
%
% 980626  hkh  Avoid feval
% 010715  hkh  Only compute X, use lower and upper bounds as input
% 020831  hkh  Use rand('state'), Matlab 5.x random generator
% 040306  hkh  Avoid setting the state, done in ego and rbfSolve 
% 040307  hkh  Add extra parameter M, to generate a large set of points
% 050117  med  mlint revision