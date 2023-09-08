%		L2_H.m
%
% function H=L2_H(x, Prob, varargin)
%
% L2_H computes the 2nd derivative of the rewritten sparse least squares 
% objective function f in the point x
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Apr 13, 2002. Last modified Apr 13, 2002.

function H=L2_H(x, Prob, varargin)

m = Prob.L2.m;

n = length(x)-m;

H = [spalloc(n,n+m,0);spalloc(m,n,0),speye(m,m)];

%
% 020413 hkh Written

