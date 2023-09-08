%
% function J = lls_J(x, Prob)
%
% Compute Jacobian in Linear Least Squares Problem
%
% x      Parameter vector x 
% Prob   Problem structure
%        ||Cx-d|| is minimized where
%        C = Prob.LS.C
%        d = Prob.LS.y
%        The Jacobian is the matrix C
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Sep 11, 1999.      Last modified Nov 5, 2000.
%

function J = lls_J(x, Prob)

J = Prob.LS.C;

