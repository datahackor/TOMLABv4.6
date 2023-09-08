function [f,J]=evalJ(fun,x,Extra,m,JPI,verb,stepsize)
%
%evalJ	Compute Sparse Jacobian Via AD
%       evalJ is part of ADMIT toolbox. It computes sparse Jacobian
%       of a general nonlinear mapping.       
%
%       The user's function is identified by a unique identifier
%       fun. See ADMIT users manual for detail.
%
%       f=evalJ(fun,x)  computes just the function value at point x, the 
%                       is assume to return same number of output variables as in x.
%       [f,J]=evalJ(fun,x) Also computes the n x n sparse Jacobian at x.
%
%       Please refer to ADMIT users manual for complete reference. 
%       
%       JPI argument can be used to efficient execution. SEE GetJPI, dispJPI
%       Extra contains extra matrix argument for user's function. See users manual
%       for reference
%
%	ALSO SEE evalH
%
%
%       ******************************************************************
%       *              ADMIT Project   10-1-96                           *
%       *              Copyright (c) 1996 Cornell University.            *
%       ******************************************************************

nargin;

disp('DUMMY for TOMLAB if ADMAT/ADMIT not installed');
f=[];
J=[];

