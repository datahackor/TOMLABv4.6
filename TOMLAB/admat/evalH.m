function [val,grad,H]=evalH(fun,x,Extra,HPI,verb)
%
%evalH  Compute Sparse Hessians Via AD
%       evalH is part of ADMIT toolbox. It computes sparse Hessian
%       of a general scalar nonlinear mapping.        
%
%       The user's function is identified by a unique identifier
%       fun. See ADMIT users manual for detail.
%
%       val=evalH(fun,x)  computes just the function value at point x. 
%       [val,grad]=evalH(fun,x) Also computes the n x 1 gradient at x.
%       [val,grad,H]=evalH(fun,x) Also computes the n x n
%                    sparse Hessian matrix at x.
%
%       Please refer to ADMIT users manual for complete reference. 
%       
%       HPI argument can be used to efficient execution. SEE GetHPI,dispHPI 
%       Extra contains extra matrix argument for user's function. See users manual
%       for reference
%
%       ALSO SEE  evalJ
%
%
%       ******************************************************************
%       *              ADMIT Project   10-1-96                           *
%       *              Copyright (c) 1996 Cornell University.            *
%       ******************************************************************
 

nargin;
disp('DUMMY for TOMLAB if ADMAT/ADMIT not installed');
val=[];
grad=[];
H=[];

