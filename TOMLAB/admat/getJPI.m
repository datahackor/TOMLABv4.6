%       getJPI 
%       Compute sparsity and coloring information for sparse Jacobians.
%       This is really useful for efficiency of the  computing sparse
%       Jacobians with evalJ, because once computed sparsity and coloring
%       Information can be saved once for all. 
%
%       The user's function is identified by a unique identifier
%       fun. See ADMIT users manual for detail.
%
%       PartInfo = GetJPI(fun,m,n)  computes  sparsity + coloring Info for m xn J.
%
%       Please refer to ADMIT users manual for complete reference. 
%       
%       Extra contains extra matrix argument for user's function. See users manual
%       for reference
%
%       ALSO SEE evalJ, dispJPI
%
%
%       ******************************************************************
%       *              ADMIT Project   10-1-96                           *
%       *              Copyright (c) 1996 Cornell University.            *
%       ******************************************************************

function [PartInfo,SPJ]= getJPI(fun, m, n, Extra, method, SPJ)

nargin;

disp('DUMMY for TOMLAB if ADMAT/ADMIT not installed');
PartInfo=[];
SPJ=[];