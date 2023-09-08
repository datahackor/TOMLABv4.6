%       GetHPI
%
%       Compute sparsity and coloring information for sparse Hessians.
%       This is really useful for efficiency of the  computing sparse
%       Hessians with evalH, because once computed sparsity and coloring
%       Information can be saved once for all.
%
%       The user's function is identified by a unique identifier
%       fun. See ADMIT users manual for detail.
%
%       PartInfo = GetHPI(fun,n)  computes  sparsity + coloring Info for nxn H.
%
%       Please refer to ADMIT users manual for complete reference.
%
%       Extra contains extra matrix argument for user's function. See users manual
%       for reference
%
%       ALSO SEE evalH, dispHPI
%
%
%       ******************************************************************
%       *              ADMIT Project   10-1-96                           *
%       *              Copyright (c) 1996 Cornell University.            *
%       ******************************************************************

function [PartInfo,SPH]= getHPI(fun, n,Extra, method, SPH)

nargin;

disp('DUMMY for TOMLAB if ADMAT/ADMIT not installed');
PartInfo=[];
SPH=[];