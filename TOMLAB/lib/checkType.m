%
% checkType returns true if the input type (probType) is the same
% as the wished type given as a string (strType)
%
% It may also return the string name of the type of optimization, 
% given the corresponding number, or the number given the string name of 
% the type of optimization 
%
% function isType = checkType(strType,numType)
%
% INPUT: 
% strType    The problem type given as a string, or empty
% numType    The problem type given as a number, or empty (or not given)
%
% OUTPUT: 
% isType  True if problem area defined by strType is the same as by numType
%         If isempty(strType), the string Type corresponding to numType(1) is 
%         given
%         If isempty(numType), the number corresponding to strType is given
% 
% Example:
%         isType = checkType('con',1);  gives isType == 0  (false)
%         isType = checkType('con',3);  gives isType == 1  (true)
%         isType = checkType([],3);     gives isType == 'con'
%         isType = checkType('con');    gives isType == 3
%
% Current types of optimization problems in TOMLAB:
%
%      1.  uc    Unconstrained Optimization
%      2.  qp    Quadratic Programming
%      3.  con   Constrained Optimization (nonlinear programming)
%      4.  ls    Nonlinear Least Squares
%      5.  lls   Linear Least Squares
%      6.  cls   Constrained Nonlinear Least Squares
%      7.  mip   Mixed-Integer Programming
%      8.  lp    Linear Programming
%      9.  glb   Box-bounded Global Optimization
%      10. glb   Constrained Global Optimization, also integer constraints
%      11. miqp  Mixed-Integer Quadratic Programming (MIQP)
%      12. minlp Mixed-Integer Nonlinear Programming (MINLP)
%      13. sdp   Semidefinite Programming, Linear SDP with LMI constraints
%      14. bmi   Linear SDP with BMI Constraints
%      15. exp   Parameter estimation in exponential models
%      16. nts   Nonlinear Time Series
%      17. lcp   Standard Linear Complementarity Problem (LCP)
%      18. mcp   Polyhedrally constrained variational inequality Problem or
%                Mixed Complementarity Problem(MCP)
%      19. miqq  Mixed-Integer Quadratic Programming with Quadratic 
%                constraints
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written Nov 5, 1998.    Last modified May 17, 2004.

function isType = checkType(strType,numType)

if nargin < 2
   numType = [];
end

if isempty(numType)
   if nargin == 2
      isType = 0;
      return
   end
   switch lower(strType)
      case 'uc' 
         isType=1;
      case 'qp' 
         isType=2;
      case 'con' 
         isType=3;
      case 'ls' 
         isType=4;
      case 'lls' 
         isType=5;
      case 'cls' 
         isType=6;
      case 'mip' 
         isType=7;
      case 'lp' 
         isType=8;
      case 'glb' 
         isType=9;
      case 'glc' 
         isType=10;
      case 'miqp' 
         isType=11;
      case 'minlp' 
         isType=12;
      case 'sdp' 
         isType=13;
      case 'bmi' 
         isType=14;
      case 'exp' 
         isType=15;
      case 'nts' 
         isType=16;
      case 'lcp' 	
         isType=17;
      case 'mcp' 	
         isType=18;
      case 'miqq' 	
         isType=19;   
      otherwise
         isType=[];
   end
   return
end

for i = 1:length(numType)
   switch numType(i)
      case 1 
         Type='uc';
      case 2
         Type='qp';
      case 3
         Type='con';
      case 4
         Type='ls';
      case 5
         Type='lls';
      case 6
         Type='cls';
      case 7
         Type='mip';
      case 8
         Type='lp';
      case 9 
         Type='glb';
      case 10
         Type='glc';
      case 11
         Type='miqp';
      case 12
         Type='minlp';
      case 13
         Type='sdp';
      case 14
         Type='bmi';
      case 15
         Type='exp';
      case 16
         Type='nts';
      case 17 	
         Type='lcp';
      case 18 	
         Type='mcp';
      case 19 	
         Type='miqq';   
      otherwise
         Type=' ';
   end

   if isempty(strType)
      isType = Type;
      return
   end

   isType = strcmpi(Type,strType);
 
   if isType, return; end
end

% MODIFICATION LOG
%
% 001105 hkh Written
% 020701 hkh Adding four new types, miqp, minlp, sdp, miqq.
% 021010 hkh Set isType=0 if probType is empty (may occur in tomGUI)
% 030117 hkh Change miqpp to bmi
% 040419 med Added lcp and mcp
% 040517 med Added miqq