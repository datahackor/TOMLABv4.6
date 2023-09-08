% probInit.m :
%
% General initialization routine for TOMLAB optimization problems, 
% when a file in the TOMLAB Init File format is used
%
% function [Prob] = probInit (probFile, P, ask, Prob);
%
% INPUT:  
% probFile Name of problem definition Init File, as a string, e.g. 'con_prob'.
% P        Problem number
%          If P=0, call Init File to get the defined problems and let the user
%          choose a problem by displaying a menu
%          if isempty(P), just return the number of problems
% ask      1:  ask questions; 0: use defaults; -1: use values in prob.uP
%          11: ask questions in the GUI window
%          If isempty(ask) then [if isempty(prob.uP),ask=1, else ask=-1];
% Prob     Structure with all parameters defining the problem
%
% OUTPUT: 
% Prob     Problem structure or 
%          Number of problems, if isempty(P)
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written June 16, 1999. Last modified Oct 12, 2000.
%

function [Prob] = probInit (probFile, P, ask, Prob);

if nargin < 4
   Prob=[];
   if nargin < 3
      ask=[];
      if nargin < 2
         P=[];
         if nargin < 1
            probFile='con_prob';
         end
      end
   end
end

if isempty(P)   % Only return the number of predefined problems
   probList=feval(probFile);
   nProbs=size(probList,1);
   Prob=nProbs;
   return
end

if isempty(ask) 
   if isstruct(Prob)
      if isempty(Prob.uP) 
         ask=1; 
      else
         ask=-1; 
      end
   else
      ask=-1; 
   end
end

if P(1)<=0  % Give menu to choose problem from
   P=strmenu('Choice of test function ',feval(probFile,[]));
end


[probList, Prob] = feval(probFile, P, ask, Prob);

global probType

probType=Prob.probType;

Prob.probFile=probFile;

% MODIFICATION LOG:
%
% 980825  hkh  Changed error in comments. Changed call to usr_prob-routine.
% 980910  mbk  Do not set Prob to [] on line 72.
% 981026  hkh  Changed ProbFile to probFile. Set field Prob.probFile to probFile
% 990617  hkh  Change comments
% 990622  hkh  Add first argument (empty) to initial call to probFile,usr_prob
% 001012  hkh  Delete usr_prob possibility, not needed any longer

