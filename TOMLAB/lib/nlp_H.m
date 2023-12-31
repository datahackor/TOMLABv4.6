%
% function H = nlp_H(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the Hessian matrix
%
% H(i,j) is d2f/dx_i,dx_j
% Currently the symmetry is not taken into account
%
% nlp_H calls the routine Prob.USER.H either as 
%           H=feval(Prob.USER.H, x) or
%           H=feval(Prob.USER.H, x, Prob) or
%           H=feval(Prob.USER.H, x, Prob, varargin{:})
% depending on the number of inputs
%
% The global counter variable n_H is incremented
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.1$
% Written Oct 10, 1998.   Last modified May 26, 2004.
%

function H = nlp_H(x, Prob, varargin)

global n_H NARG PartSep

% Communication nlp_g/H 
global NLP_g NLP_xg NLP_pSepIndex

n_H=n_H+1;

Func = Prob.USER.H;

x=x(:);
NumDiff =  abs(Prob.NumDiff);

if Prob.ADObj == -1
   global mad_g
   if all(NLP_xg == x)
      H=getinternalderivs(mad_g);
   else
      nlp_g(x, Prob, varargin);
      H=getinternalderivs(mad_g);
   end
elseif Prob.ADObj == -2
   global HPI ad_f ad_g ad_H
   H=full(ad_H);
elseif isempty(Func) | NumDiff > 0
   % Here we call a numerical difference routine
   % Check if we have computed g(x) already (the normal case)
   gx=[];
   if PartSep
      if ~isempty(NLP_xg)
         if NLP_pSepIndex == Prob.PartSep.index
            if all(x == NLP_xg) 
               gx = NLP_g;
            end
         end
      end
   elseif ~isempty(NLP_xg)
      if all(x == NLP_xg) 
         gx = NLP_g;
      end
   end
   if NumDiff > 1 & NumDiff < 5
      H=FDHess2(x, Prob, gx, varargin{:}); % Using splines
      %if exist('csapi') & exist('csaps') & exist('spaps')
      %   H=FDHess2(x, Prob, gx, varargin{:}); % Using splines
      %else
      %   fprintf('Can not find directory SPLINES\n')
      %   fprintf('Have you got a license for Spline Toolbox?\n')
      %   fprintf('Using finite differece routine FDHess.m\n')
      %   H=FDHess(x, Prob, gx, varargin{:}); % FD Algorithm
      %end
   else
      H = FDHess(x, Prob, gx, varargin{:}); % FD Algorithm
   end   
else
   N = min(length(x),Prob.N);
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(3);
   end
   if p>2
      H=feval(Func, x(1:N), Prob, varargin{:});
   elseif p==2
      H=feval(Func, x(1:N), Prob);
   else
      H=feval(Func, x(1:N));
   end
   if Prob.CheckNaN ~= 0
      [iN,jN,HNew] = find(isnan(H));
      if ~isempty(iN)    % There are elements set to NaN, to be estimated
         % Set HessPattern to the NaN elements
         Prob.HessPattern = sparse(iN,jN,HNew,size(H,1),size(H,2));
         Prob.HessIx = findpatt(Prob.HessPattern);
         gx=[];
         if PartSep
            if ~isempty(NLP_xg)
               if NLP_pSepIndex == Prob.PartSep.index
                  if all(x == NLP_xg) 
                     gx = NLP_g;
                  end
               end
            end
         elseif ~isempty(NLP_xg)
            if all(x == NLP_xg) 
               gx = NLP_g;
            end
         end
         if (abs(Prob.CheckNaN) > 1) & (abs(Prob.CheckNaN) < 5)
            HNew=FDHess2(x(1:N), Prob, gx, varargin{:}); % Using splines
            %if exist('csapi') & exist('csaps') & exist('spaps')
            %   HNew=FDHess2(x, Prob, gx, varargin{:}); % Using splines
            %else
            %   fprintf('Can not find directory SPLINES\n')
            %   fprintf('Have you got a license for Spline Toolbox?\n')
            %   fprintf('Using finite differece routine FDHess.m\n')
            %   HNew=FDHess(x, Prob, gx, varargin{:}); % FD Algorithm
            %end
         else
            HNew = FDHess(x(1:N), Prob, gx, varargin{:}); % FD Algorithm
         end   
         % Merge analytic and numerical dc
         H(isnan(H)) = HNew(isnan(H));
      end
   end
end
%xprint(x,'x:')
%H

% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation 
% 981124   mbk   Added calls to FDHess and FDHess2 for FD and
%                spline approximation ov the Hesian
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
% 981127   hkh   Safeguard test on field PartSep.
% 990615   hkh   global NLP_pSepFunc NLP_pSepIndex was not defined
% 990909   hkh   Avoid structural information
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020528   hkh   Efficient handling of partially separable functions
% 030127   hkh   Check NaN elements if Prob.CheckNaN ~= 0, estimate numerically
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 040407   hkh   Always call FDHess2, no check on spline TB routines
% 040526   hkh   Use x(1:N) in function calls

