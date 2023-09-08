%
% function J = nlp_J(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the Jacobian matrix
%
% nlp_J calls the routine Prob.USER.J either as 
%           J=feval(Prob.USER.J, x) or
%           J=feval(Prob.USER.J, x, Prob) or
%           J=feval(Prob.USER.J, x, Prob, varargin{:})
% depending on the number of inputs in Prob.USER.J
%
% J(i,j) is J/dx_i
%
% The global counter variable n_J is incremented
%
% Weighting is computed using the global wLS (computed in nlp_r.
%
% J and x are stored in globals LS_J and LS_xJ, to avoid recomputing.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.1$
% Written Oct 10, 1998.  Last modified May 26, 2004.
%

function J = nlp_J(x, Prob, varargin)

global wLS
global n_J NARG
global LS_x LS_r LS_xJ LS_J

x=x(:);
N=Prob.N; % N is the number of variables to be sent to user routines

if ~isempty(LS_xJ)
   if Prob.LS.SepAlg
      % LS_xJ or x may have extra variables, if separable NNLS, therefor (1:n)
      n=min(length(x),length(LS_xJ));
   else
      n=length(x);
   end
   if all(x(1:n)==LS_xJ(1:n))
      J=LS_J;
      return
   end
end

n_J=n_J+1; 

Func = Prob.USER.J;

if Prob.ADObj == 1
   global mad_r
   J=getinternalderivs(mad_r);
elseif Prob.ADObj == -1  
   global mad_J
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(8);
   end
   if isa(x,'fmad')
      x=getvalue(x);
   end
   if p > 2
      mad_J=feval(Func, fmad(x(1:N),speye(N)), Prob, varargin{:});
   elseif p > 1
      mad_J=feval(Func, fmad(x(1:N),speye(N)), Prob);
   else
      mad_J=feval(Func, fmad(x(1:N),speye(N)));
   end
   J=getvalue(mad_J);
elseif Prob.ADObj == 2
   global JPI ad_r ad_J ad_rM
   J=full(ad_J);
elseif isempty(Func) | Prob.NumDiff > 0
   if isempty(Prob.USER.r)
      LS_x=[];
      LS_r=[];
      LS_J=[];
      J=[];
      return
   else
      % Numerical difference computation of J
      if ~isempty(LS_x)
         if all(x==LS_x)
            rx=LS_r;
         else
            rx=[];
         end
      else
         rx=[];
      end
      if any(Prob.NumDiff == [2 3 4])
         [J,rx]=FDJac2(x(1:N), Prob, 'nlp_r', rx, varargin{:} ); % Using splines
         %if exist('csapi') & exist('csaps') & exist('spaps')
         %   [J,rx]=FDJac2(x, Prob, 'nlp_r', rx, varargin{:} ); % Using splines
         %else
         %   fprintf('Can not find directory SPLINES\n')
         %   fprintf('Have you got a license for Spline Toolbox?\n')
         %%   fprintf('Using finite differece routine FDJac.m\n')
         %   [J,rx]=FDJac(x, Prob, 'nlp_r', rx, varargin{:} ); % FD Algorithm
         %end
      elseif Prob.NumDiff == 5
         % Using complex number technique
         [J,rx]=FDJac3(x(1:N), Prob, 'nlp_r', rx, varargin{:} ); 
      else
         [J,rx]=FDJac(x(1:N), Prob, 'nlp_r', rx, varargin{:} ); % FD Algorithm
      end   
      LS_x = x;
      LS_r = rx;
   end
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(8);
   end

   if p > 2
      J=feval(Func, x(1:N), Prob, varargin{:} );
   elseif p==2
      J=feval(Func, x(1:N), Prob);
   else
      J=feval(Func, x(1:N));
   end
   if Prob.CheckNaN > 0
      [iN,jN,J2] = find(isnan(J));
      if ~isempty(iN)    % There are elements set to NaN, to be estimated
         Prob.JacPattern = sparse(iN,jN,J2,size(J,1),N);
         if ~isempty(LS_x)
            if all(x==LS_x)
               rx=LS_r;
            else
               rx=[];
            end
         else
            rx=[];
         end
         if any(Prob.CheckNaN == [2 3 4])
            [J2,rx]=FDJac2(x(1:N), Prob, 'nlp_r', rx, varargin{:} );%Using splines
            %if exist('csapi') & exist('csaps') & exist('spaps')
            %   [J2,rx]=FDJac2(x, Prob, 'nlp_r', rx,varargin{:});% Using splines
            %else
            %   fprintf('Can not find directory SPLINES\n')
            %   fprintf('Have you got a license for Spline Toolbox?\n')
            %   fprintf('Using finite differece routine FDJac.m\n')
            %   [J2,rx]=FDJac(x,Prob,'nlp_r',rx, varargin{:}); % FD Algorithm
            %end
         elseif Prob.CheckNaN == 5
            % Using complex number technique
            [J2,rx]=FDJac3(x(1:N), Prob, 'nlp_r', rx, varargin{:} ); 
         else
            [J2,rx]=FDJac(x(1:N), Prob, 'nlp_r', rx, varargin{:} ); % FD Algorithm
         end   
         % Merge analytic and numerical dc
         J(isnan(J)) = J2(isnan(J));
         LS_x = x;
         LS_r = rx;
      end
   end
end
if ~isempty(wLS)
   if any(size(wLS)==1) & any(size(wLS)==size(J,1))
      % Scaling matrix given as a vector, length number of residuals
      wLS=wLS(:);
      % J=(wLS*ones(1,size(J,2))).*J;
      if Prob.LargeScale
         % Avoid forming a scaling matrix
         for i = 1:size(J,1)
             J(i,:) = J(i,:)*wLS(i);
         end
      else
         J=diag(wLS)*J;
      end
   elseif size(wLS,2)==size(J,2) & size(wLS,1)==size(J,2) 
      % Must be square weighting matrix, size: number of variables
      J=J*wLS;
   elseif size(wLS,2)==size(J,2) & size(wLS,1)==size(J,1) 
      % Must be weighting matrix with same size as Jacobian
      J=J.*wLS;
   else  % Trouble!!!, ignore weights
      fprintf('nlp_J: Error!!! Given weight vector/matrix has wrong size!\n');
      size(wLS)
      size(J)
      wLS
      wLS=[];
   end
end

% Save computed Jacobian in global, along with the x.
LS_xJ=x;
LS_J=J;

% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation 
% 981018   hkh   Added weighted least squares handling
% 981023   hkh   Make the routine more general. Use globals.
%                Return J if already computed (used by g and H routines).
% 981024   hkh   Wrong logic, must also save x for J computation, as LS_xJ
% 981025   hkh   Make separable nonlinear least squares work
% 981027   hkh   Numerical differences possible
% 981028   hkh   Change to use general code FDJac2 for numerical differences
% 981029   hkh   Check that rx = r(x) is computed in the correct point x
% 981123   mbk   Check if spline routines is in path. Call to new routine
%                FDJac.m     
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
% 990523   hkh   Dangerous change for minimax. Use N=Prob.N to determine how
%                many variables to send to user routine (to allow for extra
%                variables in x vector).
% 990909   hkh   Avoid structural information handling
% 011203   hkh   Check for Prob.NumDiff > 5 (internal derivatives)
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020416   hkh   Do not set J empty if NumDiff == 6, will not be called by SOL
% 030127   hkh   Check for NaN elements, estimate numerically
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 031213   hkh   Revise weight handling, avoid forming diag(wLS) if LargeScale
% 040331   hkh   Use 2nd parameter rx from FDJac-routines, J could get []
% 040331   hkh   Set LS_x, LS_r using rx output from FDJac-routines
% 040407   hkh   Always call FDNum2, no check on spline TB routines
% 040526   hkh   Use x(1:N) in all function calls
% 040901   med   getvalue lower case