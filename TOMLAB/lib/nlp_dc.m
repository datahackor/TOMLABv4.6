%
% TOMLAB gateway routine for computation of the constraint Jacobian matrix
%
% nlp_dc calls the routine Prob.USER.dc either as 
%           dc=feval(Prob.USER.dc, x) or
%           dc=feval(Prob.USER.dc, x, Prob) or
%           dc=feval(Prob.USER.dc, x, Prob, varargin{:})
% depending on the number of inputs in Prob.USER.dc
%
% function dc = nlp_dc(x, Prob, varargin)
%
% dc is a mL x n matrix, where mL = number of nonlinear constraints
%                               n = number of variables
%
% The global counter variable n_dc is incremented
% If Prob.ConsDiff > 0, dc is estimated numerically
% If Prob.CheckNaN > 0, NaN elements in dc are estimated numerically
%
% Numerical method is dependent on Prob.ConsDiff and Prob.CheckNaN
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.1$
% Written Oct 10, 1998.   Last modified May 26, 2004.
%

function dc = nlp_dc(x, Prob, varargin)

%nargin;

global n_dc NARG

global NLP_xc NLP_xdc NLP_c NLP_dc  % Communication nlp_c/dc

if Prob.simType > 0
   [g,dc] = sim_gdc(x, Prob, varargin{:} );
   return
end

x=x(:);

N = min(length(x),Prob.N);

if ~isempty(NLP_xdc)
   if length(x)~=length(NLP_xdc)
      NLP_xdc=[];
   elseif all(x==NLP_xdc)
      dc=NLP_dc;
      return
   end
end

n_dc=n_dc+1;

Func = Prob.USER.dc;

if Prob.ADCons == 1
   global mad_c
   if ~isempty(NLP_xc)
      if all(x==NLP_xc)
         cx=NLP_c;
      else
         cx=[];
      end
   else
      cx=[];
   end
   if isempty(cx)
      if isempty(NARG)
         p = xnargin(Func);
      else
         p = NARG(4);
      end
      %cx = nlp_c(x, Prob, varargin);
      %dc=getinternalderivs(mad_c);
      Func  = Prob.USER.c;
      fdvar = Prob.FDVar;
      if fdvar == 0 | length(fdvar) == N
         Z = speye(N);
      else
         z = zeros(N,1);
         z(fdvar) = 1;
         Z = spdiags(z,0,N,N);
      end
      if p > 2
         %madc=feval(Func, fmad(x(1:N),speye(length(x(1:N)))),Prob,varargin{:});
         %madc=feval(Func, fmad(x(1:N),speye(N)), Prob,varargin{:});
         madc=feval(Func, fmad(x(1:N),Z), Prob,varargin{:});
      elseif p > 1
         madc=feval(Func, fmad(x(1:N),Z), Prob);
      else
         madc=feval(Func, fmad(x(1:N),Z));
      end
      %c=getvalue(mad_c);
      dc=getinternalderivs(madc);
   else
      dc=getinternalderivs(mad_c);
   end
   %NLP_xc = x;
   %NLP_c  = cx;
elseif Prob.ADCons == -1  
   global mad_dc
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(5);
   end
   if p > 2
      mad_dc=feval(Func, fmad(x(1:N),speye(N)), Prob, varargin{:});
   elseif p > 1
      mad_dc=feval(Func, fmad(x(1:N),speye(N)), Prob);
   else
      mad_dc=feval(Func, fmad(x(1:N),speye(N)));
   end
   dc=getvalue(mad_dc);

elseif Prob.ADCons == 2
   global cJPI ad_c ad_dc ad_cM
   dc=full(ad_dc);
   % elseif Prob.ADCons == -2 Maybe not possible to handle with ADMAT
elseif isempty(Func) | Prob.ConsDiff > 0
   if isempty(Prob.USER.c)
      dc=zeros(0,N);
      NLP_xdc = [];
      NLP_dc  = [];
      return
   else
      % Numerical difference computation of dc
      if ~isempty(NLP_xc)
         if all(x==NLP_xc)
            cx=NLP_c;
         else
            cx=[];
         end
      else
         cx=[];
      end
      if any(Prob.ConsDiff == [2 3 4])
         [dc,cx]=FDJac2(x(1:N), Prob, 'nlp_c', cx, varargin{:}); % Using splines
         %if exist('csapi') & exist('csaps') & exist('spaps')
         %   [dc,cx]=FDJac2(x, Prob, 'nlp_c', cx, varargin{:}); % Using splines
         %else
         %   fprintf('Can not find directory SPLINES\n')
         %   fprintf('Have you got a license for Spline Toolbox?\n')
         %   fprintf('Using finite differece routine FDJac.m\n')
         %   [dc,cx]=FDJac(x, Prob, 'nlp_c', cx, varargin{:}); % FD Algorithm
         %end
      elseif Prob.ConsDiff == 5
         % Using complex number technique
         [dc,cx]=FDJac3(x(1:N), Prob, 'nlp_c', cx, varargin{:}); 
      else
         [dc,cx]=FDJac(x(1:N), Prob, 'nlp_c', cx, varargin{:}); % FD Algorithm
      end   
      NLP_xc = x;
      NLP_c  = cx;
   end
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(5);
   end

   if p>2
      dc = feval(Func, x(1:N), Prob, varargin{:});
   elseif p==2
      dc = feval(Func, x(1:N), Prob);
   else
      dc = feval(Func, x(1:N));
   end
   if Prob.CheckNaN > 0
      [iN,jN,dcN] = find(isnan(dc));
      if ~isempty(iN)    % There are elements set to NaN, to be estimated
         Prob.ConsPattern = sparse(iN,jN,dcN,size(dc,1),N);
         if ~isempty(NLP_xc)
            if all(x==NLP_xc)
               cx=NLP_c;
            else
               cx=[];
            end
         else
            cx=[];
         end
         if any(Prob.CheckNaN == [2 3 4])
            [dcN,cx]=FDJac2(x(1:N), Prob, 'nlp_c', cx, varargin{:}); 
            %if exist('csapi') & exist('csaps') & exist('spaps')
            %   % Using splines
            %   [dcN,cx]=FDJac2(x, Prob, 'nlp_c', cx, varargin{:}); 
            %else
            %   fprintf('Can not find directory SPLINES\n')
            %   fprintf('Have you got a license for Spline Toolbox?\n')
            %   fprintf('Using finite differece routine FDJac.m\n')
            %   [dcN,cx]=FDJac(x, Prob, 'nlp_c',cx,varargin{:}); % FD Algorithm
            %end
         elseif Prob.CheckNaN == 5
            % Using complex number technique
            [dcN,cx]=FDJac3(x(1:N), Prob, 'nlp_c', cx, varargin{:}); 
         else
            [dcN,cx]=FDJac(x(1:N), Prob, 'nlp_c', cx, varargin{:});%FD Algorithm
         end   
         % Merge analytic and numerical dc
         dc(isnan(dc)) = dcN(isnan(dc));
         NLP_xc = x;
         NLP_c  = cx;
      end
   end
end
%numnan=sum(isnan(dc(:)));
%numnan
%if numnan > 0
%   fprintf('NaN Elements in dc %d\n',numnan);
%end
%if any(any(isnan(dc)))
%        keyboard
%end

NLP_xdc=x;
NLP_dc=dc;
% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation 
% 981028   hkh   Add code for numerical differences using splines
% 981102   hkh   Just using NLP_c and NLP_xc. Change empty setting of dc
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
% 001114   hkh   Use Prob.ConsDiff instead of NumDiff for explicit control
% 001212   hkh   Use Prob.ConsDiff > 5 for solver estimate of gradient
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020413   hkh   Just send the first Prob.N x variables to Func
% 020416   hkh   Do not set dc empty if ConsDiff == 6, will not be called by SOL
% 020429   hkh   length(x) < Prob.N if nnJac < nnObj for SOL solvers, use min
% 030127   hkh   Check for NaN elements, estimate numerically
% 030524   hkh   Add handling of simulation problems
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 040331   hkh   Use 2nd parameter cx from FDJac-routines, dc could get []
% 040331   hkh   Save NLP_dc, and use if x==NLP_xdc
% 040407   hkh   Always call FDJac2, no check on spline TB routines
% 040409   hkh   More efficent handling of repeated nlp_d2c calls for AD-num.der
% 040526   hkh   Use x(1:N) in all function calls
% 040901   med   getvalue lower case