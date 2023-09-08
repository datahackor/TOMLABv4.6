%		optim_fgH.m
%
% function f = optim_fgH(x, Prob)
%
% optim_fgH is used to implement the OPT TB 2.0 interface
%
% The function f is returned 
%
% If the gradient is computed, it is stored in the global variable
% NLP_g (and the corresponding x value in NLP_xg) 
%
% If the Hessian is computed, it is stored in the global variable
% NLP_H (and the corresponding x value in NLP_xH) 
%
% optim_fgH is called from the TOMLAB gateway function nlp_f.
%
% Kenneth Holmstrom, Tomlab Optimization AB, E-mail: hkh@tomlab.biz.
% Copyright (c) 1999-2004 by Tomlab Optimization AB, Sweden. $Release: 4.6.0$
% Written July 29, 1999.     Last modified Dec 22, 2004.
%

function f = optim_fgH(xNew, Prob)

global NLP_xg NLP_g NLP_xH NLP_H

Func    = Prob.OPTTB.f;
aF      = Prob.OPTTB.funArgIn;
oF      = Prob.OPTTB.funArgOut;
NumDiff = Prob.NumDiff;
x       = Prob.OPTTB.x;
x(:)    = xNew;                   % To handle cases when x is a matrix

if Prob.OPTTB.M7 == 0

   if NumDiff >0 | Prob.ADObj > 0 
      % Gradient should not be computed
      if oF < 0
         f=eval(Func);
      else
         if aF > 1
            f=feval(Func, x, Prob.varargin{:} );
         else
            f=feval(Func, x);
         end
      end
      NLP_xg=[];
      NLP_g=[];
      NLP_xH=[];
      NLP_H=[];
   elseif NumDiff | Prob.ADObj < 0       
      % Hessian should not be computed
      if oF < 0
         f=eval(Func);
         NLP_g = eval(Prob.OPTTB.g);
      else
         if aF > 1
            [f,NLP_g] = feval(Func, x, Prob.varargin{:});
         else
            [f,NLP_g] = feval(Func, x);
         end
      end
      NLP_xg=x(:);
   else
      if oF < 0
         f=eval(Func);
         NLP_g = eval(Prob.OPTTB.g);
         if oF == -2
            NLP_H = eval(Prob.OPTTB.H);
        else
            NLP_H = [];
        end
      else
         if oF > 2
            if aF > 1
               [f,NLP_g,NLP_H] = feval(Func, x, Prob.varargin{:});
            else
               [f,NLP_g,NLP_H] = feval(Func, x);
            end
         else
            if aF > 1
               [f,NLP_g] = feval(Func, x, Prob.varargin{:});
            else
               [f,NLP_g] = feval(Func, x);
            end
            NLP_H = [];
         end
         HFunc = Prob.HessMult;
         if ~isempty(HFunc)
            n = Prob.N;
            y = zeros(n,1);
            Hinfo = NLP_H;
            NLP_H = sparse(n,n);
            % Sick way to obtain the full Hessian - works, but slow 
            for i=1:n
                y(i) = 1;
                NLP_H(:,i) = feval(HFunc, Hinfo, y, Prob.varargin{:});
                y(i) = 0;
            end
         end
      end
      NLP_xg=x(:);
      NLP_xH=x(:);
   end

else
   if NumDiff >0 | Prob.ADObj > 0 
      % Gradient should not be computed
      if aF > 1
         f=Func( x, Prob.varargin{:} );
      else
         f=Func( x);
      end
      NLP_xg=[];
      NLP_g=[];
      NLP_xH=[];
      NLP_H=[];
   elseif NumDiff | Prob.ADObj < 0       
      % Hessian should not be computed
      if aF > 1
         [f,NLP_g]=Func( x, Prob.varargin{:} );
      else
         [f,NLP_g]=Func( x);
      end
      NLP_xg=x(:);
   else
      if oF > 2
         if aF > 1
            [f,NLP_g,NLP_H] = Func( x, Prob.varargin{:} );
         else
            [f,NLP_g,NLP_H] = Func( x);
         end
      else
         if aF > 1
            [f,NLP_g]=Func( x, Prob.varargin{:} );
         else
            [f,NLP_g]=Func( x);
         end
         NLP_H = [];
      end
      HFunc = Prob.HessMult;
      if ~isempty(HFunc)
         n = Prob.N;
         y = zeros(n,1);
         Hinfo = NLP_H;
         NLP_H = sparse(n,n);
         % Sick way to obtain the full Hessian - works, but slow 
         for i=1:n
             y(i) = 1;
             NLP_H(:,i) = feval(HFunc, Hinfo, y, Prob.varargin{:});
             y(i) = 0;
         end
      end
      NLP_xg=x(:);
      NLP_xH=x(:);
   end
end

% MODIFICATION LOG:
%
% 030114 hkh Major revision, avoid check of nargin(Func) here
% 030115 hkh Add handling of cell, inline and function_handle input
% 030126 hkh Handle the case when x is a matrix, and not a vector
% 030128 hkh Handle the HessMult option, generate Hessian
% 031101 hkh Change AutoDiff to new field ADObj, add test for Hessian for MAD
% 041119 hkh Test oF, allow both two or three outputs from user function
% 041222 hkh Handle new type of function handle in Matlab 7.x

