%
% TOMLAB gateway routine 
% Used as interface to filterSQP and MINLP
%
% nlp_d2L computes the Hessian 
%    
%       d2f(x) - lam' * d2c(x)
%
% to the Lagrangian function,
%
%   L(x,lam) =   f(x) - lam' * c(x)
%
% If input argument phase == 1, then only the second part lam'*d2c
% is calculated. If phase == 2, the full Hessian is calculated
%
% function d2L=nlp_d2L(x, lam, phase, Prob)
%
% Anders Göran, Tomlab Optimization Inc, E-mail: anders@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.4.0$
% Written Nov 11, 2002.   Last modified Oct 17, 2004.
%

function d2L=nlp_d2L(x, lam, phase, Prob)

global n_H n_d2c NARG

x=x(:);
d2H = [];

if(isempty(lam))
   % If no nonlinear c/s
   d2L = nlp_H(x,Prob);
else
   
   if phase==1 % Only second part of Hessian
      d2L = nlp_d2c(x,lam,Prob);
   elseif phase==2 % Hessian of full Lagrangian
      d2H = nlp_H(x,Prob);
      if isempty(d2H)
         d2L = - nlp_d2c(x,lam,Prob);
      else
         d2L = d2H - nlp_d2c(x,lam,Prob);
      end
   end
end

% filterSQP and MINLP needs transpose: d2L = d2L';

% if Prob.LargeScale
%    d2L = sparse(d2L');
%    %if ~issparse(d2L), d2L = sparse(d2L); end
% else
%    d2L = full(d2L');
%    %if issparse(d2L), d2L = full(d2L); end
% end

d2L = sparse(d2L);

if Prob.Warning
   if(any(find(d2L-d2L')))
    if Prob.ADCons == -1 | Prob.ADObj == -1
        d2L = 0.5*(d2L+d2L'); % ONLY IF MAD IS USED.
    else
      fprintf('\n\n')
      if isempty(lam)
         fprintf('The Hessian H of the objective is not symmetric!!!')
         fprintf('\n\n')
         fprintf('You must either correct your computations so that ')
         fprintf('the difference H-H'' is of the order of 1E-16 at most')
         fprintf('\n')
         fprintf('or make the Hessian symmetric by doing the following trick:')
         fprintf('\n')
         fprintf('H = 0.5*(H+H'');')
         fprintf('\n\n')
      elseif phase==1 | isempty(d2H)
         fprintf('The weighted Hessian d2L of the constraints ')
         fprintf('is not symmetric!!!')
         fprintf('\n\n')
         fprintf('You must either correct your computations so that ')
         fprintf('the difference d2L-d2L'' is of the order of 1E-16 at most')
         fprintf('\n')
         fprintf('or make d2L symmetric by doing the following trick:')
         fprintf('\n')
         fprintf('d2L = 0.5*(d2L+d2L'');')
         fprintf('\n\n')
      else
         if(any(find(d2H-d2H')))
            fprintf('The Hessian H of the objective is not symmetric!!!')
            fprintf('\n\n')
            fprintf('You must either correct your computations so that ')
            fprintf('the difference H-H'' is of the order of 1E-16 at most or')
            fprintf('\n')
            fprintf('make the Hessian symmetric by doing the following trick:')
            fprintf('\n')
            fprintf('H = 0.5*(H+H'');')
            fprintf('\n\n')
         end
         if(any(find((d2L-d2H)-(d2L'-d2H'))))
            fprintf('The weighted Hessian d2L of the constraints ')
            fprintf('is not symmetric!!!')
            fprintf('\n\n')
            fprintf('You must either correct your computations so that ')
            fprintf('the difference d2L-d2L'' is of the order of 1E-16 at most')
            fprintf('\n')
            fprintf('or make d2L symmetric by doing the following trick:')
            fprintf('\n')
            fprintf('d2L = 0.5*(d2L+d2L'');')
            fprintf('\n\n')
         end
      end
      error('Cannot proceed, Second order Lagrangian must be symmetric!!!');
      %error('nonsymmetric')
    end
   end
end

% if( all( ~ ( isinf(x) | isnan(x) ) ) & any(isinf(d2L)|isnan(d2L) ) )
% disp('nlp_d2L: x has no Inf/NaN but d2L does!')
% keyboard;
% end

% MODIFICATION LOG:
%
% 021111 ango Wrote file
% 021112 ango Added handling of sparse Hessians
% 021223 hkh  Faster to avoid tests before doing sparse/full/transpose
% 021229 hkh  Change comments
% 021231 ango Always return sparse Hessian and transpose removed
% 030127 hkh  Display warning if not symmetric
% 040412 hkh  Only check if symmetric if Prob.Warning true
% 041016 hkh  Write long text and stop computations if unsymmetric d2L
% 041017 hkh  Correcting symmetry of d2L if MAD is used
