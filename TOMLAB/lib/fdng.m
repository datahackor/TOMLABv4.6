%
% Finite-Difference Numerical approximation of Gradient. (FDNG)
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% function g = fdng(x, Prob, g, fx, varargin)
%
% If g is nonempty, estimate any elements of g set to NaN.
%
% The numerical step is given by setting Prob.GradTolg
% It is either a scalar or a vector of length(x)
% If isempty(Prob.GradTolg), then the scalar Prob.optParam.DiffInt is used
%
% If both parameters are empty, then tolerances are estimated. This
% automated algorithm takes both forward and backward steps, and therefore
% might step outside simple bounds on parameters, if the initial x is chosen
% to be on the boundary for some component.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.5.0$
% Written Oct 2, 1998.     Last modified Nov 27, 2004.
%

function g = fdng(x, Prob, g, fx, varargin)

n = Prob.N;
if isempty(g)
   ALLg = 1;
   g = zeros(n,1);
else
   ALLg = 0;
   ix = find(isnan(g));
   if isempty(ix), return, end
end

if nargin < 4, fx = []; end

x = x(:);

if isempty(fx)
   fx = nlp_f(x, Prob, varargin{:});
end

% gTol is a vector of relative intervals, See Pract. Opt. page 130

global gTol

if isempty(gTol)
   gTol = Prob.GradTolg;

   if isempty(gTol)
      h = Prob.optParam.DiffInt;
      if h > 0
         gTol = h*ones(n,1);
      end
   elseif gTol(1) <= 0
      gTol = [];
   else
      gTol = gTol(:);
   end
end
if length(gTol) >= n
   x_L  = Prob.x_L(:);
   if isempty(x_L)
      x_L = -inf*ones(n,1);
   elseif length(x_L) < n
      x_L(length(x_L)+1:n,1) = -inf;
   end
   x_U  = Prob.x_U(:);
   if isempty(x_U)
      x_U =  inf*ones(n,1);
   elseif length(x_U) < n
      x_U(length(x_U)+1:n,1) =  inf;
   end
   % If close to a bound, make the step in the opposite direction, changing sign
   ixx=find(x(1:n)+gTol(1:n) > x_U | x(1:n)+gTol(1:n) < x_L);
   gTol(ixx)=-gTol(ixx);
end

if isempty(gTol)

   % First time this routine is called in the current optimization.
   % Run the FD algoritmh to determine the relative intervals.
   % rho = 10;
   eps_R = eps;
   eps_A = eps_R*(1 + abs(fx));
   eps_phi = 1E-14; % zero toleranse for phi_F, phi_B and PHI. 
   M = 1E10; % Large number.
   K = 6;
   PriLev = Prob.PriLevOpt;
   gTol =zeros(n,1);
   if ALLg
      ix = [1:n]; % Must generate indices for all elements
   end
   
   for i = 1:length(ix)
      dim        = ix(i);
      Prob.FDVar = dim;
      z          = x(dim);
      % Flags:
      FD3 = 1; FD4 = 1; FD5 = 1; FD6 = 1;
   
      % FD1: Initialization
      % h_line  = 2*(eta + abs(z))*sqrt(eps_A/(omega+abs(fx)));
      % eta = 1, omega = 1
      h_line  = 2*(1 + abs(z))*sqrt(eps_A/(1+abs(fx)));
      h(1)    = 10*h_line; % h(i+1) corresponds to h(i) in Pract. Opt.
      k       = 0;
      x(dim)  = z + h(1);
      fx_ph   = nlp_f(x, Prob, varargin{:}); % f(x+h)
      x(dim)  = z - h(1);
      fx_mh   = nlp_f(x, Prob, varargin{:}); % f(x-h)
      x(dim)  = z;
      phi_F   = (fx_ph-fx) / h(1);
      phi_B   = (fx-fx_mh) / h(1);
      phi_C   = (fx_ph-fx_mh) / (2*h(1));
      PHI     = (fx_ph - 2*fx + fx_mh) / ( h(1)^2 );
      if abs(phi_F) <= eps_phi
         C_phi_F = M;
      else
         C_phi_F = 2*eps_A / ( h(1) * abs(phi_F) );
      end   
      if abs(phi_B) <= eps_phi
         C_phi_B = M;
      else
         C_phi_B = 2*eps_A / ( h(1) * abs(phi_B) );
      end   
      if abs(PHI) <= eps_phi
         C_PHI = M;
      else
         C_PHI   = 4*eps_A / ( h(1)^2 * abs(PHI) );
      end   
      h_s = -1;
      
      
      % FD2: Decide whether to accept the initial interval
      if max(C_phi_F,C_phi_B) <= 0.1
         h_s = h(1);
      end
      if ( C_PHI >= 0.001 ) & ( C_PHI <= 0.1 )
         h_phi = h(1);
         FD3   = 0; FD4 = 0;  % GOTO STEP FD5
      elseif C_PHI < 0.001
         FD3   = 0;           % GOTO STEP FD4
      end
      
      % FD3: Increase h
      while FD3 & (k < K)
         k = k + 1;
         h(k+1) = 10*h(k); 
         % Compute the associated finite-difference estimates and their
         % relative errors:
         x(dim)  = z + h(k+1);
         fx_ph   = nlp_f(x, Prob, varargin{:}); % f(x+h)
         x(dim)  = z - h(k+1);
         fx_mh   = nlp_f(x, Prob, varargin{:}); % f(x-h)
         x(dim)  = z;
         phi_F   = (fx_ph-fx) / h(k+1);
         phi_B   = (fx-fx_mh) / h(k+1);
         PHI     = (fx_ph - 2*fx + fx_mh) /  h(k+1)^2 ;
         if abs(phi_F) <= eps_phi
            C_phi_F = M;
         else
            C_phi_F = 2*eps_A / ( h(k+1) * abs(phi_F) );
         end
         if abs(phi_B) <= eps_phi
            C_phi_B = M;
         else
            C_phi_B = 2*eps_A / ( h(k+1) * abs(phi_B) );
         end   
         if abs(PHI) <= eps_phi
            C_PHI = M;
         else
            C_PHI   = 4*eps_A / ( h(k+1)^2 * abs(PHI) );
         end
         if (h_s < 0) & (max(C_phi_F,C_phi_B) <= 0.1)
            h_s = h(k+1);
         end
         if C_PHI <= 0.1
            h_phi = h(k+1);
            FD3 = 0; FD4 = 0; % GOTO STEP FD5
         end
      end
      if k == K
         FD4 = 0; FD5 = 0; % GOTO STEP FD6
      end
      
      while FD4 & (k < K) 
         % FD4: Decrease h
         k       = k + 1;
         h(k+1)  = h(k)/10;
         % Compute the associated finite-difference estimates and their
         % relative errors:
         x(dim)  = z + h(k+1);
         fx_ph   = nlp_f(x, Prob, varargin{:}); % f(x+h)
         x(dim)  = z - h(k+1);
         fx_mh   = nlp_f(x, Prob, varargin{:}); % f(x-h)
         x(dim)  = z;
         phi_F   = (fx_ph-fx) / h(k+1);
         phi_B   = (fx-fx_mh) / h(k+1);
         PHI     = (fx_ph - 2*fx + fx_mh) / ( h(k+1)^2);
         if abs(phi_F) <= eps_phi
            C_phi_F = M;
         else
            C_phi_F = 2*eps_A / ( h(k+1) * abs(phi_F) );
         end
         if abs(phi_B) <= eps_phi
            C_phi_B = M;
         else
            C_phi_B = 2*eps_A / ( h(k+1) * abs(phi_B) );
         end
         if abs(PHI) <= eps_phi
            C_PHI = M;
         else
            C_PHI = 4*eps_A / ( h(k+1)^2 * abs(PHI) );
         end
         if C_PHI > 0.1
            h_phi = h(k);
            FD4   = 0; % GOTO STEP FD5
         else
            if max(C_phi_F,C_phi_B) < 0.1
               h_s = h(k+1);
            end
            if ( C_PHI >= 0.001 ) & ( C_PHI <= 0.1 )
               h_phi = h(k+1);
               FD4 = 0; % GOTO STEP FD5
            end
         end
      end
      if k == K
         FD5 = 0; % GOTO STEP FD6
      end
      
      if FD5
         % FD5: Compute the estimate of the optimal interval
         % HKH 020404,020616 - Avoid division with 0 

         if PHI == 0 | isinf(PHI)
            h_F = 200 * h(k+1);
         else
            h_F = 2 * sqrt( eps_A / abs(PHI) );
         end
         x(dim)  = z + h_F;
         fx_ph   = nlp_f(x, Prob, varargin{:}); % f(x+h)
         phi     = (fx_ph-fx) / h_F;
         x(dim)  = z;
         
         % Compute the estimated error bound
         E_F = h_F*abs(PHI)/2 + 2*eps_A/h_F;
         
         % Compute the difference between phi and phi_C(h_phi)
         x(dim)  = z + h_phi;
         fx_ph   = nlp_f(x, Prob, varargin{:}); % f(x+h)
         x(dim)  = z - h_phi;
         fx_mh   = nlp_f(x, Prob, varargin{:}); % f(x-h)
         phi_C   = (fx_ph-fx_mh) / (2*h_phi);
         E_line  = abs(phi - phi_C); 
         x(dim)  = z;
         
         if max(E_F,E_line) <= 0.5*abs(phi)
            FD6 = 0; % the algorithm terminates successfully
         else
            FD6 = 0;
            if PriLev > 0
               disp('Termination with an error condition in FD5 in fdng.m')
            end
         end
      end
   
      if FD6
         % FD6: Check unsatisfactory cases
         if PriLev > 2
            disp(' FD6: Check cases when no satisfactory interval was found');
         end   
         if h_s < 0
            if PriLev > 2
               disp(' f appear to be nearly constant')
            end
            h_F   = h_line;
            phi   = 0;
            PHI   = 0;
            E_F   = 0;
         elseif ( C_PHI > 0.1 ) & ( h_s > 0 )
            if PriLev > 2
               disp(' f appears to be odd or nearly linear')
            end
            h_F    = h_s;
            x(dim) = z + h_F;
            fx_ph  = nlp_f(x, Prob, varargin{:}); % f(x+h)
            phi    = (fx_ph-fx) / h_F;
            PHI    = 0;
            E_F    = 2*eps_A/h_F;
            x(dim) = z;
         else
            if PriLev > 2
               disp(' d2f/dx2 appears to be increasing rapidly as h decreases')
            end
            h_F    = h(K+1);
            x(dim) = z + h_F;
            fx_ph  = nlp_f(x, Prob, varargin{:}); % f(x+h)
            phi    = (fx_ph-fx) / h_F;
            x(dim) = z - h_F;
            fx_mh  = nlp_f(x, Prob, varargin{:}); % f(x-h)
            PHI    = (fx_ph - 2*fx + fx_mh) / ( h_F^2 );
            E_F    = h_F*abs(PHI)/2 + 2*eps_A/h_F;
            x(dim) = z;
         end
      end
      
      g(dim)       = phi;
      gTol (dim)   = h_F/(1+abs(z));
   end % for
   
elseif ALLg
   % Do not run the FD algorithm, compute the gradient by using the
   % relative intervals computed for x_0.
   % All elements are estimated
   if length(gTol) < n
      gTol = gTol(1);
      for dim = 1:n
         Prob.FDVar = dim;
         z          = x(dim);
         h          = gTol * (1+abs(z));
         x(dim)     = z + h;
         fx_ph      = nlp_f(x, Prob, varargin{:}); % f(x+h)
         g(dim)     = (fx_ph-fx) / h;
         x(dim)     = z;
      end   
   else
      for dim = 1:n
         Prob.FDVar = dim;
         z          = x(dim);
         h          = gTol(dim) * (1+abs(z));
         x(dim)     = z + h;
         fx_ph      = nlp_f(x, Prob, varargin{:}); % f(x+h)
         g(dim)     = (fx_ph-fx) / h;
         x(dim)     = z;
      end
   end   
else
   % Do not run the FD algorithm, compute the gradient by using the
   % relative intervals computed for x_0.
   % Only NaN elements are estimated

   if length(gTol) < n
      gTol = gTol(1);
      for i = 1:length(ix)
         dim        = ix(i);
         Prob.FDVar = dim;
         z          = x(dim);
         h          = gTol * (1+abs(z));
         x(dim)     = z + h;
         fx_ph      = nlp_f(x, Prob, varargin{:}); % f(x+h)
         g(dim)     = (fx_ph-fx) / h;
         x(dim)     = z;
      end   
   else
      for i = 1:length(ix)
         dim        = ix(i);
         Prob.FDVar = dim;
         z          = x(dim);
         h          = gTol(dim) * (1+abs(z));
         x(dim)     = z + h;
         fx_ph      = nlp_f(x, Prob, varargin{:}); % f(x+h)
         g(dim)     = (fx_ph-fx) / h;
         x(dim)     = z;
      end
   end   
end   

% MODIFICATION LOG
%
% 981026  hkh  Changed name to fdng. Changed global name to gTol   
%              Use isempty(gTol ) as flag if first attempt
% 981027  hkh  Return empty g in case of bound error on x
% 981104  mbk  Rewritten by following the algorithm description in P.O. 1993
% 981104  hkh  Changed PriLev > 0 to PriLev > 2
% 981123  mbk  Missing semicolon added.
% 990306  hkh  Safeguard against x slightly outside bounds
% 990626  hkh  Avoid feval
% 000910  hkh  Written more efficiently, avoid arrays and zero multiplications
% 000911  hkh  Using new step size choice
% 020404  hkh  Fix labled HKH to avoid division with 0 for weird EGO problem
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 020616  hkh  Another fix for label HKH, avoid 0-division in EGO, PHI is NaN
% 030127  hkh  Added g as input, estimate NaN elements
% 030131  hkh  Initialize g as column vector when empty
% 030523  hkh  x(dim) = z not set if gTol length < n
% 040411  hkh  Send Prob.FDVar with the variable index changed
% 040413  hkh  Check length of x_L and x_U and adjust gTol
% 040413  hkh  gTol check only if Prob.GradTolg has length n
% 040505  hkh  ix used for two different purposes
% 040610  med  Print out bug removed
