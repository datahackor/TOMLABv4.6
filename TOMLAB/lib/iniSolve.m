% init routine. Called before solving an optimization problem
%
% Initializes global variables
%
% function Prob = iniSolve(Prob,solvType,ObjDers,ConsDers);
%
% solvType   Solver type
% ObjDers    Level of derivatives needed for the objective, 0,1,2
% ConsDers   Level of derivatives needed for the constraints, 0,1,2
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.3.1$
% Written June 1, 1997.  Last modified June 3, 2004.

function Prob = iniSolve(Prob,solvType,ObjDers,ConsDers);

if nargin < 4
   ObjDers = 2;
   if nargin < 3
      ConsDers = 2;
      if nargin < 2
         if isfield(Prob,'solvType')
            solvType=Prob.solvType;
         else
            solvType=7;
         end
      end
   end
end
if ~isfield(Prob,'optParam') | isempty(Prob.optParam)
   fprintf('\n\n')
   fprintf('Field optParam is not correctly defined in the Prob structure\n');
   fprintf('Use tomRun, or first call:\n');
   fprintf('\n');
   fprintf('Prob.optParam = optParamDef(Solver,Prob.probType);\n');
   fprintf('\n');
   fprintf('where Solver is the name of the solver, e.g. ''snopt'' and\n');
   fprintf('Prob.probType is the Tomlab problem type if the Prob structure\n');
   fprintf('is correctly defined.\n');
   fprintf('It is also possible to check the full Prob structure and call: \n');
   fprintf('\n');
   fprintf('Prob = ProbCheck(Prob,Solver);\n');
   fprintf('\n')
   error('iniSolve: Field optParam is not correctly defined')
end

% Call to set up catching of search directions and line search steps

% ---------
% inibuild
% ---------
% Initialization to build search directions, step lengths and 
% variable limits for solver routines. 
% Initialization of global variables sent to pbuild.m.
%
% Set flag BUILDP as true.
%
%

global p_dx alphaV X_min X_max F_X BUILDP PartSep
p_dx=[]; alphaV=[]; X_min=[]; X_max=[];
BUILDP=double(0);

PartSep=0;
if isfield(Prob.PartSep,'pSepFunc') & isfield(Prob.PartSep,'index')
   if Prob.PartSep.pSepFunc > 1 
      PartSep = 1;
   end
end
% PRESOLVE
if Prob.optParam.PreSolve
   Prob = preSolve(Prob);
end
smallA = DefPar(Prob,'smallA',0);
if smallA
   % Find small elements in A and remove them
   v         = max(max(abs(Prob.A)));
   [i,j,s]   = find(sparse(Prob.A));
   ix        = abs(s) < eps*v;
   A0        = sum(ix);
   if A0 > 0
      [m,n]  = size(Prob.A);
      if Prob.Warning | Prob.optParam.IterPrint | Prob.PriLevOpt > 0 
         fprintf('iniSolve: Delete %d small elements in A matrix\n',A0);
      end
      ix1    = find(ix == 0);
      % Create new A matrix
      Prob.A = sparse(i(ix1),j(ix1),s(ix1),m,n);
   end
end

% End of code from inibuild script

% Avoid building search directions if not using the GUI or MENU system
if Prob.MENU > 0
   if any(solvType==[2 5 7 8 9 10]) | Prob.LargeScale
   % Avoid collection the search steps for MIP,GLB,GLC and for large probs
      BUILDP=double(0);
   else
      switch Prob.Solver.Name
      case {'ucSolve','conSolve','MINOS','SNOPT','NPSOL','fmincon','constr'...
           ,'fminunc','clsSolve','leastsq','lsqnonlin','fminu'}
         BUILDP=double(1);
      otherwise
         BUILDP=2;
      end
   end
end
% Initialization for CPU time and real time elapsed
%global TIME0 TIME1

Prob.TIME0 = cputime;
Prob.TIME1 = clock;
Prob.FDVar = double(0);  % Init 0 value for indices in numerical derivatives

global n_f n_g n_H    % Count of function evals, gradient evals, Hessian evals
global n_c n_dc n_d2c % Count of constraint evals, constraint gradient evals
global n_J n_r n_d2r  % Count of Jacobian evals, residual evals, 2nd der evals

% Init counters to 0
n_f   = double(0); n_g=double(0); n_H=double(0); n_c=double(0); n_dc=double(0); 
n_d2c = double(0); n_r=double(0); n_J=double(0); n_d2r=double(0);

global LS_x LS_r LS_xJ LS_J LS_A wLS
LS_x=[]; LS_r=[]; LS_xJ=[]; LS_J=[]; LS_A=[]; wLS=[]; 

% User communication f,g,H,c,dc
global US_A US_B US_C US_D
US_A=[]; US_B=[]; US_C=[]; US_D=[];

% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f NLP_xc NLP_c NLP_xdc NLP_dc NLP_pSepIndex
NLP_x=[]; NLP_f=[]; NLP_xc=[]; NLP_c=[]; NLP_xdc = []; NLP_dc = [];

global NLP_g NLP_xg
NLP_g = []; NLP_xg = [];

global SepIndex SEP_z SEP_Jz
NLP_pSepIndex=[]; NLP_z=[]; NLP_Jz=[];

global gTol JTol 
gTol=[]; JTol=[]; 

global cJPI HPI JPI
cJPI=[]; HPI=[]; JPI=[];

global NARG NARGO
NARG=zeros(11,1);
NARGO=zeros(1,1);
z= Prob.USER.f;
if ~isempty(z), NARG(1)  = xnargin(z); end
if ~isempty(z), NARGO(1) = xnargout(z); end
z= Prob.USER.g;
if ~isempty(z), NARG(2)  = xnargin(z); end
z= Prob.USER.H;
if ~isempty(z), NARG(3)  = xnargin(z); end
z= Prob.USER.c;
if ~isempty(z), NARG(4)  = xnargin(z); end
z= Prob.USER.dc;
if ~isempty(z), NARG(5)  = xnargin(z); end
z= Prob.USER.d2c;
if ~isempty(z), NARG(6)  = xnargin(z); end
z= Prob.USER.r;
if ~isempty(z), NARG(7)  = xnargin(z); end
z= Prob.USER.J;
if ~isempty(z), NARG(8)  = xnargin(z); end
z= Prob.USER.fc;
if ~isempty(z), NARG(10) = xnargin(z); end
z= Prob.USER.gdc;
if ~isempty(z), NARG(11) = xnargin(z); end

% Check Automatic differentiation

% Might have to check if Prob.NumDiff or ConsDiff == 6, possible conflict
if Prob.ADObj ~= 0 | Prob.ADCons ~=0
   switch ObjDers
   case 0
     Prob.ADObj = 0; 
   case 1
     if Prob.ADObj < 0, Prob.ADObj = 0; end 
   case 2
     % Must use Numerical differentiation for 2nd order if MAD
     if Prob.ADObj == 1 & Prob.NumDiff == 0, Prob.NumDiff = -1; end 
   otherwise
     Prob.ADObj = 0; 
   end
   switch ConsDers
   case 0
     Prob.ADCons = 0; 
   case 1
     if Prob.ADCons < 0, Prob.ADCons = 0; end 
   case 2
     % Must use Numerical differentiation for 2nd order if MAD
     if Prob.ADCons == 1 & Prob.ConsDiff == 0, Prob.ConsDiff = -1; end 
   otherwise
     Prob.ADCons = 0; 
   end
   if abs(Prob.ADObj == 1) | abs(Prob.ADCons) == 1
      madOK = checkMAD(1);
      if madOK
         global mad_x mad_f mad_g mad_c mad_dc mad_r mad_J
         mad_x  = []; mad_f  = []; mad_g  = []; mad_c  = []; mad_dc = [];
         mad_r  = []; mad_J  = [];
      else
         Prob.ADObj  = double(0);
         Prob.ADCons = double(0);
         error('Fatal error in iniSolve, MAD not installed, cannot proceed!!!')
      end
   elseif abs(Prob.ADObj == 2) | abs(Prob.ADCons) == 2
      admatOK = checkADMAT(1);
      if admatOK
         global HPI JPI cJPI ad_f ad_g ad_H ad_c ad_dc ad_cM
         ad_f = []; ad_g = []; ad_H = []; ad_c = []; ad_dc = []; ad_cM = [];
         cJPI = [];
         % Check if HPI, JPI is noempty???
      else
         Prob.ADObj  = double(0);
         Prob.ADCons = double(0);
         error('Fatal error in iniSolve, ADMAT not installed, cannot proceed!')
      end
   end
end
if ConsDers > 0 & Prob.LargeScale
   if isempty(Prob.ConsPattern) & Prob.mNonLin > 0
      if Prob.Warning
         fprintf('===============\n'); 
         fprintf('TOMLAB WARNING: Large Scale problem, ');
         fprintf('but Prob.ConsPattern not defined'); 
         fprintf('\n'); 
         fprintf('===============\n'); 
         fprintf('Calling estConsPattern to estimate Prob.ConsPattern');  
         fprintf('\n'); 
         fprintf('\n'); 
      end
      Prob.ConsPattern = estConsPattern(Prob);
   end
end
if ObjDers > 1 & Prob.LargeScale
   if isempty(Prob.HessPattern)
      if Prob.Warning
         fprintf('===============\n'); 
         fprintf('TOMLAB WARNING: Large Scale problem, ');
         fprintf('but Prob.HessPattern not defined'); 
         fprintf('\n'); 
         fprintf('===============\n'); 
         fprintf('Calling estHessPattern to estimate Prob.HessPattern');  
         fprintf('\n'); 
         fprintf('\n'); 
      end
      Prob.HessPattern = estHessPattern(Prob);
   end
end
if ObjDers > 0 & Prob.LargeScale & NARG(7) > 0
   if isempty(Prob.JacPattern)
      if Prob.Warning
         fprintf('===============\n'); 
         fprintf('TOMLAB WARNING: Large Scale NLLS problem, ');
         fprintf('but Prob.JacPattern not defined'); 
         fprintf('\n'); 
         fprintf('===============\n'); 
         fprintf('Calling estJacPattern to estimate Prob.JacPattern');  
         fprintf('\n'); 
         fprintf('\n'); 
      end
      Prob.JacPattern = estJacPattern(Prob);
   end
end
if Prob.ConsDiff > 10 
   if ConsDers > 0
      if isempty(Prob.ConIx)
         Prob.ConIx    = findpatt(Prob.ConsPattern);
      end
   end
   Prob.ConsDiff = min(5,Prob.ConsDiff-10);
elseif Prob.ConsDiff < -10 
   %Prob.ConsDiff = min(-1,Prob.ConsDiff+10);
   Prob.ConsDiff = -1;
   if ConsDers > 1
      if isempty(Prob.ConIx)
         Prob.ConIx    = findpatt(Prob.ConsPattern);
      end
   end
end
if ObjDers > 1
   if ~isempty(Prob.d2LPattern) & isempty(Prob.HessPattern)
      Prob.HessPattern = Prob.d2LPattern;
   elseif isempty(Prob.d2LPattern) & ~isempty(Prob.HessPattern) ...
                                   & ~isempty(Prob.ConsPattern)
      Prob.d2LPattern  = Prob.HessPattern | maked2cPatt(Prob.ConsPattern);
      %d2L= full(Prob.d2LPattern)
   elseif isempty(Prob.d2LPattern) & ~isempty(Prob.HessPattern) ...
                                   & Prob.mNonLin == 0
      Prob.d2LPattern  = Prob.HessPattern;
   end
end
if Prob.NumDiff > 10 
   if ObjDers > 0
      if isempty(Prob.JacIx) & ~isempty(Prob.JacPattern)
         Prob.JacIx    = findpatt(Prob.JacPattern);
      end
   end
   if ObjDers > 1
      if isempty(Prob.HessIx) & ~isempty(Prob.HessPattern)
         [ix,iy,iv]       = find(Prob.HessPattern);
         iz               = find(ix <= iy);
         Prob.HessPattern = sparse(ix(iz),iy(iz),iv(iz),Prob.N,Prob.N);
         Prob.HessIx      = findpatt(Prob.HessPattern);
      end
   end
   Prob.NumDiff = min(5,Prob.NumDiff-10);
elseif Prob.NumDiff < -10 
   if ObjDers > 1
      if isempty(Prob.HessIx) & ~isempty(Prob.HessPattern)
         [ix,iy,iv]       = find(Prob.HessPattern);
         iz               = find(ix <= iy);
         Prob.HessPattern = sparse(ix(iz),iy(iz),iv(iz),Prob.N,Prob.N);
         Prob.HessIx      = findpatt(Prob.HessPattern);
      end
   end
   if Prob.NumDiff < -14 
      Prob.NumDiff = -1; 
   else
      Prob.NumDiff = min(-1,Prob.NumDiff+10);
   end
elseif Prob.NumDiff ~= 0
   if ObjDers > 1 & ~isempty(Prob.HessPattern)
      % Must ensure that HessPattern is upper triangular, if given
      [ix,iy,iv]       = find(Prob.HessPattern);
      iz               = find(ix <= iy);
      Prob.HessPattern = sparse(ix(iz),iy(iz),iv(iz),Prob.N,Prob.N);
      % Prob.HessIx      = findpatt(Prob.HessPattern);
   end
end

% MODIFICATION LOG
%
% 031201  hkh  Change to use double(0) and double(1)
% 031201  hkh  ADD AD handling
% 031203  hkh  Use checkMAD och checkADMAT to check if AD installed
% 040101  hkh  Check AD, and change in Prob, use input ObjDers, ConsDers
% 040102  hkh  Do presolve, if set, in this routine
% 040124  hkh  Test if optParam field is correct, otherwise display error text
% 040206  hkh  Better error text if MAD and ADMAT not installed properly
% 040402  hkh  Use fields TIME0 and TIME1 in Prob, instead of globals
% 040406  hkh  Add call to findpatt to compute Prob.ConIx, if ConsDiff > 6
% 040411  hkh  Set HessPatten, if d2LPattern set
% 040411  hkh  Set d2LPattern using maked2cPatt(ConsPattern) & HessPattern
% 040413  hkh  Call estConsPattern or estHessPattern if LargeScale
% 040414  hkh  Call estJacPattern, change order of code, 1st define NARG
% 040425  hkh  Add input smallA, if 1 detect and remove small A elements
% 040603  ango Message related to Prob.smallA changed. 
% 040828  ango Calls xnargout, for safety with function handles.
