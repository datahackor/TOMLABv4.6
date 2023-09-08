%
% Routine to check preset value in Prob and set all undefined values in Prob. 
%
% function Prob = ProbCheck(Prob, Solver, solvType, probType);
%
% INPUT:
%  Prob     Problem structure
%  Solver   Solver name
%  solvType Solver type (or optType, if solver not known)
%  probType Problem type
%
% OUTPUT:
%  Prob     Problem structure
%               
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written Aug 5, 1999.    Last modified May 6, 2004.
%

function Prob = ProbCheck(Prob, Solver, solvType, probType)

if isfield(Prob,'CHECK')
   if Prob.CHECK==1
      % ProbCheck has already checked this structure
      return
   end
end

if nargin < 4
   probType=[];
   if nargin < 3
      solvType=[];
   end
end

if isfield(Prob,'probType'), probType=Prob.probType; end

if isempty(solvType) & isfield(Prob,'solvType')
   solvType=Prob.solvType;
end
% Guess that solvType the same as probType if not set
if isempty(solvType), solvType=probType; end
% Guess that solvType is constrained if not either probType or solvType set
if isempty(solvType), solvType=3; end

Prob.solvType=solvType;


% Guess that probType the same as solvType if not set
if isempty(probType),        probType=solvType; end
Prob.probType=probType;

if ~isfield(Prob,'LineParam')
   Prob.LineParam=LineParamDef;
else
   Prob.LineParam=LineParamSet(Prob.LineParam);
end




% Check the struct and define any undefined items

if ~isfield(Prob,'A'),            Prob.A=[]; end
if ~isfield(Prob,'ADObj'),        Prob.ADObj=double(0); end
if ~isfield(Prob,'ADCons'),       Prob.ADCons=double(0); end
% Safeguard AD: ADMAT only handles option 2, not -2
if Prob.ADObj  < -1, Prob.ADObj  = 0; end
if Prob.ADCons < -1, Prob.ADCons = 0; end

if ~isfield(Prob,'b_L'),          Prob.b_L=[]; end
if ~isfield(Prob,'b_U'),          Prob.b_U=[]; end
if ~isfield(Prob,'c_L'),          Prob.c_L=[]; end
if ~isfield(Prob,'c_U'),          Prob.c_U=[]; end
if ~isfield(Prob,'CheckNaN'),     Prob.CheckNaN=0; end
if ~isfield(Prob,'cName'),        Prob.cName=[]; end
if ~isfield(Prob,'cols'),         Prob.cols=[]; end
if ~isfield(Prob,'ConIx'),        Prob.ConIx=[]; end
if ~isfield(Prob,'ConsDiff'),     Prob.ConsDiff=double(0); end
if ~isfield(Prob,'ConsIdx'),      Prob.ConsIdx=[]; end
if ~isfield(Prob,'ConsPattern'),  Prob.ConsPattern=[]; end
if ~isfield(Prob,'d2LPattern'),   Prob.d2LPattern=[]; end
if ~isfield(Prob,'f_Low'),        Prob.f_Low=-realmax; end
if ~isfield(Prob,'f_opt'),        Prob.f_opt=[]; end
if ~isfield(Prob,'GradTolg'),     Prob.GradTolg=[]; end
if ~isfield(Prob,'GradTolH'),     Prob.GradTolH=[]; end
if ~isfield(Prob,'GradTolJ'),     Prob.GradTolJ=[]; end
if ~isfield(Prob,'HessIx'),       Prob.HessIx=[]; end
if ~isfield(Prob,'HessPattern'),  Prob.HessPattern=[]; end
if ~isfield(Prob,'JacIx'),        Prob.JacIx=[]; end
if ~isfield(Prob,'JacPattern'),   Prob.JacPattern=[]; end
if ~isfield(Prob,'LargeScale'),   Prob.LargeScale=0; end
if ~isfield(Prob,'MaxCPU'),       Prob.MaxCPU=inf; end
if ~isfield(Prob,'MENU'),         Prob.MENU=0; end
if ~isfield(Prob,'Mode'),         Prob.Mode=2; end
Name=[];
if isfield(Prob,'Name'), Name = Prob.Name; end
if isempty(Name),        Name = 'User Problem 1'; end
Prob.Name=Name;
if ~isfield(Prob,'nState'),       Prob.nState=double(1); end
if ~isfield(Prob,'NumDiff'),      Prob.NumDiff=double(0); end
% N is set below
P=[];
if isfield(Prob,'P'), P = Prob.P; end
if isempty(P),        P = double(1); end
Prob.P=P;
if ~isfield(Prob,'plotLine'),     Prob.plotLine=0; end
if ~isfield(Prob,'PriLev'),       Prob.PriLev=double(1); end
if ~isfield(Prob,'PriLevOpt'),    Prob.PriLevOpt=double(0); end
if ~isfield(Prob,'probFile'),     Prob.probFile=[]; end
% probType set above
if ~isfield(Prob,'rows'),         Prob.rows=[]; end
if ~isfield(Prob,'simType'),      Prob.simType=double(0); end
if ~isfield(Prob,'smallA'),       Prob.smallA=double(0); end
if ~isfield(Prob,'SolverDLP'),    Prob.SolverDLP=[]; end
if ~isfield(Prob,'SolverFP'),     Prob.SolverFP=[]; end
if ~isfield(Prob,'SolverLP'),     Prob.SolverLP=[]; end
if ~isfield(Prob,'SolverQP'),     Prob.SolverQP=[]; end
if ~isfield(Prob,'uP'),           Prob.uP=[]; end
if ~isfield(Prob,'uPName'),       Prob.uPName=[]; end
if ~isfield(Prob,'WarmStart'),    Prob.WarmStart=0; end
if ~isfield(Prob,'Warning'),      Prob.Warning=1; end
if ~isfield(Prob,'x_0'),          Prob.x_0=[]; end
if ~isfield(Prob,'x_L'),          Prob.x_L=[]; end
if ~isfield(Prob,'x_U'),          Prob.x_U=[]; end
if ~isfield(Prob,'x_min'),        Prob.x_min=[]; end
if ~isfield(Prob,'x_max'),        Prob.x_max=[]; end
if ~isfield(Prob,'x_opt'),        Prob.x_opt=[]; end
if ~isfield(Prob,'xName'),        Prob.xName=[]; end
if ~isfield(Prob,'PartSep'),      Prob.PartSep=[]; end
if ~isfield(Prob,'Solver'),       Prob.Solver=[]; end
if ~isfield(Prob.Solver,'Alg'),   Prob.Solver.Alg=double(0); end
if ~isfield(Prob.Solver,'Name'),  Prob.Solver.Name=Solver; end
if ~isfield(Prob.Solver,'Method'),Prob.Solver.Method=double(0); end

if ~isfield(Prob,'QP'),           Prob.QP=[]; end
if ~isfield(Prob,'MIP'),          Prob.MIP=[]; end
if ~isfield(Prob,'GO'),           Prob.GO=[]; end
if ~isfield(Prob,'CGO'),          Prob.CGO=[]; end
if ~isfield(Prob,'ExpFit'),       Prob.ExpFit=[]; end
if ~isfield(Prob,'NTS'),          Prob.NTS=[]; end
if ~isfield(Prob,'USER'),         Prob.USER.f=[]; end
if ~isfield(Prob.USER,'fc'),      Prob.USER.fc=[]; end
if ~isfield(Prob,'DUNDEE'),       Prob.DUNDEE=[]; end
if ~isfield(Prob,'PENOPT'),       Prob.PENOPT=[]; end

if ~isfield(Prob,'SOL') 

    Prob.SOL = struct( 'SpecsFile', [], 'PrintFile', [], 'SummFile', [], ...
        'xs', [], 'hs', [], 'nS', double(0), 'hElastic', [], ...
        'iState', [], 'cLamda', [], 'R', [], 'optPar', -999*ones(1,63), ...
        'optParN', 63); 

   %Prob.SOL.SpecsFile = []; % File to read SPECS parameters from
   %Prob.SOL.PrintFile = []; % Log Print File
   %Prob.SOL.SummFile  = []; % Summary File (standard output, but not in Matlab)
   %Prob.SOL.start     = 'Cold';
   %% xs,hs,nS is used for warm start of SNOPT and SQOPT
   %Prob.SOL.xs        = [];
   %Prob.SOL.hs        = [];
   %Prob.SOL.nS        = 0;
   %Prob.SOL.hElastic  = [];
   %% iState is used for warm start of QPOPT, NPSOL
   %Prob.SOL.iState    = [];
   %% cLamda,R is used for warm start of NPSOL
   %Prob.SOL.cLamda    = [];
   %Prob.SOL.R         = [];
   %% Vector with optimization parameters. 
   %% The meaning of each index is described in the m-file interface routines:
   %% minos.m (minosTL.m), snopt.m sqopt.m qpopt.m npsol.m etc.
   %Prob.SOL.optPar(1:63) = -999; 
   %Prob.SOL.optParN   = 63; % Number of parameters in optPar 
%else
   %if ~isfield(Prob.SOL,'optParN'),   Prob.SOL.optParN=63; end
   %if ~isfield(Prob.SOL,'SpecsFile'), Prob.SOL.SpecsFile=[]; end
   %if ~isfield(Prob.SOL,'SummFile'),  Prob.SOL.SummFile=[]; end
end


if checkType('qp',[probType;solvType]) | ... 
   checkType('lp',[probType;solvType]) | ...
   checkType('mip',[probType;solvType])
   if ~isfield(Prob.QP,'F'),         Prob.QP.F=[]; end
   if ~isfield(Prob.QP,'c'),         Prob.QP.c=[]; end
   if ~isfield(Prob.QP,'B'),         Prob.QP.B=[]; end
   if ~isfield(Prob.QP,'y'),         Prob.QP.y=[]; end
   if ~isfield(Prob.QP,'Q'),         Prob.QP.Q=[]; end
   if ~isfield(Prob.QP,'R'),         Prob.QP.R=[]; end
   if ~isfield(Prob.QP,'E'),         Prob.QP.E=[]; end
   if ~isfield(Prob.QP,'Ascale'),    Prob.QP.Ascale=[]; end
   if ~isfield(Prob.QP,'DualLimit'), Prob.QP.DualLimit=[]; end
   if ~isfield(Prob.QP,'UseHot'),    Prob.QP.UseHot=[]; end
   if ~isfield(Prob.QP,'HotFile'),   Prob.QP.HotFile=[]; end
   if ~isfield(Prob.QP,'HotFreq'),   Prob.QP.HotFreq=[]; end
   if ~isfield(Prob.QP,'HotN'),      Prob.QP.HotN=[]; end
end

if checkType('mip',[probType;solvType])
   if ~isfield(Prob.MIP,'IntVars'),   Prob.MIP.IntVars=[]; end
   if ~isfield(Prob.MIP,'VarWeight'), Prob.MIP.VarWeight=[]; end
   if ~isfield(Prob.MIP,'fIP'),       Prob.MIP.fIP=[]; end
   if ~isfield(Prob.MIP,'xIP'),       Prob.MIP.xIP=[]; end
   if ~isfield(Prob.MIP,'PI'),        Prob.MIP.PI=[]; end
   if ~isfield(Prob.MIP,'SC'),        Prob.MIP.SC=[]; end
   if ~isfield(Prob.MIP,'SI'),        Prob.MIP.SI=[]; end
   if ~isfield(Prob.MIP,'sos1'),      Prob.MIP.sos1=[]; end
   if ~isfield(Prob.MIP,'sos2'),      Prob.MIP.sos2=[]; end
   if ~isfield(Prob.MIP,'xpcontrol'), Prob.MIP.xpcontrol=[]; end
   if ~isfield(Prob.MIP,'callback'),  Prob.MIP.callback=[]; end
   if ~isfield(Prob.MIP,'KNAPSACK'),  Prob.MIP.KNAPSACK=0; end
end


if checkType('exp',probType)
   if ~isfield(Prob.ExpFit,'p'),       Prob.ExpFit.p=[]; end
   if ~isfield(Prob.ExpFit,'wType'),   Prob.ExpFit.wType=[]; end
   if ~isfield(Prob.ExpFit,'eType'),   Prob.ExpFit.eType=double(1); end
   if ~isfield(Prob.ExpFit,'infCR'),   Prob.ExpFit.infCR=[]; end
   if ~isfield(Prob.ExpFit,'dType'),   Prob.ExpFit.dType=[]; end
   if ~isfield(Prob.ExpFit,'geoType'), Prob.ExpFit.geoType=[]; end
   if ~isfield(Prob.ExpFit,'qType'),   Prob.ExpFit.qType=[]; end
   if ~isfield(Prob.ExpFit,'sigType'), Prob.ExpFit.sigType=[]; end
   if ~isfield(Prob.ExpFit,'lambda'),  Prob.ExpFit.lambda=[]; end
   if ~isfield(Prob.ExpFit,'alpha'),   Prob.ExpFit.alpha=[]; end
   if ~isfield(Prob.ExpFit,'beta'),    Prob.ExpFit.beta=[]; end
   if ~isfield(Prob.ExpFit,'x0Type'),  Prob.ExpFit.xoType=[]; end
   if ~isfield(Prob.ExpFit,'sumType'), Prob.ExpFit.sumType=[]; end
end


%if checkType('nts',probType)
%   Prob.NTS.SepAlg     =SepAlg;
%   Prob.NTS.ntsModel   =ntsModel;
%   Prob.NTS.p          =p;
%   Prob.NTS.pL         =pL;
%   Prob.NTS.pA         =pA;
%   Prob.NTS.LambdaScale=LambdaScale;
%   Prob.NTS.ntsSeed    =ntsSeed;
%   Prob.NTS.N          =N;
%   Prob.NTS.t1         =t1;
%   Prob.NTS.tN         =tN;
%   Prob.NTS.gamma      =gamma;
%   Prob.NTS.alphaArt   =alphaArt;
%   Prob.NTS.lambdaArt  =lambdaArt;
%
%   Prob.NTS.lambda=lambda;
%   Prob.NTS.alpha=alpha;
%end

if checkType('glb',probType) | checkType('glc',probType)
   % Global Optimization
   if ~isfield(Prob,'N')
      Prob.N=max(length(Prob.x_L),length(Prob.x_U));
   end
else
   Prob.N=max([length(Prob.x_L),size(Prob.A,2),length(Prob.x_0), ...
               length(Prob.x_U)]);
end

if isfield(Prob,'mLin') 
   mLin    = Prob.mLin;
else
   Prob.mLin=[]; 
end
if isfield(Prob,'mNonLin')
   mNonLin = Prob.mNonLin;
else
   Prob.mNonLin=[]; 
end

N       = Prob.N;
if isempty(mLin)
   mLin = size(Prob.A,1);
   Prob.mLin = mLin;
end
if isempty(mNonLin)
   mNonLin = max(length(Prob.c_L),length(Prob.c_U));
   Prob.mNonLin = mNonLin;
end

M = mLin + mNonLin;

if ~isfield(Prob,'LS'), Prob.LS=[]; end

% Cases LS, LLS, CLS, EXP, NTS
if any( checkType('ls',probType)  | checkType('lls',probType)  ...
      | checkType('cls',probType) | checkType('exp',probType)  ...
      | checkType('nts',probType))
   if ~isfield(Prob.LS,'weightType'), Prob.LS.weightType=double(0); end
   if ~isfield(Prob.LS,'weightY'),    Prob.LS.weightY=[]; end
   if ~isfield(Prob.LS,'t'),          Prob.LS.t=[]; end
   if ~isfield(Prob.LS,'y'),          Prob.LS.y=[]; end
   if ~isfield(Prob.LS,'C'),          Prob.LS.C=[]; end
   if ~isfield(Prob.LS,'damp'),       Prob.LS.damp=[]; end
   if ~isfield(Prob.LS,'L'),          Prob.LS.L=[]; end
   if ~isfield(Prob.LS,'yUse'),       Prob.LS.yUse=1; end
   if ~isfield(Prob.LS,'SepAlg'),     Prob.LS.SepAlg=0; end


   % If empty USER.f and USER.fc, assume least squares problem

   if isempty(Prob.USER.f) & isempty(Prob.USER.fc)
      Prob.USER.f = 'ls_f';
      Prob.USER.g = 'ls_g';
      if checkType('ls',probType)
         Prob.USER.H = 'lls_H';
      else
         Prob.USER.H = 'ls_H';
      end
   end
end

if ~isfield(Prob.USER,'g'),   Prob.USER.g=[]; end
if ~isfield(Prob.USER,'H'),   Prob.USER.H=[]; end
if ~isfield(Prob.USER,'c'),   Prob.USER.c=[]; end
if ~isfield(Prob.USER,'dc'),  Prob.USER.dc=[]; end
if ~isfield(Prob.USER,'d2c'), Prob.USER.d2c=[]; end
if ~isfield(Prob.USER,'r'),   Prob.USER.r=[]; end
if ~isfield(Prob.USER,'J'),   Prob.USER.J=[]; end
if ~isfield(Prob.USER,'d2r'), Prob.USER.d2r=[]; end
if ~isfield(Prob.USER,'gdc'), Prob.USER.gdc=[]; end

if size(Prob.A,1) > 0 & N ~= size(Prob.A,2)
   fprintf('\n\n');
   fprintf('Number of variables to optimize %d\n',N);
   fprintf('Rows in matrix A                %d\n',size(Prob.A,1));
   fprintf('Columns in matrix A             %d\n',size(Prob.A,2));
   error('ProbCheck: Columns in A must match number of variables')
end
if length(Prob.b_L) ~= size(Prob.A,1) & length(Prob.b_U) ~= size(Prob.A,1)
   fprintf('\n\n');
   fprintf('Length of b_L    %d\n',length(Prob.b_L));
   fprintf('Length of b_U    %d\n',length(Prob.b_U));
   fprintf('Rows in matrix A %d\n',size(Prob.A,1));
   error('ProbCheck: One of lengths of b_L and b_U must match rows in A')
end

if ~isfield(Prob,'optParam')
   Prob.optParam = optParamDef(Solver,probType,N,N,M);
else
   % Check the optParam structure
   Prob.optParam = optParamSet(Prob.optParam,Solver,probType,N,N,M);
end

Prob.CHECK=1;

% MODIFICATION LOG:
%
% 990812 hkh  Expand QP fields.
% 990824 hkh  Expand QP fields with solver selections.
% 990909 hkh  Adding HessPattern, JacPattern and ConsPattern
% 000830 hkh  Fix for SOL subfields
% 000922 hkh  Correct the logic for global optimization, fields must be defined
% 001014 hkh  Added test on ConsDiff
% 010903 hkh  Now 63 parameters in optPar
% 011110 hkh  Add two new fields, GO and RBF, for global optimization
% 011112 hkh  Bug in definition of Prob.N, if a global opt problem defined as
%             being of another type
% 011213 hkh  Changed Prob.MIP.SC into SC, SI, semi-continuous and semi-integer
% 011226 hkh  Changed Prob.MIP.SOS1 and SOS2 to sos1 and sos2.
%             Initialize Prob.MIP.KNAPSACK = 0; not as empty
% 020103 hkh  Change field RBF to CGO
% 020409 hkh  Add field Mode and nState
% 020512 hkh  Check field DUNDEE for Fletcher and Leyffer Dundee solvers
% 020630 hkh  Check field PENSDP for Tomlab /PENSDP
% 021216 hkh  Add damp and L in field LS
% 030117 hkh  Change field PENSDP to PENOPT, for /PENSDP and /PENBMI
% 030117 hkh  Add check on field CheckNaN, default 0
% 030522 hkh  Add check on field USER.fc and in check of least squares
% 030524 hkh  Add check on field USER.gdc. Use double(0),double(1)
% 031129 hkh  Change fields AutoDiff to ADObj, ADCons
% 031201 hkh  Safeguard for ADObj <-1, ADCons <-1, ADMAT do not handle -2 
% 040102 hkh  Compute fields mLin,mNonLin; Add checks on linear constraints
% 040115 hkh  Avoid working with mLin, mNonLin, N twice, set fields of missing
% 040115 hkh  Set fields mLin and mNonLin if missing
% 040303 hkh  Move definition of Prob.USER.fc outside if-then-end-block
% 040425 hkh  Add field Prob.smallA + more, change order more like ProbDef
% 040506 hkh  conIx should be ConIx
% 040728 med  Keyboard removed
