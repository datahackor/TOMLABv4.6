% TOMLAB BAR/SQP QP Interface
% 
% Don't call this routine directly. Use barqpTL and sprqpTL instead.
%

function Result = bosqpTL(Prob, Solver);

if nargin < 2
   error(['barqpTL needs two inputs, the Prob structure and the ' ...
          'solver to run.']);
end

if isempty(Prob.QP.F) | isempty(find(Prob.QP.F))
   pType = checkType('lp');
else
   pType = checkType('qp');
end

Prob   = iniSolve(Prob,pType,0,0);
Result = ResultDef(Prob);
Result.Solver = 'BARQP';
Result.SolverAlgorithm = 'Sparse Barrier QP Code';

[lwr, upr, n, m1, m2] = defblbu(Prob, Inf, 1);

if isempty(Prob.QP.c)
  Prob.QP.c = zeros(n,1);
end


% Bounds and starting point
xl = lwr(1:n);
xu = upr(1:n);
x0 = DefPar(Prob,'x_0',(xl+xu)/2);

x0idx = find(xl <= -inf | xu >= inf);
x0(x0idx) = 0;

% Constraint bounds 
bl = lwr(n+1:n+m1);
bu = upr(n+1:n+m1);

%PriLev = DefPar(Prob,'PriLevOpt',1);
PriLev = DefPar(Prob,'PriLev',0);

BOS       = DefPar(Prob,'BOS',[]);
options   = DefPar(BOS,'options',[]);
moremem   = DefPar(BOS,'moremem',[0,0]);


validlist = {'CONTOL', 'OBJTOL', 'PGDTOL', 'MAXNFE', 'NITMAX', ...
             'NITMIN', 'IT1MAX', 'ALFLWR', 'ALFUPR', 'LYNFNC', ...
             'LYNPLT', 'LYNPNT', 'LYNVAR', ...             
             'IOFLAG', 'IOFLIN', 'IOFMFR', 'IOFPAT', ...
             'MAXLYN', 'TOLFIL', 'TOLKTC', 'TOLPVT', ...
             'ALGOPT', 'KTOPTN', 'IPOSTO'};
if(Solver<10)
  appvalidlist = {'NEWTON', 'IRELAX', 'BIGCON', 'FEATOL', ...
                  'PMULWR', 'PTHTOL', 'RHOLWR', 'IMAXMU', 'MUCALC', ...
                  'MXQPIT', 'IOFSHR'};
else
  appvalidlist = {'SLPTOL', 'SFZTOL', 'IOFSHR', 'IOFSRC', ...
                  'JACPRM', 'QPOPTN'};
end

validlist = {validlist{1:end}, appvalidlist{1:end}};

bosopt = tom2bosopt(options, validlist);

%
% Assumptions in MEX call:
%
% The MEX handles dense or sparse Prob.QP.F matrices
% automatically. 
%
% The constraint matrix must be sparse on input.  
% 

if pType == checkType('lp') % Linear
   H = [];
   
   [ier,x_k,f_k,g_k,Ax,lamv,lamc,statv,statc,cnd,needed] = bos(...
       Solver, Prob, BOS, bosopt, PriLev,...
       n,m1,xl,xu,x0, ...
       sparse(Prob.A), bl, bu,...
       H, full(Prob.QP.c), ...
       moremem);
   
else % QP
   
  [ier,x_k,f_k,g_k,Ax,lamv,lamc,statv,statc,cnd,needed] = bos(...
      Solver, Prob, BOS, bosopt, PriLev,...
      n,m1,xl,xu,x0, ...
      sparse(Prob.A), bl, bu,...
      Prob.QP.F, full(Prob.QP.c), ...
      moremem);
   
end

% Solution data 
Result.x_k = x_k;
Result.f_k = f_k;
Result.g_k = g_k;
Result.v_k = [lamv(:) ; lamc(:)];
Result.Ax = Ax;

[txt,flg] = ExitText(ier);
Result.ExitText = txt;
Result.ExitFlag = flg;
Result.Inform   = ier;

Result = StateDef(Result, x_k, Ax, [], ...
   Prob.optParam.xTol, Prob.optParam.bTol, [], lwr, upr, 1);

% Condition number, states
Result.BOS.cndnum  = cnd;
Result.BOS.statc   = statc;
Result.BOS.statv   = statv;

% Memory consumption (useful if re-solving??)
Result.BOS.moremem = needed;

Result = endSolve(Prob,Result);


% -------------------------------------------------
% Exit text and flag depending on 'ier' value
% -------------------------------------------------

function [T,F] = ExitText(ier)

F = [];

switch(ier)
   case 0,   T='Normal termination'; F = 0;
   case 101, T='Weak solution found'; F = 0;
   case 102, T=['Number of equality constraints = number of free ' ...
                'variables and ALGOPT ~= F'];
   case 104, T='Maximum number of function evaluations'; F = 1;
   case 105, T=['Small step in optimization phase; suboptimal feasible ' ...
                'point found'];
   case 106, T='Maximum number of iterations'; F = 1;
   case 108, T='Feasible point not found'; F = 4;
   case 109, T='Maximum number of interval halves in line search', ...
             F = 1;
   case 110, T=['Hessian diagonal reached its maximum value; suboptimal ' ...
                'feasible point found'];
   case 111, T=['Projected gradient calculation failed; constraints ' ...
                'may be degenerate'];
   case 112, T='Calculation of first order multiplier estimates failed; constraints may be degenerate';
   case 113, T='Suboptimal feasible point found';
   case 114, T='Barrier NLP failed with unexpected error';
   case 115, T='CONTOL > OBJTOL; convergence tolerances may be inappropriate';
   case 116, T='Uphill direction detected in line search';
   case 117, T='Reduced objective function is linear';
   case 118, T='One or more constraints were ignored';
      
      % Some (but not all) of the negative errors are also interesting
      
   case -103, T='QPOPTN incorrect value'; F = 10;
   case -104, T='|IHESHN| > 3'; F = 10;
   case -105, T='Max. number of iterations is too small'; F = 10;
   case -106, T='Number of line search steps is less than one';
              F = 10;
   case -109, T='ALGOPT invalid value'; F = 10;
   case -110, T='OBJTOL too small'; F = 10;
   case -111, T='PGDTOL too small'; F = 10;
   case -112, T='CONTOL too small'; F = 10;
   case -113, T='Constraint status invalid input'; F = 10;
   case -114, T='Invalid Jacobian data'; F = 10;
   case -115, T='Invalid Jacobian column value'; F = 10;
   case -116, T='Invalid Jacobian row value'; F = 10;
   case -117, T='Constraint bound(s) crossed over'; F = 10;
   case -124, T='Variable bound(s) crossed over'; F = 10;
   case -129, T=['Function error at initial point or during gradient ' ...
                 'evaluation']; F = 10;
   case -130, T='BIGCON < CONTOL'; F = 10;
   case -133, T=['Rank deficient Jacobian detected on successive ' ...
                 'iterations.'];
   case -134, T='Invalid option passed to solver.'; F = 10;
   case -136, T='FEATOL < CONTOL'; F = 10;
   case -137, T='Print files not closed properly.';
   case -138, T='NEWTON ~= 0, 1, 2; invalid NEWTON input'; F = 10;
   case -146, T='Option PMULWR too small'; F = 10;
   case -147, T='Unexpected error'; F = 10;
   case -148, T='|IRVCOM| > 1'; F = 10;
   case -149, T='PTHTOL too small'; F = 10;
   case -150, T='RHOLWR too small'; F = 10;
   case -152, T='IMAXMU < 1'; F = 10;
   case -153, T='I/O Error, too little disk space?'; F = 10;
   case -154, T='TOLKTC < 1'; F = 10;
   case -155, T='TOLPVT < 0 or TOLKTC > 0.5'; F = 10;
   case -156, T='IRELAX < 0 or IRELAX > 2'; F = 10;
   case -157, T='MUCALC = 0 or |MUCALC| > 3'; F = 10;
   case -158, T='MXQPIT < 1'; F = 10;

    
   case -1001, T=['(NDIM < 1); the number of variablesisles sth an ' ...
                  'one.']; F = 10;
   case -1002, T=['(MCON < 0); the number of constraints is ' ...
                  'negative.']; F = 10;
   case -1003, T=['Warm or Hot Start not preceded by a cold start. ' ...
                  'ISTART must be positive on first call.']; F = 10;
   case -1004, T=['(JSTRH(1) <= 0) or (JSTRH(k) >= JSTRH(k+1)); invalid ' ...
                  'input for Hessian column start array.']; F = 10;
   case -1005, T=['(NONZH <= 0) or (NONZH > NDNSH) or (NONZH > NZHDIM) ' ...
                  'where NONZH =JSTRH(NDIM+1)-1, and NDNSH =NDIM*(NDIM+1)/2; ' ...
                  'the number of Hessian nonzeros is less than one, ' ...
                  'exceeds the size of a dense Hessian, or is not ' ...
                  'consistent with the input dimension NZHDIM.']; F = 10;
   case -1006, T=['(IROWH(i) <= IDIAG) or (IROWH(i) > NDIM) where IDIAG ' ...
                  '= row index of diagonal; invalid value for Hessian ' ...
                  'row index. Hessian must be lower triangular.']; F = 10;
   case -1007, T=['Invalid input for Jacobian column index array. ' ...
                  'When ISTART = 1, valid input requires( JCOLA(1) = ' ...
                  '1) and (JCOLA(i) <= JCOLA(i+1)). When ISTART = 2, ' ...
                  'valid input requires (JCOLA(i) > 0) and (JCOLA(i) ' ...
                  ' <= NDIM).']; F = 10;
   case -1008, T=['(NONZA <= 0) or (NONZA > NDIM*MCON) or (NONZA > ' ...
                  'NZADIM) where NONZA = JCOLA(NDIM+1)-1; the number ' ...
                  'of Jacobian nonzeros is either less than one, ' ...
                  'exceeds the size of a dense matrix, or is ' ...
                  'inconsistent with the input dimension NZADIM.']; ...
                  F = 10;
   case -1009, T=['(IROWA(i) <= 0) or (IROWA(i) > MCON); invalid input ' ...
                  'for Jacobian row index.']; F = 10;
   case -1010, T=['(BUPR(i) < BLWR(i)) or (BUPR(i) >= BIGBND) and ' ...
                  '(BLWR(i) <= -BIGBND); constraint upper bound is ' ...
                  'less than lower bound, or both bounds are ' ...
                  'infinite.']; F = 10;
   case -1011, T=['(XUPR(i) < XLWR(i)); the variable upper bound is ' ...
                  'less than the lower bound.']; F = 10;
   case -1012, T=['(MEQUAL + NFIXVR > NDIM) where MEQUAL is the ' ...
                  'number of equality constraints, and NFIXVR is ' ...
                  'the number of fixed variables; there are no ' ...
                  'degrees of freedom.']; F = 10;
   case -1013, T=['(IPC < 0) and (IPU > 0); invalid output control ' ...
                  'flag.']; F = 10;
   case -1014, T=['NHOLD too small; insufficient real storage detected ' ...
                  'in HDSQSH. The required storage is specified in ' ...
                  'NEEDED.']; F = 10;
   case -1016, T=['A row of the constraint matrix is structurally ' ...
                  'zero.']; F = 10;
   case -1017, T=['Unexpected error in CNVRTR.'];
   case -1018, T=['Unexpected error in DCNVRT.'];
   case -1019, T=['NIHOLD too small; insufficient integer storage ' ...
                  'detected in HDSQSH. The required storage is specified ' ...
                  'in NEEDED.']; F = 10;
   case -1020, T=['Inconsistent linear constraints.']; F = 10;
   case -1021, T=['Incorrect value for KTOPTN. Defaultsm ay need to ' ...
                  'be set.']; F = 10;
   case -1103, T=['Singular KT matrix detected by multifrontal algorithm ' ...
                  '(IER = -503).'];
   case -1104, T=['Incorrect inertia of KT matrix.'];
   case -1105, T=['Unexpected error from multifrontal algorithm.'];
   case -1106, T=['Condition number of KT matrix isto o large (IRETCD ' ...
                  '= 12 or 14).'];
   case -1107, T=['Incorrect inertia of KT matrix after active set ' ...
                  'change with warm start (IRETCD = 11).'];
   case -1108, T=['Excessive fill during numeric factorization of KT ' ...
                  'matrix presumably because of ill-conditioning. A ' ...
                  'larger value of TOLFIL may be required (see ' ...
                  'subroutine HHSNLP).']; F = 10;
   case -1109, T=['Incorrect inertia for Schur-complement with relaxed ' ...
                  'feasibility tolerance.'];
   case -1111, T=['I/O error, probably caused by insufficient disk ' ...
                  'space.'];
 
   otherwise,
      if(ier<0) 
         T=sprintf('Status code %d',ier);
      else
         T=sprintf('Unknown solver status value %d',ier);
      end
      
end

if(ier<0)
   T = ['Error: ' T];
end


% ExitFlag, not implemented. 
if isempty(F)
  F = -1;
end

% MODIFICATION LOG
%
% 041123 ango Wrote file
% 041202 hkh  Clean up, skip return, MOD LOG at the end of the file
% 041202 hkh  Revise call to defblbu and StateDef
% 050124 frhe Renamed from barqpTL.m to bosqpTL.m
% 050124 frhe Modified to call the bos.dll-file, and added spr
%             error codes.


