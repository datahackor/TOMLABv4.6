% TOMLAB BAR/SQP NLP Interface
% 
% Don't call this routine directly. Use barnlpTL and sqpnlpTL instead.
%

function Result=bosnlpTL(Prob, Solver)

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

if nargin < 1, error('bosnlpTL needs the Prob structure as input'); end

Prob.solvType = 3; % NLP (CON) solver

Prob = iniSolve(Prob,3,2,2);

Result = ResultDef(Prob);
if(Solver<10)
  Result.Solver='BARNLP';
  Result.SolverAlgorithm='MEX-interface to sparse BARNLP code';
else
  Result.Solver='SPRNLP';
  Result.SolverAlgorithm='MEX-interface to sparse SPRNLP code';
end

Result.Prob = Prob;

%PriLev = DefPar(Prob,'PriLevOpt',1);
PriLev = DefPar(Prob,'PriLev',0);

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

BIG=1E20;
[bl,bu,n,m1,m2] = defblbu(Prob, BIG);
m=m1+m2;

xl = bl(1:n);
xu = bu(1:n);
cl = bl(n+1:end);
cu = bu(n+1:end);

x_0=DefPar(Prob,'x_0',zeros(n,1));

% The row & col indices for the entire Jacobian, including linear constraints:
if m2>0 & isempty(Prob.ConsPattern)
   ConsPattern = spones( [Prob.A ; ones(m2,n) ] );
else
   ConsPattern = spones( [Prob.A ; Prob.ConsPattern ] );
end
[rowG,colG] = find( ConsPattern );
nonzG = length(rowG);

% Linear index from multiple subscripts for nonzero pattern (including linear matrix A)
if ~isempty(ConsPattern)
   Prob.ConsIdx = sub2ind(size(ConsPattern),rowG,colG);
else
   Prob.ConsIdx = [];  % Unconstrained 
end
 
   
% Hessian sparsity information and linear index, needed for callbacks. 
d2LPattern = tril(spones(DefPar(Prob,'d2LPattern',[])));
if isempty(d2LPattern)
   d2LPattern = tril(spones(ones(n,n)));
end
[rowH,colstH]=tomsol(31,d2LPattern,1);
nonzH = length(rowH);

% Should give the same rowH as above, but we need the column indices too. 
[rowH,colH]=find(d2LPattern);
Prob.d2LIdx = sub2ind(size(d2LPattern),rowH,colH);

% Solver specifics
BOS       = DefPar(Prob,'BOS',[]);
PrintFile = DefPar(BOS,'PrintFile','');
options   = DefPar(BOS,'options',struct([]));

validlist = {'CONTOL', 'OBJTOL', 'PGDTOL', 'MAXNFE', 'NITMAX', ...
             'NITMIN', 'IT1MAX', 'ALFLWR', 'ALFUPR', 'LYNFNC', ...
             'LYNPLT', 'LYNPNT', 'LYNVAR', ...             
             'IOFLAG', 'IOFLIN', 'IOFMFR', 'IOFPAT', ...
             'MAXLYN', 'TOLFIL', 'TOLKTC', 'TOLPVT', ...
             'ALGOPT', 'KTOPTN', 'IPOSTO'};
if(Solver<10)
  appvalidlist = {'NEWTON', 'IRELAX', 'BIGCON', 'FEATOL', ...
                  'PMULWR', 'PTHTOL', 'RHOLWR', 'IMAXMU', 'MUCALC', ...
                  'MXQPIT'};
else
  appvalidlist = {'SLPTOL', 'SFZTOL', 'IOFSHR', 'IOFSRC', ...
                  'JACPRM', 'QPOPTN'};
end

validlist = {validlist{1:end}, appvalidlist{1:end}};

bosopt = tom2bosopt(options, validlist);

%optPar    = DefPar(BOS,'optPar',NaN*ones(43,1));
%optParam  = DefPar(Prob,'optParam',[]);
%optPar    = optParBAR(options,optPar,optParam);

morereal = DefPar(BOS,'morereal',0);
moreint  = DefPar(BOS,'moreint',0);

moremem = [morereal moreint];

[Inform,x_k,f_k,c_k,rcv,rcc,nfev]=bos(Solver, Prob, BOS, bosopt, ...
                                      PriLev, ...
                                      n,m1,m2, ...
                                      xl,xu,x_0,cl,cu,nonzG,rowG,colG,...
                                      nonzH,rowH,colstH,moremem);

Result.v_k = [rcv;rcc];

Result.x_k = x_k;
Result.x_0 = x_0;
Result.f_k = f_k;

Result.Ax   = c_k(1:m1);
Result.c_k  = c_k(m1+1:m1+m2);

Result.cJac = nlp_dc(x_k,Prob);
global n_f n_g n_H n_c 
Result.FuncEv=n_f;
Result.GradEv=n_g;
Result.HessEv=n_H;
Result.ConstrEv=n_c;

switch(Inform)
   case {103,104,106}
      ExitFlag = 1; % Too many iterations (or function evals/errors)
   case 108
      ExitFlag = 4; % Infeasible
   otherwise
      ExitFlag = 0;
end
      
if Inform<0,ExitFlag = 10; end % Input errors 

switch(Inform)
   case {103,104,106}
      ExitFlag = 1; % Too many iterations (or function evals/errors)
   case 108
      ExitFlag = 4; % Infeasible
   otherwise
      ExitFlag = 0;
end
      
if Inform<0,ExitFlag = 10; end % Input errors 

% ExitText
switch(Inform)
   case 0,   Text='Normal termination.';
      % Less than optimal. 
   case 101, Text='Weak solution found (relaxation required and/or multipliers near zero).';
   case 102, Text='Number of equality constraints=number of free variables and ALGOPT ~= ''F''.';
   case 103, Text='Maximum number of consecutive function errors.';
   case 104, Text='Maximum number of function evaluations.';
   case 105, Text='Small step termination in optimization phase; suboptimal feasible point found.';
   case 106, Text='Maximum number of iterations.';
   case 108, Text='Feasible point not found.';
   case 109, Text='Maximum number of interval halves in line search.';
   case 110, Text='Hessian diagonal reached its maximum value; suboptimal feasible point found.';
   case 111, Text='Projected gradient calculation failed; constraints may be degenerate.';
   case 112, Text='Calculation of first order multiplier estimates failed; constraints may be degenerate.';
   case 113, Text='Suboptimal feasible point found.';
   case 114, Text='Barrier NLP failed with unexpected error.';
   case 115, Text='CONTOL > OBJTOL; convergence tolerances may be inappropriate.';
   case 116, Text='Uphill direction detected in line search.';
   case 117, Text='Reduced objective function is linear.';
   case 118, Text='Both cLi = -.01/HDMCON(5) and cUi = .01/HDMCON(5); constraints ignored.';
   case 119, Text='Terminate after diagnostic line search.';
   case 121, Text='Terminate after postoptimality analysis.';
   case 704, Text='Sparse or dense parameter defaults may cause poor algorithm performance.';
      % ERRORS
   case -101, Text='(MCON < 0); the number of constraints is negative.';
   case -102, Text='(NDIM < 1); the number of variables is less than one.';
   case -103, Text='Incorrect value for QPOPTN.';
   case -104, Text='|IHESHN| > 3.';
   case -105, Text='(NITMAX < max(1,NITMIN)); the maximum number of iterations is either less than 1 or less than the minimum number of iterations.';
   case -107, Text='(IT1MAX < 1); the number of line search steps is less than 1.';
   case -108, Text='(IOFLAG < 0) or (IOFLAG > 30); invalid input for the output control flag.';
   case -109, Text='Invalid input for ALGOPT.';
   case -110, Text='(OBJTOL <= 10*HDMCON(5)); objective function tolerance is too small.';
   case -111, Text='(PGDTOL <= [HDMCON(5)]^(1/2) or (PGDTOL > 10E?2); the projected gradient tolerance is too small or too large.';
   case -112, Text='CONTOL < [HDMCON(5)]^(1/2) ; the constraint tolerance is too small.';
   case -113, Text='(ISTATC(i) < 0) or (ISTATC(i) > 4); invalid input for constraint status.';
   case -114, Text='(NONZG <= 0) or (NONZG > NDIM*MCON); the number of Jacobian nonzeros is either less than one, or exceeds the size of a dense matrix.';
   case -115, Text='(JCOLG(i) <= 0) or (JCOLG(i) > NDIM); invalid input for Jacobian column index array.';
   case -116, Text='(IROWG(i) <= 0) or (IROWG(i) > MCON); invalid input for Jacobian row index.';
   case -117, Text='(CUPR(i) < CLWR(i)); constraint upper bound is less than lower bound.';
   case -118, Text='Either (CUPR(i) = CLWR(i)) and (ISTATC(i)!=3) or (CUPR(i)!=CLWR(i)) and (ISTATC(i) = 3); constraint status array is not consistent with bounds.';
   case -119, Text='(CUPR(i)!=CLWR(i)) and (|CUPR(i)-CLWR(i)| < CONTOL); constraint bounds are not equal, but di?er by less than the constraint tolerance.';
   case -120, Text='(ISTATV(i) < 0) or (ISTATV(i) > 3); invalid input for variable status array.';
   case -121, Text='(JSTRH(1)!=1) or (JSTRH(k)>=JSTRH(k+1)); invalid input for Hessian column start array.';
   case -122, Text='(NZHDIM < NDIM) or (NZHDIM > NDNSH) or (NZHDIM!=NONZH) where NZHDIM = JSTRH(NDIM+1)-1, and NDNSH = NDIM*(NDIM+1)/2; the number of Hessian nonzeros is less than one, exceeds the size of a dense Hessian, or is not consistent with the input value NONZH.';
   case -123, Text='(IROWH(i) <= 0) or (IROWH(i) > NDIM); invalid value for Hessian row index.';
   case -124, Text='(XUPR(i) < XLWR(i)); the variable upper bound is less than the lower bound.';
   case -125, Text='Either (XUPR(i) = XLWR(i)) and (ISTATV(i)!=3) or (XUPR(i)!=XLWR(i)) and (ISTATV(i) = 3); the variable status array is not consistent with the bound values.';
   case -126, Text='(XUPR(i) != XLWR(i)) and (|XUPR(i)-XLWR(i)| < CONTOL); variable bounds are not equal but differ by less than the constraint tolerance.';
   case -127, Text='Double precision hold array too small; insufficient storage detected in HDBNLP interface. The required storage is specified in NEEDED.';
   case -128, Text='Integer hold array too small; insufficient storage detected in HDBNLP interface. The required storage is specified in NEEDED.';
   case -129, Text='Function error at initial point or during gradient evaluation.';
   case -130, Text='BIGCON < CONTOL.';
   case -131, Text='NHOLD too small; insufficient double precision storage detected in algorithm. The required storage is specified in NEEDED.';
   case -132, Text='NIHOLD too small; insufficient integer storage detected in algorithm. The required storage is specified in NEEDED.';
   case -133, Text='Rank deficient Jacobian detected on successive iterations.';
   case -134, Text='HHSNLP input error; invalid character string displayed.';
   case -135, Text='(IFERR < 0) or (IFERR > 1); invalid value for function error flag.';
   case -136, Text='FEATOL < CONTOL; the initial feasibility tolerance is less than the constraint tolerance.';
   case -137, Text='Conflict between user and multifrontal file number; check IPUMF1,. . .,IPUMF6 IPUDRF, IPUFZF, IPUSTF.';
   case -138, Text='NEWTON != 0, 1, 2; invalid input for Newton method option flag.';
   case -145, Text='(IROWH(JSTRH(i)) != i) for some value of i; Hessian diagonal elements are incorrect.';
   case -146, Text='PMULWR < HDMCON(5).';
   case -147, Text='Unexpected error; check storage allocation.';
   case -148, Text='|IRVCOM| > 1.';
   case -149, Text='PTHTOL < [HDMCON(5)]^(1/2).';
   case -150, Text='RHOLWR < HDMCON(5).';
   case -151, Text='A constraint is inconsistent (either ci(x) < cLi or cUi < ci(x)) and cannot be changed because the corresponding row of the Jacobian is zero.';
   case -152, Text='IMAXMU < 1.';
   case -153, Text='I/O error, probably caused by insufficient disk space.';
   case -154, Text='TOLKTC < 1.';
   case -155, Text='TOLPVT < 0 or TOLPVT > .5.';
   case -156, Text='IRELAX < 0 or IRELAX > 2.';
   case -157, Text='MUCALC = 0 or |MUCALC| > 4.';
   case -158, Text='MXQPIT < 1.';
      %
   otherwise, Text='Unknown IERNLP return value';
end

optParam  = DefPar(Prob,'optParam',[]);
Result = StateDef(Result, x_k(1:n), Result.Ax, Result.c_k, ...
                  optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 1);

Result.ExitFlag = ExitFlag;
Result.ExitText = Text;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG
%
% 031023 ango Created file
% 040122 hkh  Revision for v4.2, add iniSolve and endSolve calls
% 050125 frhe Renamed to bosnlpTL.m from barnlpTL.m and modified to
%             call the new single bos.dll.
