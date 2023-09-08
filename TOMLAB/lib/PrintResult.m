%
% Print results, stored in the Result structure
%
% function PrintResult(Result,PriLev);
% 
% Called from the driver routine, tomRun.m.
% May also be called from command line, if Result and Result.Prob is defined
%
% ............................................................................
% PriLev - Printing level from 0-3
% PriLev = 0   Silent
% PriLev = 1   Problem Number and Name. f_*, f_0(if computed), f_opt (if given)
% PriLev = 2   Optimal x, Starting value x_0, 
%              ExitFlag, Inform, number of m-file evaluations  
% PriLev = 3   Lagrange multipliers, both returned and TOMLAB estimate, 
%              TOMLAB computes reduced gradient (projected gradient)
%              Distance from start to solution
%              LS Residual, gradient and projected gradient.
% PriLev = 4   LS Jacobian, Hessian, Quasi-Newton Hessian approximation
%
% ............................................................................
% If Result.Prob.PrintLM = 0, avoid TOMLAB computation of Lagrange multipliers
% and reduced (projected) gradient (may take some CPU time)
% 
% If Result.Prob.PrintLM = 1 or not set, compute Lagrange multipliers
% and reduced Hessian, and display statistics about constraint violation
%
% See help LagMult.m for the results computed
% ............................................................................
%
% To avoid too many variables, constraints and residuals in the output,
% three global variables are limiting the number printed:
%
% MAX_x   Maximal number of variables
% MAX_c   Maximal number of constraints
% MAX_r   Maximal number of residuals in least squares problems
%
% Example: To increase the number of variables to 50, do
%
% global MAX_x
% MAX_x = 50;
% PrintResult(Result);
%              
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomlab.biz
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.3.0$
% Written May 19, 1998.   Last modified Apr 7, 2004.
%

function PrintResult(Result,PriLev);

if nargin < 2
   PriLev=3; % Maximal printing level - 1
end
if PriLev < 1, return; end

global n_f n_g n_H    % Count of function evals, gradient evals, Hessian evals
global n_c n_dc n_d2c % Count of constraint evals, constraint gradient evals
global n_J n_r n_d2r  % Count of Jacobian evals, residual evals, 2nd der evals
global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

if isfield(Result,'Prob')
   Prob       = Result.Prob;
else
   if isempty(Result)
      fprintf('\nPrintResult: Empty Result. Nothing to print!\n\n');
   else
      Result
      fprintf('\n\n\nPrintResult: System Error!!! ');
      fprintf('Result.Prob not defined\n\n\n');
      fprintf('Can not print results.\n');
   end
   return
end

probType   = Prob.probType;
if isfield(Prob,'solvType')
   solvType  = Prob.solvType;
else
   solvType = Result.solvType;
end
Solver     = Result.Solver;
wait       = Prob.optParam.wait;
probNumber = Prob.P;
xTol       = Prob.optParam.xTol;
PrintLM = DefPar(Prob,'PrintLM',1);

x_0   = Result.x_0;
x_k   = Result.x_k;
% x_0
% x_k
% pause
n     = length(x_k);
v_k   = Result.v_k;
f_k   = Result.f_k;
f_0   = Result.f_0;
g_k   = Result.g_k;
r_k   = Result.r_k;

f_opt = Prob.f_opt;
x_opt = Prob.x_opt;

% Avoid computing r_k. If user is using extra params, it will not work.
%if probType > 3 & probType < 8
%   r_k = nlp_r(x_k, Prob);
%end
if any(size(x_opt)==1)
   x_opt=x_opt(:)';  % make it a row vector
end

if PriLev > 1, fprintf('\n'); end

   % x_k could be empty if no feasible point is found in glcSolve
   fprintf('===== * * * =================================================================== * * *\n');
   fprintf('%s',tomlablic(10));
   fprintf('\n');
   fprintf('=====================================================================================\n');
   s='Problem:';
   if ~isempty(Prob.probFile)
      if isstr(Prob.probFile)
         s2=deblank(Prob.probFile);
         ss=sprintf('%%%ds -',1+length(s2));
         s=[s sprintf(ss,s2)];
      else
         %s=[s ' No Init File    -'];
         s=[s ' ---'];
      end
   end 
   s=[s sprintf(' %2.0f: %s',probNumber,Prob.Name)];
   if length(s) < 46
      s=[s blanks(46-length(s))];
   end
   fprintf('%s',s);
   if abs(f_k) > 1E30 
      fprintf('  f_k %26.18e',f_k);
   else
      fprintf('  f_k %26.18f',f_k);
   end
   fprintf('\n');
   if ~isempty(f_opt)
      for i = 1:length(f_opt)
          fprintf('%s',blanks(length(s)-11)); 
          fprintf('User given f(x_*)%26.18f\n',f_opt(i)); 
      end
   end 
if ~isempty(x_k) 
   if size(Prob.A,1) > 0 
      Ax=DefPar(Result,'Ax',[]);
      if isempty(Ax) & size(x_k,1) == size(Prob.A,2)
         Ax=Prob.A*x_k(:,1); 
      end
   else
      Ax=zeros(0,1);
   end
   c_k = Result.c_k;

   if isempty(c_k)
      zz=[x_k(:,1);Ax]; 
      bl=[Prob.x_L;Prob.b_L];
      bu=[Prob.x_U;Prob.b_U];
   else
      zz=[x_k(:,1);Ax;c_k(:,1)]; 
      mL=length(c_k);
      if isempty(Prob.c_L)
         bl=[Prob.x_L;Prob.b_L;-Inf*ones(mL,1)];
      else
         bl=[Prob.x_L;Prob.b_L;Prob.c_L];
      end
      if isempty(Prob.c_U)
         bu=[Prob.x_U;Prob.b_U;Inf*ones(mL,1)];
      else
         bu=[Prob.x_U;Prob.b_U;Prob.c_U];
      end
   end
   zz = real(zz);
   if length(bl)==length(zz) & length(bu)==length(zz) 
      h_k = norm(min(0,zz-bl),1) + norm(min(0,bu-zz),1);    
   else
      fprintf('Conflicting dimensions on bounds and constraint values\n');
      fprintf('Length(bl) %d Length(bu) %d Length([x,Ax,c_k]) %d\n',...
               length(bl),length(bu),length(zz));
      NN=min([length(bl),length(bu),length(zz)]);
      h_k = norm(max(0,bl(1:NN)-zz(1:NN)),1)+norm(max(0,zz(1:NN)-bu(1:NN)),1);
   end
   if ~isempty(h_k) & ~(h_k==0)
      fprintf('%s',blanks(length(s)-16));
      fprintf('         sum(|constr|)%26.18f',h_k);
      fprintf('\n');
      fprintf('%s',blanks(length(s)-16));
      fprintf('f(x_k) + sum(|constr|)%26.18f',f_k+h_k);
      fprintf('\n');
   end

   if ~isempty(f_0) & PriLev > 1
      fprintf('%s',blanks(length(s)));
      fprintf('f(x_0)%26.18f\n',f_0);
   end
   vL=sum(Prob.x_L > x_k(:,1)+xTol);
   vU=sum(Prob.x_U < x_k(:,1)-xTol);
   if vL+vU > 0
      fprintf(' ==>  Number of variables violating lower bound%4d. ',vL);
      fprintf(' Number of variables violating upper bound%4d',vU);
      fprintf('\n');
   end
end

   if ~isempty(Prob.uP) & PriLev > 1
      fprintf('User parameters uP:\n');
      z=Prob.uP(1:min(length(Prob.uP),15));
      disp(num2str(z(:)'));
   end

   if PriLev > 0
      fprintf('\nSolver: %s',Solver);
      if isstr(Result.ExitFlag)
         fprintf('.  EXIT=%s',Result.ExitFlag);
      else
         fprintf('.  EXIT=%d',Result.ExitFlag);
      end
      if ~isempty(Result.Inform)
         if isstr(Result.Inform)
            fprintf('.  INFORM=%s', Result.Inform);
         else
            fprintf('.  INFORM=%d', Result.Inform);
         end
      end
      fprintf('.\n');
      for i = 1:size(Result.SolverAlgorithm,1)
         fprintf('%s',Result.SolverAlgorithm(i,:));
         fprintf('\n');
      end
    
      for i = 1:size(Result.ExitText,1)
         fprintf('%s',Result.ExitText(i,:));
         fprintf('\n');
      end
     if strcmpi(Solver,'glcFast') | strcmpi(Solver,'glcSolve') | ...
         strcmpi(Solver,'glbFast') | strcmpi(Solver,'glbSolve') ...
     else
      NumDiff = Prob.NumDiff;

      if any(Prob.ADObj == [-2 2]) | any(Prob.ADCons == [-2 2])
         fprintf('ADMAT TB Automatic Differentiation estimating:');   
         if Prob.ADObj == 2 
            fprintf(' gradient, Hessian');   
            if Prob.ADCons == 2
               fprintf(' and');   
            end
         end
         if Prob.ADCons == 2
            fprintf(' constraint gradient');   
         end
         fprintf('\n');   
      end
      if any(Prob.ADObj == [-1 1]) | any(Prob.ADCons == [-1 1])
         fprintf('MAD TB Automatic Differentiation estimating:');   
         if Prob.ADObj ==1 
            fprintf(' gradient');   
         elseif Prob.ADObj ==-1
            fprintf(' Hessian');   
         end
         if any(Prob.ADObj == [-1 1]) & any(Prob.ADCons == [-1 1])
            fprintf(' and');   
         end
         if Prob.ADCons ==1
            fprintf(' constraint gradient');   
         elseif Prob.ADCons ==-1
            fprintf(' constraint Hessian');   
         end
         fprintf('\n');   
      end
      if NumDiff > 0
         fprintf('Gradient estimated with ');
         switch NumDiff
         case 1
            fprintf('finite differencing');
         case 2
            fprintf('Matlab spline routine ');
         case 3
            fprintf('spline routine csaps:');
            fprintf(' Smoothing  %f',Prob.optParam.splineSmooth);
         case 4
            splineTol = Prob.optParam.splineTol;
            if splineTol < 0
               fprintf('spline routine csapi:');
            else
               fprintf('spline routine spaps:');
               fprintf(' Tolerance %f',splineTol);
            end
         case 5
            fprintf('complex variable method');
         case 6
            fprintf('internal solver method');
         end
         fprintf('\n');
      elseif NumDiff < 0
         if strcmpi(Solver,'minlpbb') | strcmpi(Solver,'filterSQP') ...
            | strcmpi(Solver,'conSolve') | strcmpi(Solver,'nlpSolve') ...
            | strcmpi(Solver,'knitro') | strcmpi(Solver,'conopt') ...
            | strcmpi(Solver,'ucSolve') 
            fprintf('Hessian estimated with ');
            switch abs(NumDiff)
            case 1
               fprintf('finite differencing');
            case 2
               fprintf('Matlab spline routine ');
            case 3
               fprintf('spline routine csaps:');
               fprintf(' Smoothing  %f',Prob.optParam.splineSmooth);
            case 4
               splineTol = Prob.optParam.splineTol;
               if splineTol < 0
                  fprintf('spline routine csapi:');
               else
                  fprintf('spline routine spaps:');
                  fprintf(' Tolerance %f',splineTol);
               end
            case 5
               fprintf('complex variable method');
            case 6
               fprintf('internal solver method');
            end
            fprintf('\n');
         end
      end
      ConsDiff=Prob.ConsDiff;
      if ConsDiff > 0
         fprintf('Constraint Jacobian estimated with ');
         switch ConsDiff
         case 1
            fprintf('finite differencing');
         case 2
            fprintf('Matlab spline routine ');
         case 3
            fprintf('spline routine csaps:');
            fprintf(' Smoothing  %f',Prob.optParam.splineSmooth);
         case 4
            splineTol = Prob.optParam.splineTol;
            if splineTol < 0
               fprintf('spline routine csapi:');
            else
               fprintf('spline routine spaps:');
               fprintf(' Tolerance %f',splineTol);
            end
         case 5
            fprintf('complex variable method');
         case 6
            fprintf('internal solver method');
         end
         fprintf('\n');
      elseif ConsDiff < 0
         if strcmpi(Solver,'minlpbb') | strcmpi(Solver,'filterSQP') ...
            | strcmpi(Solver,'nlpSolve') ...
            | strcmpi(Solver,'knitro') | strcmpi(Solver,'conopt')
            fprintf('Hessian of constraints estimated with ');
            switch abs(ConsDiff);
            case 1
               fprintf('finite differencing');
            case 2
               fprintf('Matlab spline routine ');
            case 3
               fprintf('spline routine csaps:');
               fprintf(' Smoothing  %f',Prob.optParam.splineSmooth);
            case 4
               splineTol = Prob.optParam.splineTol;
               if splineTol < 0
                  fprintf('spline routine csapi:');
               else
                  fprintf('spline routine spaps:');
                  fprintf(' Tolerance %f',splineTol);
               end
            case 5
               fprintf('complex variable method');
            case 6
               fprintf('internal solver method');
            end
            fprintf('\n');
         end
      end
      fprintf('\n');
     end
   end
     
   if PriLev > 0
      if checkType('ls',solvType)  | checkType('cls',solvType) | ...
         checkType('lls',solvType) | checkType('exp',solvType)
         fprintf('ResEv  %4.0f JacEv  %4.0f ',n_r,n_J);
         if checkType('cls',solvType)  | checkType('exp',solvType)
             if ~isempty(Result.ConstrEv)
                fprintf(' ConstrEv %4.0f ',Result.ConstrEv)
             end
         end
      else
         if any([Result.FuncEv,Result.GradEv,Result.ConstrEv] > 0)
            if ~isempty(Result.FuncEv)
               fprintf('FuncEv %4.0f ',Result.FuncEv);
            end
            if ~isempty(Result.GradEv)
               fprintf('GradEv %4.0f ',Result.GradEv);
            end
            if Result.HessEv > 0
               fprintf('HessEv %4.0f ',Result.HessEv);
            end
            if checkType('qp',solvType)  | checkType('con',solvType) | ...
               checkType('glc',solvType) | checkType('minlp',solvType)
               if ~isempty(Result.ConstrEv)
                  fprintf('ConstrEv %4.0f ',Result.ConstrEv)
               end
            end
            if Result.ConJacEv > 0
               fprintf('ConJacEv %4.0f ',Result.ConJacEv);
            end
            if Result.ConHessEv > 0
               fprintf('ConHessEv %4.0f ',Result.ConHessEv);
            end
         end
      end
      %if isempty(n_f), n_f=Result.FuncEv; end
      %if isempty(n_g), n_g=Result.GradEv; end
      %if isempty(n_c), n_c=Result.ConstrEv; end
      %if (n_f~=Result.FuncEv) | (n_g~=Result.GradEv) | ...
      %   (n_c~=Result.ConstrEv & n_c > 0)
      %   if n_f+n_c+n_g > 0
      %      if any([Result.FuncEv,Result.GradEv,Result.ConstrEv] > 0)
      %         fprintf('\nTOMLAB Global Variable Counters give:\n'); 
      %      end
      %      if solvType == 9 | solvType == 10 % Global solver
      %         fprintf(' FuncEv %4.0f ',n_f);
      %      else
      %         fprintf(' FuncEv %4.0f GradEv %4.0f',n_f,n_g)
      %      end
      %      if (solvType==2 | solvType ==3 | solvType==5 | solvType == 6) ...
      %         fprintf(' ConstrEv %4.0f ConGradEv %4.0f',n_c,n_dc);
      %      end
      %   end
      %end
      if Result.Iter > 0, fprintf('Iter %4.0f ',Result.Iter); end
      if Result.MinorIter > 0, fprintf('MinorIter %4.0f',Result.MinorIter); end
      fprintf('\n');
   end
   if PriLev > 0 & any(solvType == [ 9 10])
      if ~isempty(Result.maxTri)
         fprintf('Maximum rectangle size %25.15f\n',Result.maxTri); 
      end
      if isfield(Result,'minPnts')
         if ~isempty(Result.minPnts)
            fprintf('Optimal points found at f(x) evaluation:'); 
           
            if length(Result.minPnts) > 10, fprintf('\n'); end
            for i = 1:length(Result.minPnts)
                fprintf(' %d ',Result.minPnts(i)); 
                if mod(i,15) == 0, fprintf('\n'); end
            end
            fprintf('\n');
         end
      end
   end

   if PriLev > 0
      priline=0;
      if Result.CPUtime > 0
         priline=1;
         fprintf('CPU time: %f sec. ',Result.CPUtime);
      end
      if Result.REALtime > 0
         priline=1;
         fprintf('Elapsed time: %f sec. ',Result.REALtime);
      end
      if priline
         fprintf('\n');
      end
   end


   if PriLev > 1 & ...
      ~any([checkType('glb',solvType),checkType('glc',solvType)]) ...
      & ~isempty(x_0)
      fprintf('Starting vector x');
      if length(x_0) > MAX_x
         fprintf(' (Length %d. Print first %d)',length(x_k), MAX_x);
      end
      fprintf(':\n');
      xprint(x_0(1:min(length(x_0),MAX_x)),'x_0:');
   end
if ~isempty(x_k) 
   if PriLev > 1
      if checkType('glb',solvType) | checkType('glc',solvType)
         fprintf('Optimal vector(s) x');
         if size(x_k,2) > MAX_r 
            fprintf('.  Number of minima %d - display %d of them', ...
               size(x_k,2),min(size(x_k,2),MAX_r));
         end
         if size(x_k,1) > MAX_x 
            fprintf('.   Length(x) %d - display %d first elements', ...
               size(x_k,1),min(size(x_k,1),MAX_x));
         end
         fprintf(':\n');
         mPrint(x_k(1:min(size(x_k,1),MAX_x),1:min(size(x_k,2),MAX_r))',...
                'x_k',' %13.6f',5);
      else
         fprintf('Optimal vector x');
         if length(x_k) > MAX_x
            fprintf(' (Length %d. Print first %d)',length(x_k),MAX_x);
         end
         fprintf(':\n');
         xprint(x_k(1:min(length(x_k),MAX_x)),'x_k:');
      end
      if ~isempty(x_opt)
         for j=1:size(x_opt,1)
             if size(x_opt,2) > n
                minPnt=x_opt(j,n+1); 
             else
                minPnt=0;
             end
             fprintf('User given stationary point x_* (%d):',j);
             if minPnt==0
                fprintf(' - Minimum point'); 
             elseif minPnt==1
                fprintf(' - Saddle point');
             elseif minPnt==2
                fprintf(' - Maximum point');
             end
             fprintf('\n');
             xprint(x_opt(j,1:min([n,size(x_opt,2),MAX_x])),'x_*:');
         end
      end
   end
   if PriLev > 2 &  ...
      (checkType('qp',solvType) | checkType('con',solvType) | ...
       checkType('cls',solvType) | checkType('lp',solvType)) 
      fprintf('Lagrange multipliers v. Vector length %d',length(v_k));
      if length(v_k) > MAX_c
         fprintf(' (Print first %d)',MAX_c);
      end
      fprintf(':\n');
      xprinte(v_k(1:min(length(v_k),MAX_c)),'v_k:');
   end
   if PriLev > 3 & ...
      ~any([checkType('glb',solvType),checkType('glc',solvType)]) & ...
      ~isempty(x_0) & ~isempty(x_0)
      if length(x_k)==length(x_0)
         x_diff=x_k-x_0;
         fprintf('Diff x-x0:\n');
         xprinte(x_diff(1:min(length(x_diff),MAX_x)),'    ');
      end
   end
   if checkType('ls',solvType)  | checkType('cls',solvType) | ...
      checkType('lls',solvType) | checkType('exp',solvType)
      %% Must recompute residual in the correct point
      %r_k = nlp_r(x_k, Prob); 
      %Result.J_k = nlp_J(x_k, Prob); 
      %g_k = Result.J_k'*r_k;

      if n_f<=0, n_f=n_r; end
      if n_g<=0, n_g=n_J; end
   end

%
% Compute the Hessian
%
   if PriLev > 2 & (checkType('ls',solvType)  | checkType('cls',solvType) | ...
         checkType('lls',solvType) | checkType('exp',solvType))
      fprintf('Residual r_k');
      if length(r_k) > MAX_r
         fprintf('(Length %d. Print %d first elements)',length(r_k),MAX_r);
      end
      fprintf(':\n');
      xprinte(r_k(1:min(length(r_k),MAX_r)),'r_k:');

      global wLS
      r_U=r_k;
      if ~isempty(wLS)
         if any(wLS~=1)
            xprinte(wLS(1:min(length(wLS),MAX_r)),'w_*:');
            fprintf('UNWEIGHTED Residual r_U: \n');
            ix=find(wLS~=0);
            r_U(ix)=r_k(ix)./wLS(ix);
            xprinte(r_U(1:min(length(r_k),MAX_r)),'r_U:');
         end
         if ~isempty(Prob.LS.y)
            r_R=nan*ones(length(r_U),1);
            ix=find(Prob.LS.y~=0);
            if ~isempty(ix)
               r_R(ix)=r_U(ix)./Prob.LS.y(ix);
               fprintf('Relative error in unweighted residual, r_R \n');
               fprintf('(y(t)==0 ==> NaN):\n');
               xprinte(r_R(1:min(length(r_k),MAX_r)),'r_R:');
            end
         end
      end
   end
   if PriLev > 2 & ~any(solvType == [7 9 10 11])
      % Must recompute f_k, in case of global vars used 
      % (if not initialized as empty). Applies to ls and cls probType.
      %if solvType < 4 | solvType > 6
      %   f_k = nlp_f(x_k, Prob);
      %   g_k = nlp_g(x_k, Prob);
      %end
      if ~isempty(g_k)
         fprintf('Gradient g_k');
         if length(g_k) > MAX_x
            fprintf(' Largest abs(gradient)%30.17f',...
                      max(abs(g_k)));
            fprintf(' (Print first %d)',MAX_x);
         end
         fprintf(':\n');
         xprinte(g_k(1:min(length(g_k),MAX_x)),'g_k:');
      end
   end
   %if (probType==3 | probType==2 | probType==6 | probType==8) & PriLev > 1
   if PrintLM
    if ((PriLev > 2 & n < 200) | (PriLev > 3 & n < 1000)) ...
      & ~any(solvType == [5 7 9 11])

      [v_k, Zv, P, Z_L, cErr, ceq, cineq, gProj] = LagMult(Prob,Result);

      m=length(c_k);

      if m > 0
         cTol=Prob.optParam.cTol;

         if length(c_k) > MAX_c
            fprintf('Number of nonlinear constraints = %d',length(c_k));
            fprintf('. Print the %d first error deviations.\n',MAX_c);
         end
         xprinte(cErr(1:min(size(c_k,1),MAX_c)),'cErr');
         if ceq > 0
            fprintf('%d equalities off more than cTol = %15.3e\n',ceq,cTol);
         end
         if cineq > 0
            fprintf('%d inequalities off more than cTol = %15.3e\n',...
                     cineq,cTol);
         end
         if ceq > 0 | cineq > 0 
            [m1,i1]= max(abs(cErr));
            fprintf('Worst constraint validation = %15.3e',m1);
            fprintf(' for constraint # %d\n',i1);
         end
      end

      % Display projected gradient, should get descent

      if size(Z_L,2)~=0 & ~isempty(g_k)
         % If [U S V]=svd(Z_L); then gProj=(eye(n)-U*S*V'*V*pinv(S)*U')*g_k;
         % if Z_L = Q * R * E'; then gProj=g_k-Q*Q'*g_k;
         
         fprintf('Projected gradient gPr: ');
         if length(gProj) > MAX_x
            fprintf(' Largest abs(Projected gradient)%30.17f',...
                     max(abs(gProj)));
            fprintf(' (Print first %d)',MAX_x);
         end
         fprintf(':\n');
         xprinte(gProj(1:min(length(gProj),MAX_x)),'gPr:');
      else
         gProj=g_k;
      end

      j=find(abs(gProj) > 1E-5);
      if ~isempty(j)
         fprintf('*** WARNING: %d reduced gradient values > 1E-5 ***',...
                 length(j));
         fprintf(' Worst value: %e \n',max(abs(gProj)));
      end

      % Lagrange multipliers

      % v = Z_L \ g_k
      m2=min(size(Z_L));
      if m2 < size(Z_L,2)
         fprintf('TOMLAB found %d active constraints. ',size(Z_L,2))
         fprintf('The Lagrange multipliers are not unique.')
         if ~isempty(g_k)
            fprintf('\n');
            fprintf('Estimate Lagrange multipliers for the %d first',m2);
            fprintf(' and fix the others as zero.');
         end
         fprintf('\n');
      end
      if m2 > 0 
         if ~isempty(g_k) & ~any(probType==[2 7 8])
            jj=6;
            fprintf('    ');
            ii=0;
            for i=2:min(MAX_x,size(Zv,1))
                if ii == jj 
                   fprintf('\n    ');
                   ii=0;
                end
                ii=ii+1;
                fprintf('%s',Zv(i,:));
                fprintf(' ');
            end
            fprintf('\n');
            ix = find(P~=0);
            xprinte(v_k(ix(1:min(MAX_c,length(ix)))),'v:');
            %if S_inv(m2)==0 & ~any(solvType==[4 5 6 9 10])
            %   disp('Warning: Rank problems estimating Lagrange multipliers')
            %end
         end
      elseif size(Prob.A,1) > 0 | m > 0
         fprintf('TOMLAB found no active constraints.\n');
      end
    end
   end
   if PriLev > 3 & ~isempty(Prob.USER.H) & ...
      ~any(solvType==[4 5 6 8 9 10]) & n < 1000
      if isempty(Result.H_k)
         %
         % Try to compute the Hessian
         %
         Result.H_k = nlp_H(x_k, Prob);
      end
      if ~isempty(Result.H_k) 
         fprintf('Eigenvalues of Hessian at x_k\n');
         lambda=eig(full(Result.H_k));
         if min(lambda) < 1E-5
            xprinte(lambda(1:min(length(lambda),MAX_x)),'eig:');
         else
            xprint(lambda(1:min(length(lambda),MAX_x)),'eig:');
         end
      end
   end
   if PriLev > 3 & ~isempty(Result.B_k) & ~any(solvType == [8 9 10]) & ...
      isempty(Result.H_k)
      fprintf('\n');
      fprintf('The Quasi-Newton Hessian approximation\n');
      fprintf('\n');
      mPrint(Result.B_k(1:min(length(x_k),MAX_x),1:min(length(x_k),MAX_x)),...
            'B_k',' %10.6e',6);
      if PriLev > 3 & solvType ~= 9 & solvType ~= 10
          if wait, pause; end
      end
   end
   if PriLev > 3  & ...
      (checkType('ls',solvType)  | checkType('cls',solvType) | ...
       checkType('lls',solvType) | checkType('exp',solvType))
       fprintf('\n');
       [M,N] = size(Result.J_k);
       fprintf('The Jacobian J_k');
       if M > MAX_r | N > MAX_x
          fprintf('(Dim %d by %d. Print %d by %d)',M,N,MAX_r, MAX_x);
       end
       fprintf('\n');
       mPrint(Result.J_k(1:min(M,MAX_r),1:min(N,MAX_x)), 'J_k',' %13.6e',5);
   end
   if PriLev > 2  & ...
      (checkType('ls',solvType)  | checkType('cls',solvType) | ...
       checkType('lls',solvType) | checkType('exp',solvType))
       fprintf('\n');
       N = size(Result.J_k,2);
       if issparse(Result.J_k)
          %fprintf('The structural rank of the Jacobian J_k');
          %fprintf(' is %d out of %d columns\n',sprank(J_k),N);
       else
          J_k = Result.J_k;
          if any(isinf(J_k)) | any(isnan(J_k))
             disp('Some NaN or Inf element found in J_k');
             J_k
          else
             fprintf('The rank of the Jacobian J_k');
             fprintf(' is %d out of %d columns\n',rank(J_k),N);
          end
       end
   end
end
if PriLev > 2
   fprintf('\n');
   fprintf('=== * * * ================================================== * * *\n');
   fprintf('\n');
end
if wait, pause; end

function [h_k,cErr] = ConstraintError(v1,v2,alg)

cErr=max(v1,v2);

if alg < 2
   h_k = norm(max(0,cErr),Inf);  % Constr. violation, max norm
else
   h_k = norm(max(0,cErr),1);    % Constr. violation, L1 - sum of abs
end

% MODIFICATION LOG:
%
% 980916  hkh  Change from optPar to structure optParam
% 980917  hkh  Use Result.ExitFlag and other exit results.
% 981013  hkh  Changed the logic. PrintResult may be called from command line
% 981015  hkh  Fixed display of constraint errors. Terrible paranthesis bug:
%              ceq+abs(cErr(i)) > cTol; must be ceq+(abs(cErr(i)) > cTol);
% 981023  hkh  Added n_d2r in global
% 981026  hkh  Added projected gradient and check this, not ordinary gradient
% 981028  hkh  Introduce f_0 = f(x_0). Write another line with f_0.
% 981102  hkh  Print the user given f_* = f(x_*) if nonempty.
%              Less info written for PriLev==1.
% 981107  hkh  Check if Result.Inform is empty, or string.
%              Check if Result.ExitFlag is a string
% 981109  hkh  Print v_k if solvType==6, linear constraints are present
%              Change print levels to range 0-3.
% 981110  hkh  Change to MAX_x,MAX_c,MAX_r. Print H_k separate from B_k
% 981113  hkh  Add printing of v_k for LP
% 981116  hkh  Print f_k + sum of infeasibilities. Avoid reducing the global
%              counters, gives wrong answer! Change to f_k as optimal point
%             
% 981126  hkh  Add printing about the use of automatic/numeric differentiation
% 981207  hkh  Add printing of user given stationary points
% 990204  hkh  Add printing of rank of Jacobian
% 990304  hkh  More information about residuals: unweighted, relative errors
% 990312  hkh  Avoid problems when wLS is empty, due to Matlab "features"
% 990408  mbk  & solvType ~= 9 added at several places.
%              x_k changed to x_k(:,1) on line 78.
%              if PriLev > 0 & ~isempty(x_k) on line 60.
% 990810  hkh  Revised. Use LagMult, also for projected gradient.
% 000906  hkh  Print empty init file, if probFile is missing
% 000916  hkh  Add print of ExitText with solver result
% 001107  hkh  Use cTol for nonlinear feasibility test
% 010321  hkh  Output of worst constraint validation and projected gradient
% 011031  hkh  Output of maxTri (max triangle size) for global optimization
% 020111  hkh  Use exponential format for f(x) if abs(f(x)) > bound
% 020412  hkh  Write indices on many lines if many global solutions
% 020417  hkh  Print index for worst constraint violation
% 020619  hkh  Avoid calling LagMult if too large problem (n > 1000)
% 020702  hkh  If Result empty, print a new message
% 020922  hkh  Avoid some printing if solvType == 11 (miqp)
% 030110  hkh  No eig(Hessian) if n >= 1000
% 030113  hkh  Add printing for ConsPattern < 0 (Hessian of constraints)
% 030211  hkh  Avoid printing text about Hessians when solver not using them
% 030309  hkh  Avoid copies of Hessian and Jacobians
% 030309  hkh  Using Result.Ax for A*x if computed
% 030309  hkh  Prob.PrintLM = 0 => no Lagrange multipliers computation, no
%              printing of projected Hessian and other constraint information
% 030323  hkh  Correcting Result.J_k handling
% 030516  hkh  Change Reduced gradient to Projected gradient
% 030610  ango Change tomlab(10) to tomlablic(10)
% 031201  hkh  Change AutoDiff to ADObj, ADCons, add MAD TB for AD
% 031202  hkh  Fixed erroneous ||, not OK in Matlab 6.1 and before
% 031204  hkh  ADObj == -1 should be ADCons == -1 in constrained Hessian test
% 031204  hkh  max does not work if some element complex, sometimes h_k wrong
% 031206  hkh  Revision of AD printing for both MAD and ADMAT
% 040126  hkh  Safe guard A*x computation
% 040215  hkh  Check if PriLev <= 0 and return if so
% 040323  hkh  strcmpi('minlp', should be strcmpi('minlpbb'
% 040331  hkh  Split output into one additional print level
% 040331  hkh  Enable constraint output for type glc
% 040404  hkh  Avoid gradient output if glbFast,glcFast,glbSolve,glcSolve
% 040407  hkh  Revise NumDiff,ConsDiff-output, remove PriLev > 0

