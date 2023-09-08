% PrintAssign. Internal TOMLAB utility
%
% Do NOT call. Use lpAssign, qpAssign and probAssign

function [FName, Name] = PrintAssign(probType, problemName, setupFile, nProblem)

if isempty(setupFile)
   FName     = [];
   Name      = problemName;
   setupFile = [];
   return
end

if exist(setupFile)==2
   if nProblem(1) <= 1
      fprintf('Checking new setup file %s',setupFile');
      fprintf('\n');
      fprintf('The file already exists. Delete or rename the existing file.');
      fprintf(' Or use another name.\n')
      error('Can not create new setup file');
   else
      FileExists=1;
      P=nProblem(1);
   end
else
   P=1;
   FileExists=0;
end

%if isempty(findstr('MATCOM',version))

   [Path,setupName,Ext,Ver]=fileparts(setupFile);

   if isempty(Path)
       File=[setupName '.m'];
   else
       File=[Path filesep setupName '.m'];
   end
   setupFile=setupName;
%end

FName=sprintf('%s_P%d.mat',setupFile,P);

if nProblem(1)==1 & length(nProblem)==1
   Name=problemName;
else
   % Add a number for the problem
   Name=[problemName ' - ' num2str(nProblem(1))];
end

%save(FName,'setupFile','Name','F','c','A','b_L','b_U','x_0','x_L','x_U', ...
%     'x_min','x_max','x_opt','f_opt');

if P > 1
   % Only create the mat file with problem P
   return;
end

fid=fopen(File,'w');

if fid <= 0
   fprintf('Illegal file id %d when opening file %s',fid,File);
   error('Setup file procedure can not proceed!');
   return
end

fprintf(fid,'%s','%');
fprintf(fid,'\n');

fprintf(fid,'%s','%');
fprintf(fid,' TOMLAB User Defined ');

pCode = checkType([],probType);

switch pCode
  case 'uc'
    fprintf(fid,'Unconstrained Optimization');
  case 'qp'
    fprintf(fid,'Quadratic Programming');
  case 'con'
    fprintf(fid,'Constrained Optimization');
  case 'ls'
    fprintf(fid,'Nonlinear Least Squares');
  case 'lls'
    fprintf(fid,'Linear Least Squares');
  case 'cls'
    fprintf(fid,'Constrained Nonlinear Least Squares');
  case 'mip'
    fprintf(fid,'Mixed-Integer Linear Programming');
  case 'lp'
    fprintf(fid,'Linear Programming');
  case 'glb'
    fprintf(fid,'Box-Bounded Global Optimization');
  case 'glc'
    fprintf(fid,'Mixed-Integer Nonlinear Global Optimization');
  case 'exp'
    fprintf(fid,'Exponential Fitting');
end

fprintf(fid,' problems');

fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,' See xxx_prob.m for 1st part of code and comments');
fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,' Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz.');
fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,' Copyright (c) 1999-2004 by Tomlab Optimization Inc., Sweden. ');
fprintf(fid,'$Release 4.5.0$');
fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,' Written %s',date);
fprintf(fid,'.   Last modified %s',date);
fprintf(fid,'\n');
fprintf(fid,'%s','%');
fprintf(fid,'\n\n\n');

fprintf(fid,'function [probList, Prob] = %s',setupFile);
fprintf(fid,'(P, ask, Prob);');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'if nargin < 3 \n')
fprintf(fid,'   Prob=[]; \n')
fprintf(fid,'   if nargin < 2 \n')
fprintf(fid,'      ask=1; \n')
fprintf(fid,'      if nargin < 1 \n')
fprintf(fid,'         P=[]; \n')
fprintf(fid,'      end \n')
fprintf(fid,'   end \n')
fprintf(fid,'end \n')
fprintf(fid,'\n');
if nProblem(1) == 1 & length(nProblem)==1
   fprintf(fid,'   probList=str2mat(...\n');
   fprintf(fid,'      ''%s',problemName);
   fprintf(fid,'''...\n');
   fprintf(fid,'      ); '); 
   fprintf(fid,'%s','%');
   fprintf(fid,' TO MAKE NEW PROBLEM ENTRY.'); 
   fprintf(fid,' MAKE A COPY OF THE NEXT ROW:\n'); 
   fprintf(fid,'    ');
   fprintf(fid,'    %s','%');
   fprintf(fid,',''%s',problemName);
   fprintf(fid,'''...\n');
   fprintf(fid,'\n');
   fprintf(fid,'         %s','%');
   fprintf(fid,' CHANGE TO THE NEW NAME, REMOVE ');
   fprintf(fid,'%s','%');
   fprintf(fid,'.');
   fprintf(fid,' DO NOT REMOVE THE COMMA(,)!\n');
   fprintf(fid,'         %s','%');
   fprintf(fid,' MOVE THE ROW AFTER THE ');
   fprintf(fid,'LAST PROBLEM NAME ROW ABOVE. ');
   fprintf(fid,'%s','%');
   fprintf(fid,'\n');
elseif length(nProblem) > 1
   fprintf(fid,'nProblem=%d;',nProblem(2));
   fprintf(fid,'\n');
   fprintf(fid,'probList=[''%s',problemName);
   fprintf(fid,' - '' num2str(%d)];',nProblem(1));
   fprintf(fid,'\n');
   fprintf(fid,'\n');
   fprintf(fid,'for i=%d:nProblem',nProblem(1)+1);
   fprintf(fid,'\n');
   fprintf(fid,'   probList=str2mat(probList,[''%s',problemName);
   fprintf(fid,' - '' num2str(i)]);');
   fprintf(fid,'\n');
   fprintf(fid,'end');
   fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'if isempty(P)\n');
fprintf(fid,'   return\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'if ask==-1 & ~isempty(Prob)\n');
fprintf(fid,'   if isstruct(Prob)\n');
fprintf(fid,'      if ~isempty(Prob.P)\n');
fprintf(fid,'         if P==Prob.P & strcmp(Prob.probFile,''%s',setupFile);
fprintf(fid,'''), return; end\n');
fprintf(fid,'      end\n');
fprintf(fid,'   end\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');

fprintf(fid,'if isempty(Prob)\n');
fprintf(fid,'   Prob=ProbDef;\n');
fprintf(fid,'elseif Prob.P~=P & Prob.P > 0\n');
fprintf(fid,'   Prob=ProbDef;\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');

fprintf(fid,'global probType\n');
fprintf(fid,'probType=checkType(''%s'');\n',pCode);
fprintf(fid,'Prob.probType=probType;\n');     
fprintf(fid,'\n');

switch pCode
  case 'qp'
    %QP,LP
    fprintf(fid,'[x_0, x_L, x_U, x_min, x_max, f_Low, ');
    fprintf(fid,'xName, x_opt, f_opt,...\n');
    fprintf(fid,'      cName, A, b_L, b_U, F, c, B,  ...\n');
    fprintf(fid,'      pSepFunc, uP, uPName] = qpVarDef(Prob);\n');
  case {'mip'}
    fprintf(fid,'[x_0, x_L, x_U, x_min, x_max, f_Low, ');
    fprintf(fid,'xName, x_opt, f_opt,...\n');
    fprintf(fid,'      cName, A, b_L, b_U, c, B,  ...\n');
    fprintf(fid,'      xpcontrol, KNAPSACK, VarWeight, IntVars, ...\n');
    fprintf(fid,'      pSepFunc, uP, uPName] = mipVarDef(Prob);\n');
  case {'lp'}
    fprintf(fid,'[x_0, x_L, x_U, x_min, x_max, f_Low, ');
    fprintf(fid,'xName, x_opt, f_opt,...\n');
    fprintf(fid,'      cName, A, b_L, b_U, F, c, B,  ...\n');
    fprintf(fid,'      pSepFunc, uP, uPName] = qpVarDef(Prob);\n');
end

fprintf(fid,'\n');
fprintf(fid,'\n');

if nProblem == 1 & length(nProblem)==1
   fprintf(fid,'if P == 1\n');
   fprintf(fid,'   Name=''%s',problemName);
   fprintf(fid,''';\n');
   fprintf(fid,'   load(''%s',setupFile);
   switch pCode
     case {'qp','mip','lp'}
       %QP,MIP,LP
       fprintf(fid,'_P1'',''setupFile'',''Name'',''F'',''c'',''A'',');
       fprintf(fid,'''b_L'',''b_U'', ...\n');
       fprintf(fid,'       ''x_0'',''x_L'',''x_U'',');
       fprintf(fid,'''x_min'',''x_max'',''x_opt'',''f_opt'');\n');
   end

else
   fprintf(fid,'if P > 0 & P <= nProblem\n');
   fprintf(fid,'   FName=sprintf(''%s',setupFile);
   fprintf(fid,'_P');
   fprintf(fid,'%s','%');
   fprintf(fid,'d.mat'',P);\n');
   fprintf(fid,'   load(FName,');

   switch pCode
     case {'qp','mip','lp'}
       %QP,MIP,LP
      fprintf(fid,'''setupFile'',''Name'',''F'',''c'',''A'',');
      fprintf(fid,'''b_L'',''b_U'', ...\n');
      fprintf(fid,'       ''x_0'',''x_L'',''x_U'',');
      fprintf(fid,'''x_min'',''x_max'',''x_opt'',''f_opt'');\n');
   end
end
fprintf(fid,'\n');

fprintf(fid,'\n');

% End of problem definition

if nProblem == 1 & length(nProblem)==1

   fprintf(fid,'\n');
   fprintf(fid,'%s','%');
   fprintf(fid,'elseif P == 2\n');
   fprintf(fid,'%s','%');
   fprintf(fid,'ADD THE SECOND PROBLEM IN THE SAME WAY AS THE FIRST.\n');
   fprintf(fid,'%s','%');
   fprintf(fid,'REMOVE THE COMMENT ON THE ABOVE elseif P==2 LINE.\n');
   fprintf(fid,'\n');
   fprintf(fid,'%s','%');
   fprintf(fid,'REMEMBER TO ADD THE PROBLEM NAME IN THE probList\n');
   fprintf(fid,'%s','%');
   fprintf(fid,'DEFINITION ON ROW 17-23.\n');
   fprintf(fid,'\n');
   fprintf(fid,'%s','%');
   fprintf(fid,'ADD MORE PROBLEMS IN THE SAME WAY.\n');
   fprintf(fid,'\n');
end

fprintf(fid,'else\n');
fprintf(fid,'   error(''%s',setupFile);
fprintf(fid,': Illegal problem number'');\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

switch pCode
  case 'qp'
    %QP
    fprintf(fid,'\n');

    fprintf(fid,'Prob=tomFiles(Prob,''qp_f'',''qp_g'',''qp_H'');\n');
    fprintf(fid,'\n');

    fprintf(fid,'Prob=qpProbSet(Prob, Name, P, ...\n');
    fprintf(fid,'   x_0, x_L, x_U, x_min, x_max, f_Low, ');
    fprintf(fid,'xName, x_opt, f_opt,...\n');
    fprintf(fid,'   cName, A, b_L, b_U, F, c, B, pSepFunc, uP, uPName);\n');
  case 'mip'
    %MIP
    fprintf(fid,'\n');

    fprintf(fid,'Prob=tomFiles(Prob,''lp_f'',''lp_g'',''lp_H'');\n');
    fprintf(fid,'\n');

    fprintf(fid,'Prob=mipProbSet(Prob, Name, P, ...\n');
    fprintf(fid,'   x_0, x_L, x_U, x_min, x_max, f_Low, ');
    fprintf(fid,'xName, x_opt, f_opt,...\n');
    fprintf(fid,'   cName, A, b_L, b_U, c, B, ...\n');
    fprintf(fid,'   xpcontrol, KNAPSACK, VarWeight, IntVars, ...\n');
    fprintf(fid,'   pSepFunc, uP, uPName);\n');
  case 'lp'
    fprintf(fid,'\n');

    fprintf(fid,'Prob=tomFiles(Prob,''lp_f'',''lp_g'',''lp_H'');\n');
    fprintf(fid,'\n');

    fprintf(fid,'Prob=qpProbSet(Prob, Name, P, ...\n');
    fprintf(fid,'   x_0, x_L, x_U, x_min, x_max, f_Low, ');
    fprintf(fid,'xName, x_opt, f_opt,...\n');
    fprintf(fid,'   cName, A, b_L, b_U, F, c, B, pSepFunc, uP, uPName);\n');
end


fprintf(fid,'\n');
fclose(fid);

% MODIFICATION LOG
%
% 041117  med  xxx_prob removed and code added