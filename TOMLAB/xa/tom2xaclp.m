% function xaclp = tom2xaclp(xaControl,xaclp)
%
% Function to convert the allowed TOMLAB /XA parameters to XACLP ('command
% line parameter') strings. 
%
% Can append commands to the optional input xaclp (cell array)
%
% xaControl is a structure array like this:
%
% xaControl.PRICING = 1;
% xaControl.EUROPEAN = 'Yes';    % Can also be given numerically - yes=1,no=0
% ... and so on, as desired.
%
% It is important that the parameter identifiers, such as PRICING, is in
% upper case letters. Otherwise, they will not be recognized.
%
% The resulting output would be a cell array:
% 
% xaclp={ 'Set Pricing 1', 'European Yes' }
%
%
% Notes: 
% Currently missing parameters (at least):
%                                          (Filename)
%                                          (Memory)
%                                          (Debug)
%                                          (Torcc)
%
% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomlab.biz
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 14.0 $
% Written Oct 16, 2003.      Last modified Jan 17, 2005.

function xaclp = tom2xaclp(xaControl,xaclp)

if nargin < 2
   % Cell equivalent
   xaclp={};
   
   if nargin < 1
      error('tom2xaclp requires at least one input argument, the structure array xaControl');
   end
end


if ~isempty(xaControl)
   if ~isstruct(xaControl)
      error('xaControl must be a structure');
   end
end

if ~iscell(xaclp)
   error('tom2xaclp: second argument xaclp is not a cell array');
end

if isempty(xaControl)
  return
end


fs = fieldnames(xaControl);

for k=1:size(fs,1)
   fn = deblank(fs(k));
   fi = getfield(xaControl,c2m(fn));
   
   % lpmethod
   if(strcmp(fn,'LPMETHOD'))
       if(isnumeric(fi))
           switch fi
            case 0
                % Automatic (XA 13 = primal, XA 14 = dual)
            case 1
                % Primal simplex
                xaclp{end+1} = 'Set DualSimplex No';
                xaclp{end+1} = 'Set Barrier No';
            case 2
                % Dual simplex
                xaclp{end+1} = 'Set DualSimplex Yes';
                xaclp{end+1} = 'Set Barrier No';
            case 3
                % Barrier
                xaclp{end+1} = 'Set Barrier Yes';
            case 4
                % Network primal
               xaclpwarn([fn{1} ' parameter illegal numeric value ' num2str(fi) '. Network solver not available. Ignored' ]);
            case 5
                % Network dual
               xaclpwarn([fn{1} ' parameter illegal numeric value ' num2str(fi) '. Network solver not available. Ignored' ]);
            case 6
                % Network generalized primal
               xaclpwarn([fn{1} ' parameter illegal numeric value ' num2str(fi) '. Network solver not available. Ignored' ]);
            case 7
                % Network generalized dual
               xaclpwarn([fn{1} ' parameter illegal numeric value ' num2str(fi) '. Network solver not available. Ignored' ]);
            otherwise
               xaclpwarn([fn{1} ' parameter illegal numeric value ' num2str(fi) ' ignored' ]);
           end
       else
           xaclpwarn( wnn(fn) );
       end

   % Strategy, char or numeric
   elseif(strcmp(fn,'STRATEGY'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Strategy ' num2str(floor(fi(1))) ];
      elseif(ischar(fi))
         xaclp{end+1} = [ 'Strategy ' fi ];
      else
         xaclpwarn( ['Ignoring xaControl.' char(fn) '. Invalid format.'] );
      end

   % Presolve, char
   elseif(strcmp(fn,'PRESOLVE'))
      if(ischar(fi))
         xaclp{end+1} = [ 'Presolve ' fi ];
      else
         xaclpwarn(wnc(fn));
      end         
       
   % Eliminate, yes/no/off
   elseif(strcmp(fn,'ELIMINATE'))
      if(isyes(fi))
        xaclp{end+1} = 'Set Eliminate Yes';
      elseif(isno(fi))
        xaclp{end+1} = 'Set Eliminate No';
      elseif(isoff(fi))
        xaclp{end+1} = 'Set Eliminate Off';
      else
        xaclpwarn([fn{1} ' parameter illegal value ' num2str(fi) ' ignored' ]);
      end

   % Basis, char
   elseif(strcmp(fn,'BASIS'))
      if(ischar(fi))
         xaclp{end+1} = [ 'Basis ' fi ];
      else
         xaclpwarn(wnc(fn));
      end         

   % Maximize, yes/no
   elseif(strcmp(fn,'MAXIMIZE'))
      if(isyes(fi))
         xaclp{end+1} = 'Maximize Yes';
      else
         xaclp{end+1} = 'Maximize No';
      end
       
   % Crash, integer 0,1,2,3
   elseif(strcmp(fn,'CRASH'))
      if(isnumeric(fi))
         xaclp{end+1} = ['Set Crash ' num2str(fi(1)) ];
      else
         xaclpwarn( wnn(fn) );
      end
   
   % Decimals, integer
   elseif(strcmp(fn,'DECIMALS'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Decimals ' num2str(floor(fi(1))) ];
      else
         xaclpwarn( wnn(fn) );
      end
   
   % European, yes/no
   elseif(strcmp(fn,'EUROPEAN'))
      if(isyes(fi))
         xaclp{end+1} = 'European Yes';
      else
         xaclp{end+1} = 'European No';
      end
   
   % Fieldsize, integer
   elseif(strcmp(fn,'FIELDSIZE'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Fieldsize ' num2str(floor(fi(1))) ];
      else
         xaclpwarn( wnn(fn) );
      end
   
   % Limitsearch, char
   elseif(strcmp(fn,'LIMITSEARCH'))
      if(ischar(fi))
         xaclp{end+1} = [ 'Limitsearch ' fi ];
      end
   
   % Mute, yes/no
   elseif(strcmp(fn,'MUTE'))
      if(isyes(fi))
         xaclp{end+1} = 'Mute Yes'; 
      else
         xaclp{end+1} = 'Mute No';
      end
   
   % Degeniter, integer
   elseif(strcmp(fn,'DEGENITER'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Degeniter ' num2str(floor(fi(1))) ];
      else
         xaclpwarn( wnn(fn) );
      end
   
   % Elemsize, real
   elseif(strcmp(fn,'ELEMSIZE'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set ElemSize ' num2str(fi(1)) ];
         %sprintf('Set Elemsize %24.14g',double(fi));
      else
         xaclpwarn(wnn(fn));
      end
   
   % Freqlog, char
   elseif(strcmp(fn,'FREQLOG'))
      if(ischar(fi))
         xaclp{end+1} = [ 'Set FreqLog ' fi ];
      else
         xaclpwarn(wnc(fn));
      end         
   
   % Ibounds, real
   elseif(strcmp(fn,'IBOUNDS'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set IBounds ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Intgap, real
   elseif(strcmp(fn,'INTGAP'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Intgap ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Intlimit, integer
   elseif(strcmp(fn,'INTLIMIT'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Intlimit ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Intpct, real 
   elseif(strcmp(fn,'INTPCT'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Intpct ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end
      
   % IRound, Yes/No
   elseif(strcmp(fn,'IROUND'))
      if(isyes(fi))
         xaclp{end+1} = 'Set IRound Yes';
      else
         xaclp{end+1} = 'Set IRound No';
      end
   
   % Iteration, integer
   elseif(strcmp(fn,'ITERATION'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Iteration ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Limitnodes, integer
   elseif(strcmp(fn,'LIMITNODES'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set LimitNodes ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % LTolerance, real 
   elseif(strcmp(fn,'LTOLERANCE'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set LTolerance ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end
      
   % UTolerance, real 
   elseif(strcmp(fn,'UTOLERANCE'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set UTolerance ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end

   % Markowitz, real 
   elseif(strcmp(fn,'MARKOWITZ'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Markowitz ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end

   % Maxcpu, integer (unsure how this works with Matlab ???)
   elseif(strcmp(fn,'MAXCPU'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set MaxCPU ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Maxnodes, integer
   elseif(strcmp(fn,'MAXNODES'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set MaxNodes ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % MPricing, integer
   elseif(strcmp(fn,'MPRICING'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set MPricing ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Pertubate, real 
   elseif(strcmp(fn,'PERTUBATE'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Pertubate ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end
  
   % Pricing, integer 0,1,2,3,4
   elseif(strcmp(fn,'PRICING'))
      if(isnumeric(fi))
         if(fi(1) >= 0 & fi(1) <= 4)
            xaclp{end+1} = [ 'Set Pricing ' num2str(round(fi(1))) ];
         else
            xaclpwarn( ['Ignoring ' char(fn) ' of bounds value. Allowed: 0,1,2,3 or 4' ]);
         end
      else
         xaclpwarn(wnn(fn));
      end
   
   
   % Tolerance Dual, real 
   elseif(strcmp(fn,'TOLERANCE_DUAL'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Tolerance Dual ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end

   
   % Reinvertfreq, integer
   elseif(strcmp(fn,'REINVERTFREQ'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set ReInvertFreq ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Rejpivot, integer
   elseif(strcmp(fn,'REJPIVOT'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set RejPivot ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end

   
   % Relaxed, yes/no
   elseif(strcmp(fn,'RELAXED'))
      if(isyes(fi))
         xaclp{end+1} = 'Set Relaxed Yes';
      else
         xaclp{end+1} = 'Set Relaxed No';
      end

   
   % Scale, yes,no,2... custom type:
   elseif(strcmp(fn,'SCALE'))
      if(isnumeric(fi))
         switch(fi)
            case 0
               xaclp{end+1} = 'Set Scale No';
            case 1
               xaclp{end+1} = 'Set Scale Yes';
            case 2
               xaclp{end+1} = 'Set Scale 2';
            otherwise
               xaclpwarn([fn{1} ' parameter illegal numeric value ' num2str(fi(1)) ' ignored' ]);
         end
      elseif(ischar(fi) & ~isempty(fi))
         if(lower(fi(1))=='n')
            xaclp{end+1} = 'Set Scale No';
         elseif(lower(fi(1))=='y')
            xaclp{end+1} = 'Set Scale Yes';
         elseif(fi(1)=='2')
            xaclp{end+1} = 'Set Scale 2';
         end
      end
   
   % Sprouts, integer
   elseif(strcmp(fn,'SPROUTS'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Sprouts ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end

   % Runner, yes/no
   elseif(strcmp(fn,'RUNNER'))
      if(isyes(fi))
         xaclp{end+1} = 'Set Runner Yes';
      else
         xaclp{end+1} = 'Set Runner No';
      end
   
   % Stickwithit, integer
   elseif(strcmp(fn,'STICKWITHIT'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set StickWithIt ' num2str(round(fi(1))) ];
      else
         xaclpwarn(wnn(fn));
      end

   % Tolerance_tcoefficients, real
   elseif(strcmp(fn,'TOLERANCE_TCOEFFICIENTS') | strcmp(fn,'tol_tcoef'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Tolerance TCoefficients ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end

   % Timelimit, integer
   elseif(strcmp(fn,'TIMELIMIT'))
      if(isnumeric(fi))
         if(fi == inf)
            fi = 59+59*60+999*60^2;
         end
         max_s = mod(fi,60);
         fi = floor(fi/60);
         max_m = mod(fi,60);
         fi = floor(fi/60);
         max_h = fi;
         xaclp{end+1} = [ 'Set Timelimit ' sprintf('%02i:%02i:%02i', max_h, max_m, max_s) ];
      else
         xaclpwarn(wnn(fn));
      end

   % Transentry, real
   elseif(strcmp(fn,'TRANSENTRY'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set TransEntry ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end
   
   % Tolerance primal, tolerance_primal or tol_prim, real
   elseif(strcmp(fn,'TOLERANCE_PRIMAL') | strcmp(fn,'TOL_PRIM'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set Tolerance Primal ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end

   % YPivot, real
   elseif(strcmp(fn,'YPIVOT'))
      if(isnumeric(fi))
         xaclp{end+1} = [ 'Set YPivot ' num2str(fi(1)) ];
      else
         xaclpwarn(wnn(fn));
      end

   % StopAfter, char
   elseif(strcmp(fn,'STOPAFTER'))
      if(ischar(fi))
         xaclp{end+1} = [ 'Set StopAfter ' fi ];
      else
         xaclpwarn(wnc(fn));
      end         

   % StopUnchanged, char
   elseif(strcmp(fn,'STOPUNCHANGED'))
      if(ischar(fi))
         xaclp{end+1} = [ 'Set StopUnchanged ' fi ];
      else
         xaclpwarn(wnc(fn));
      end         
   
   % Warning, yes/no
   elseif(strcmp(fn,'WARNING'))
      if(isyes(fi))
         xaclp{end+1} = 'Set Warning Yes';
      else
         xaclp{end+1} = 'Set Warning No';
      end
  else
      xaclpwarn(['Invalid command xaControl.' fn{1} ' ignored.']);
  end
end

return


% ----------------------------------------------------
% Our own warning routine

function xaclpwarn(s)
fprintf('--- tom2xaclp: Warning: %s\n',s);

% ----------------------------------------------------
function s = wnn(fn)
s = ['Ignoring non-numeric data for xaControl.' char(fn) ];

% ----------------------------------------------------
function s = wnc(fn)
s = ['Ignoring non-char data for xaControl.' char(fn) ];


% ----------------------------------------------------
function y=isyes(s) % 1

if isempty(s)
   y=0;
elseif ischar(s)
   y=double(s(1)=='y' | s(1)=='Y');
elseif isnumeric(s) | islogical(s)
   y=double(s(1)==1);
else
   y=0;
end

function n=isno(s) % 0

if isempty(s)
   n=0;
elseif ischar(s)
   n=double(s(1)=='n' | s(1)=='N');
elseif isnumeric(s) | islogical(s)
   n=double(s(1)==0);
else
   n=0;
end

function o=isoff(s) % 2

if isempty(s)
   o=0;
elseif ischar(s)
   o=double(strcmpi(s, 'off'));
elseif isnumeric(s)
   o=double(s(1)==2);
else
   o=0;
end


function [m] = c2m(c)
% c2m Combine a cell array of matrices into one matrix.
elements = prod(size(c));

if elements == 1
  m = c{1};
  return
end

if elements == 0
  m = [];
  return
end

[rows,cols] = size(c);

if (cols == 1)
  m = cell(1,rows);
  for i=1:rows
    m{i} = [c{i}]';
  end
  m = [m{:}]';
  
else
  m = cell(1,rows);
  for i=1:rows
    m{i} = [c{i,[1:cols]}]';
  end
  m = [m{1,[1:rows]}]';
end

% MODIFICATION LOG
%
% 041013 ango Case insensitive in YES detection
% 041021 frhe Added parameter lpmethod and fixed a minor bug.
% 041103 frhe Removed Barrier-example from help
% 041103 frhe All parameter identifiers are now uppercase and case 
%             sensitive. This is to avoid ambiguous parameters, 
%             and to set the default parameters correctly (in xaTL.m).
% 050117 med  mlint revision
