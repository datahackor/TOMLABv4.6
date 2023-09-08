%
% function bosopt = tom2bosopt(options,validlist,bosopt)
%
% Function to convert the allowed TOMLAB /SPRNLP /BARNLP parameters
% to BOEING OPTION ('symbol=value') strings. 
%
% validlist is a cell arry with valid options. Options in options
% but not present in validlist will be ignored.
%
% Can append commands to the optional input bosopt (cell array)
%
% options is a structure array like this:
%
% options.NITMAX = 200;
% options.MAXNFE = 10000;
% ... and so on, as desired.
%
% It is important that the parameter identifiers are in
% upper case letters. Otherwise, they will not be recognized.
%
% The resulting output would be a cell array:
% 
% bosopt={ 'NITMAX=200', 'MAXNFE=10000' }
%

function bosopt = tom2bosopt(options,validlist,bosopt)

if nargin < 3
   % Cell equivalent
   bosopt={};
   
   if nargin < 2
      error(['tom2bosopt requires at least two input arguments, the ' ...
             'structure options and the cell array validlist']);
   end
end


if ~isempty(options)
   if ~isstruct(options)
      error('tom2bosopt: options must be a structure');
   end
end

if ~iscell(validlist)
   error('tom2bosopt: second argument validlist is not a cell array');
end

if ~iscell(bosopt)
   error('tom2bosopt: third argument bosopt is not a cell array');
end

if isempty(options)
  return
end

ints = 'ijklmnIJKLMN';

fs = fieldnames(options);

for k=1:size(fs,1)
   fn = deblank(fs(k));
   fi = getfield(options,c2m(fn));

   if(isempty(fi))
     continue;
   end
   
   % Check if the options in in the list of valid options   
   if(~isempty(find(strcmp(fn, validlist))))
     % Create option string
     if(isnumeric(fi))
       fnname = c2m(fn);
       if(isempty(strfind(ints, fnname(1))))
         bosopt{end+1} = [c2m(fn) '=' num2str(fi,'%20.15f')];
       else 
         bosopt{end+1} = [c2m(fn) '=' int2str(fi)];
       end
     elseif(ischar(fi))
       bosopt{end+1} = [c2m(fn) '=' fi];
     else
       bosoptwarn(['Ignoring non-numeric or non-char option: ' ...
                  char(fn)]);
     end
   else
     bosoptwarn(['Ignoring invalid option: ' char(fn)]);
   end
end

return


% ----------------------------------------------------
% Our own warning routine

function bosoptwarn(s)
fprintf('--- tom2bosopt: Warning: %s\n',s);

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
% 041013 frhe File created, based on tom2xaclp
