function [nsteps,nelim,nstk,nle,Pc,Pr,nnzR] = sqr2A(A,nemin)
% SQRA  Analysis routine for SQR.
%       @(#)sqrA.m Version 1.12 3/21/97
%       Pontus Matstoms, University of Linkoping.
%       e-mail: pomat@math.liu.se
%
%       [nsteps,nelim,nstk,nle,Pc,Pr] = sqrA(A)  
%       [nsteps,nelim,nstk,nle,Pc,Pr] = sqrA(A,nemin)
%
%       sqrA performes the symbolic part of the factorization. It computes
%          nsteps      Number of elimination steps.
%          nelim       Number of columns to eliminate in each step.
%          nstk        Number of sons to each elimination tree node.
%          nle         Number of rows with leading entries in each column.
%          Pc          Column ordering making the tree postordered.
%          Pr          Leading entry order of APc.

[m,n]=size(A);

if nargin == 1, nemin=10; end              % Set the default value of nemin

% --- Compute the postordering.

[eparent,Pc]=sparsfun('coletree',A);       % Elimination tree, Pc postordering.

Pcinv(Pc)=1:n;                             % Here we compute
Pcinvmf=[0 Pcinv];                         % the parent vector
eparent(Pcinv)=Pcinvmf(eparent+ones(1,n)); % for the postordered matrix.

[j,i]=find([A(:,Pc)' ; ones(1,size(A,1))]);
[le,Pr]=sort(j(find(diff([0 i']))));      % Leading entry pos for each row
num_rows=full(sparse(1,le,1));            % n._r.(i) number of rows le=i
num_rows=[num_rows zeros(1,n-length(num_rows))];
idx=find(le <= size(A,2));
Pr=Pr(idx);
le=le(idx);

% --- Node amalgamation.
[nodes,parent,nnzR]=amalg(A(:,Pc),eparent,nemin);

% --- Extract the interesting information from the parent vector 'parent'.

nsteps=length(parent);                    % # of el.steps
nelim=diff(nodes);                        % # of col. to el. in each step
nstk=zeros(1,nsteps); 
dads=find(parent);
ndads=full(sparse(1,parent(dads),1));
if ndads == 0,                            % Only disconnected nodes 
   ndads=[]; 
else
   nstk=ndads;
end

% --- Compute a row ordering making A ordered by the leading entries.
%     Rows identically zero are removed by the choice of Pr.

rows=find(diff([le' inf]));            
nle=zeros(1,n);

nle(le(rows))=diff([0 rows]);             % # of rows with LE=i

nle=[0 cumsum(nle)];
a=cumsum([1 nelim(1:(nsteps-1))])';
b=cumsum(nelim)+1;
nle=nle(b)-nle(a);
