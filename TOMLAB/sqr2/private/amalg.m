function [nodes,parent,nnzR]=amalg(A,eparent,nemin)
% AMALG	Node amalgamation for sqr. 
%       @(#)amalg.m Version 1.13 3/21/97
%       Pontus Matstoms, Linkoping University.
%       e-mail: pomat@math.liu.se
%       
%       [nodes,parent]=amalg(A,eparent,nemin) computes the elimination tree
%       where nodes of 'eparent' are amalgamated into supernodes.
%
%       Method: nemin amalgamation (See Matstoms, "The multifrontal 
%               solution ...") with repeated nemin step.
%
%       nemin=0     No amalgamation.
%       nemin=1     Fundamental supernodes formed.
%       nemin>1     Following the algorithm.
%
%       The computation of fundamental supernodes is due to John Gilbert.

n=size(A,2);

count=symbfact(A,'col');                % # of elements in the rows of R.
nnzR=sum(count);                        % # of elements in T

if nemin == 0 | n<=1,
   nodes=1:(n+1);
   parent=eparent;
   disp(['Return from amalg without node amalgamation'])
   return
end

kids=find(eparent);                     % Nodes with a father node.
nkids=zeros(1,n);
nkids=full(sparse(1,eparent(kids),1));  % # of sons to each father node.

% Form fundamental supernodes.

% The vector IsNode is now set to the nodes of the amalgamated 
% elimination tree. If a supernode contains the single nodes i,...,i+t,
% only node i is marked as a node in IsNode. Elements corresponding to 
% nodes are set to one and the other elements to zero.

IsNode  = zeros(1,n);                   
f = find(nkids ~= 1);                   % Nodes without or more than one sons.
%IsNode(f) = ones(size(f));              
IsNode(f) = 1;              

% Father nodes whose corresponding rows in R not are included in
% the structure of the sons.

if ~isempty(kids),
%  f = find ( count(eparent(kids)) ~= count(kids)-ones(size(kids)) );
  f = find ( count(eparent(kids)) ~= count(kids) - 1 );
else
  f=[];
end
DadsANode = kids(f);
%IsNode(eparent(DadsANode)) = ones(size(DadsANode));
IsNode(eparent(DadsANode)) = 1;

% Supernodes now formed.

% First phase of nemin amalgamation.

GrandSon=1:n;
AmalgNodes=find(IsNode == 0);           % Amalgamated nodes
for i=AmalgNodes,
   GrandSon(i)=GrandSon(i-1);           % GrandSon(i) gives the youngest member
end                                     % in a node comtaining node i.

nodes = sort(find(IsNode));             % Real nodes in the partially amalgamated 
Nnodes = length(nodes);                 % tree.

TopAmalgNodes=AmalgNodes(find(diff([AmalgNodes 0]) ~= 1));  % eldest sons in supern.
BotAmalgNodes=nodes(find(diff([nodes n+1]) ~= 1));          % youngest sons in supern.

eparent(BotAmalgNodes)=eparent(TopAmalgNodes);
idx=find(eparent);
eparent(idx)=GrandSon(eparent(idx));

dads=sort(eparent(nodes(1:(Nnodes-1))));
dads=dads(find(diff([0 dads])));        % Father nodes in the amalgamated tree.

SNnumber=cumsum(IsNode);                
idx=find(diff([0 SNnumber]));
FSV(idx)=[1 diff(find([diff(SNnumber) 1]))];
FSV(TopAmalgNodes)=FSV(BotAmalgNodes);  % Number of fully summed variables in nodes.

GrandDad=1:nodes(Nnodes);
GrandDad(BotAmalgNodes)=TopAmalgNodes;

for DadsNode=dads,
   if max(FSV(DadsNode),FSV(DadsNode-1)) < nemin, 
      FSV(GrandDad(DadsNode))=FSV(DadsNode)+FSV(DadsNode-1);
      IsNode(DadsNode)=0;
      eparent(DadsNode-1)=eparent(GrandDad(DadsNode));
   end
end

% Second phase of nemin amalgamation.

GrandSon=1:n;
idx=find(IsNode == 0);
for i=idx,
   GrandSon(i)=GrandSon(i-1);
end

nodes = sort(find(IsNode));   
Nnodes = length(nodes);

AmalgNodes=find(IsNode == 0);
TopAmalgNodes=AmalgNodes(find(diff([AmalgNodes 0]) ~= 1));
BotAmalgNodes=nodes(find(diff([nodes n+1]) ~= 1));

eparent(BotAmalgNodes)=eparent(TopAmalgNodes);
idx=find(eparent);
eparent(idx)=GrandSon(eparent(idx));

dads=sort(eparent(nodes(1:(Nnodes-1))));
dads=dads(find(diff([0 dads]))); 

SNnumber=cumsum(IsNode);         
idx=find(diff([0 SNnumber]));
FSV(idx)=[1 diff(find([diff(SNnumber) 1]))];
FSV(TopAmalgNodes)=FSV(BotAmalgNodes);

GrandDad=1:nodes(Nnodes);
GrandDad(BotAmalgNodes)=TopAmalgNodes;

for DadsNode=dads,
   if max(FSV(DadsNode),FSV(DadsNode-1)) < nemin, 
      FSV(GrandDad(DadsNode))=FSV(DadsNode)+FSV(DadsNode-1);
      IsNode(DadsNode)=0;
      eparent(DadsNode-1)=eparent(GrandDad(DadsNode));
   end
end

% Compute the parent-vector for the amalgamated tree.

SNnumber = cumsum(IsNode);
nodes = [sort(find(IsNode)) n+1];
nnodes = length(nodes)-1;

LastVtx = nodes-1;
LastVtx = LastVtx(2:(nnodes+1));
f = find( eparent(LastVtx) );
LastVtxOfKidNodes = LastVtx(f);
parent = zeros(1,nnodes);
parent(SNnumber(LastVtxOfKidNodes)) = SNnumber(eparent(LastVtxOfKidNodes));
