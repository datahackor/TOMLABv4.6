function C=appH(Q,B)
% appH  Applies the orthogonal matrix Q on a vector
%       2000 Version 1.2
%       Mikael Adlers, University of Linkoping.
%       e-mail: miadl@mai.liu.se
%
%       C=appH(H,B) computes the product Q*B where Q is represented by the
%       householder vectors H from the sparse multifrontal qr routine sqr

% HH.nelim is the number of rows taken from the frontal matrix in step i
% HH.H represets the orthogonal manipulations on A
% HH.H(i).frontH is the householder vectors for the QR on the i:th frontal matrix
% HH.H(i).p is the rows in the i:th frontal matrix
% HH.Pr permutates A into columnleading order
% HH.rowperm is the final row permutation of A to R

% This is what happens in sqr: rowperm*(Qn'*...*Q1')Pr A = R

% Version history
% 1.1 Added handling of Householder storage
% 1.2 Added removal of machine precision fill-in

if issparse(B),
  maxB=max(max(B));
end

B(Q.rowperm,:)=B;
for i=length(Q.front):-1:1
  if (Q.storage=='Q'),
    B(Q.front(i).p,:)=Q.front(i).H * B(Q.front(i).p,:);
  else
     %b=appL(double(B(Q.front(i).p,:)),double(Q.front(i).H),Q.front(i).tau,1);
     if (issparse(B)),
       [r,c,val]=find(B(Q.front(i).p,:));
       [c1,z]=indexA(c);                 % Convert to local coordinates
       S = double(sparse(r,c1,val)); % Compact sparse matrix
     else
       S=B(Q.front(i).p,:);
     end
     b=appL(S,Q.front(i).H,Q.front(i).tau,1);
%     f=find(abs(b)< 100*eps*max(abs(b)));
%     b(f)=0;
     if issparse(B), % Remove machine precision fill-in
       b(abs(b) < 10*eps*maxB)=0;
       [r,c1,val]=find(sparse(b));
       B(Q.front(i).p,:)=sparse(r,z(c1),val,size(b,1),size(B,2));
%       B(Q.front(i).p,c(find(diff([0; c]))))=b;
     else
       B(Q.front(i).p,:)=b;    
     end

%     B(Q.front(i).p,:)=b;
  end
end
%if issparse(B),
%  C=sparse([],[],[],size(B,1),size(B,2),nnz(B));
%else
%  C=zeros(size(B));
%end
%C(Q.Pr,:)=B;
C=B;

function [A,zm]=indexA(A);
[m,n]=size(A);
z=zeros(max(A),1);
z(A)=A;
zm=find(z);
z(zm)=1:length(zm);
A=z(A);

