function C=appHT(Q,B)
% appHT Applies the orthogonal matrix Q^T on a vector
%       1999 Version 1.1
%       Mikael Adlers, University of Linkoping.
%       e-mail: miadl@mai.liu.se
%
%       C=appHT(H,B) computes the product Q^T*B where Q is represented by the
%       householder vectors H from the sparse multifrontal qr routine sqr

% HH.H represets the orthogonal manipulations on A
% HH.H(i).frontH is the householder vectors for the QR on the i:th frontal matrix
% HH.H(i).p is the rows in the i:th frontal matrix
% HH.Pr permutates A into columnleading order
% HH.rowperm is the final row permutation of A to R

% This is what happens in sqrQ: rowperm*(Qn'*...*Q1')Pr A Pc = R

if issparse(B)
  maxB=max(max(B));
  B=B'; %Using B row wise, matlab stores columnvise
end
%B=B(Q.Pr,:);
for i=1:length(Q.front)
  if (Q.storage=='Q'),
    b=Q.front(i).H'*B(:,Q.front(i).p)';
  else
    if (issparse(B)),
      [c,r,val]=find(B(:,Q.front(i).p));
      [c1,z]=indexA(c);                 % Convert to local coordinates
      S = double(sparse(r,c1,val)); % Compact sparse matrix
    else
      S=B(Q.front(i).p,:);
    end
    b=appL(S,Q.front(i).H,Q.front(i).tau,2); 
%    b=appL(double(B(Q.front(i).p,:)),Q.front(i).H,Q.front(i).tau,2); 
  end
  if issparse(B), % Remove machine precision fill-in
    f=find(abs(b) < 10*eps*max(max(abs(b))));
    b(abs(b) < 10*eps*max(max(abs(b))))=0;
%    B(Q.front(i).p,c(find(diff([0; c]))))=b;
    [r,c1,val]=find(sparse(b));
    B(:,Q.front(i).p)=sparse(z(c1),r,val,size(B,1),size(b,1));;
   else
    B(Q.front(i).p,:)=b;    
  end
end
if issparse(B),
   C=B(:,Q.rowperm)';
else
  C=B(Q.rowperm,:);
end
function [A,zm]=indexA(A);
% This function relaces indices with the local index
% [1 4 4 1 5] -> [1 2 2 1 3]
% zm contains local to global transformation ([1 4 5] in the example above)

[m,n]=size(A);
z=zeros(max(A),1);
z(A)=A;
zm=find(z);
z(zm)=1:length(zm);
A=z(A);
%A=reshape(z(A),m,n);



