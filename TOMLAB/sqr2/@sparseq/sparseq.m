function H = sparseQ(H)
%sparseQ Sparse representation class of Q from sqr
%   990127 Version 1.1
%   Mikael Adlers, Linköpings University
%   e-mail: miadl@mai.liu.se
%
%   This overloads the sparse matrix datatype for sparse Q from sqr
%   This class overload the following commands
%   ', * 
%   double (convert to dense matrix)
%   nnz
%   size
%   sparse (convert to sparse matrix)
%   spy
%   subsref

%H.front represets the orthogonal manipulations on A
%H.front(i).H is the householder vectors for the QR on the i:th frontal matrix
%H.front(i).p is the rows in the i:th frontal matrix
%H.Pr permutates A into columnleading order
%H.rowperm is the final row permutation of A to R

superiorto('double');
superiorto('sparse');

if nargin==0,
  H.size=0;
  H.rowperm=[];
  H.front=[];
  H.nelim=[];
  H.nsteps=[];
  H.nstk=[];
  H.nle=[];
  H.storage='';
  H.transpose='No';
  H=class(H,'sparseq');
elseif nargin ==1,
  if ~isa(H,'sparseq')
    H=class(H,'sparseq');
  end
end

