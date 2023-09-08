function [R,Pc,C] = sqr2(A,B)
%SQR	Sparse orthogonal-triangular decomposition.
%       %Z%%M% Version %I% %G%
%       Copyright by Mikael Adlers.
%
%       Pontus Matstoms, Linkoping University.
%       Mikael Adlers, Linköpings University
%       e-mail: pomat@mai.liu.se, miadl@mai.liu.se
%
%       [R,p]=sqr(A) computes the upper triangular matrix R in the
%       QR factorization,
%                          Q'(AP)=[R 0]'.
%       The permutation matrix P is represented by the integer vector p.
%
%       [R,p,Q]=sqr(A) computes also an implicit representation of the 
%       orthogonal matrix Q. Q is of datatype sparseQ and can be used almost
%       in tha same manor as a sparse matrix. To convert Q to the sparse
%       format use sparse(Q) or double(Q) to a full matrix, however not 
%       recomended. 
%
%       [R,p,C]=sqr(A,B) does also compute the n first components of C=Q'B,
%       where n is the column dimension of A. The matrix is here ordered
%       by rows in the same way as the multifrontal scheme orders the rows
%       of A.
%
%       If the GLOBAL variable sqr_nemin exists, it is chosen as amalgamation
%       parameter in sqrA.

% Call the analysis routine.
[m,n]=size(A);
global sqr_nemin
if isempty(sqr_nemin),
   [nsteps,nelim,nstk,nle,Pc,Pr,nnzR] = sqr2A(A);
else
   disp(['nemin is now ',num2str(sqr_nemin)])
   [nsteps,nelim,nstk,nle,Pc,Pr,nnzR] = sqr2A(A,sqr_nemin);
end
C.size=max(m,n);
% Order A by rows (Pr) and columns (Pc).

A=A(Pr,Pc);
[m,n]=size(A);

% Call the factorization routine.

if nargin == 1,            % Only factorization
  if nargout < 3,          % Compute only R
    R = sqr2B(A,nsteps,nelim,nstk,nle,Pr,nnzR);
  elseif nargout ==3,      % Compute representation of Q
    [R,front,rowperm] = sqr2B(A,nsteps,nelim,nstk,nle,Pr,nnzR);
    C.rowperm=Pr([rowperm;((length(rowperm)+1):m)']); 
    C.front=front;
    C.nelim=nelim;
    C.nsteps=nsteps;
    C.nstk=nstk;
    C.nle=nle;
    if isfield(front,'tau'),
      C.storage='H';
    else
      C.storage='Q';
    end
    C.transpose='No';
    C=sparseq(C); % Convert to sparseq format
  else
    error('Wrong number of output parameters, see help info');    
  end
else                       % Righthand side
  if nargout <=3,           % Only righthand side
     [R,C] = sqr2B(A,nsteps,nelim,nstk,nle,Pr,nnzR,B(Pr,:));
  else
    error('Wrong number of output parameters, see help info');    
  end
end






