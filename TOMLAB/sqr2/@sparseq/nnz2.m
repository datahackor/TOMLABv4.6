function nz = nnz2(H)
%sparseQ/nnz     Number of nonzero elements in sparseQ
%        Overloads nnz(A)
%        %Z%%M% Version %I% %G%
%        Mikael Adlers, Linköpings University
%        e-mail: miadl@mai.liu.se

nz=0;
for i=1:length(H.front)
  nz=nz+nnz(H.front(i).H);
  nz=nz+length(H.front(i).p);
  if (H.storage=='H')
    nz=nz+length(H.front(i).tau);
  end
end
nz=nz+length(H.rowperm);
