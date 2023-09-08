function [R,H,rowperm] = sqr2B(A,nsteps,nelim,nstk,nle,Pr,nnzR,B)
% SQRB  Factorization routine in SQR.
%       %Z%%M% Version %I% %G%
%       Pontus Matstoms, Linkoping University.
%       Mikael Adlers, Linköping University.
%       e-mail: pomat@mai.liu.se, miald@mai.liu.se
%
%       [R,C] = sqrB(A,nsteps,nelim,nstk,nle) computes the upper triangular 
%       factor R in the QR factorization of the sparse m-by-n matrix A. If a 
%       matrix B is passed to the routine, [R,C] = sqrB(A,nsteps,nelim,nstk,nle,B),
%       then the first n components of Q'B are also computed.
%

% nsteps         The number of nodes/supernodes in the elimination tree.
% nelim()        The number of rows from A to include in ith frontal matrix.
% nstk()         The number of stacked elements to merge in the current step.
% R(,)           The upper triangular matrix R.
% STK()          The stack of update matrices.
% STKcols        Column indices of stacked update matrices.
% STKm           Row dimension of stacked update matrices.
% STKn           Column dimension of stacked update matrices.
% BSTK           Right-hand side stack.

% PSTK           Stack of row indexes. Used when computing Q
% rowidx         Row index

[m,n]=size(A);
% Check the input
computeB=(nargin==8);
computeH=(nargout>=3);

if computeB, 
%  B=B(Pr,:); 
  [Bm,Bn]=size(B); 
  C=zeros(n,Bn);
  %else, 
  %  B=0;  
end

R=spalloc(n,n,ceil(nnzR*1.01));

spnm=0;          % Pointer for STKm and STKn

rp=0;
lmr=0;  

if computeH,
  rowperm=zeros(m,1);
  last=m;
end

% --- Main loop ...

for iass=1:nsteps,
  
  % -- Compute the dimension and the global column indices of the front.
  
  numorg=nelim(iass);
  numstk=nstk(iass);
  
  colflag=zeros(1,n);     % Global column indices in the frontal matrix.
  
  % - Integer data associated with the contribution from A.
  
  Am=nle(iass);   
  FNTm=Am;
  if Am>0,
    [ignore,col]=find(A((lmr+1):(lmr+Am),:)); %Slow
    colflag(col)=ones(size(col));
    %colflag=colidx(A,lmr+1,lmr+Am);
  end
  
  % - Integer data associated with the stack contribution.
 
  if ( numstk > 0 ),
    ncol=0;
    stkm=sum(STKm((spnm-numstk+1):spnm));
    FNTm=FNTm+stkm;
    ncol=sum(STKn((spnm-numstk+1):(spnm)));
    for i=(spnm-numstk+1):spnm,
%      FNTm=FNTm+size(STK(i).A,1);
%      ncol=ncol+size(STK(i).A,2);
%size(STK(i).A),i
%      colflag(STK(i).cols)=ones(size(STK(i).A,2),1);
      colflag(STK(i).cols)=ones(STKn(i),1);
    end
  else
    stkm=0; %unnecessary i think
  end
  
  %    FNTcols=findidx(colflag);
  %    FNTn=length(FNTcols);
  
  FNTcols=find(colflag); %takes a lot of time
  FNTn=size(FNTcols,2);
  colflag(FNTcols)=1:FNTn;          % colflag maps global to local column
  % indices for the frontal matrix.
  % local=colflag(global)
  
  %    FNTn=sum(colflag);
  %    colflag(colflag>0)=1:FNTn;
  
  
  % -- Move reals to the frontal matrix.

  if FNTn > 0,
    FNT=zeros(FNTm,FNTn);
    if computeB,
      BFNT=zeros(FNTm,Bn);      
    end
    
    
    % - From A ...

   if Am>0,
      FNT(1:Am,:)=A((lmr+1):(lmr+Am),FNTcols); 
      Fm=Am;
      
      if computeB,                     % Right-hand side front.
	BFNT(1:Am,:)=B((lmr+1):(lmr+Am),:);
      end
      
      if computeH,
	rowidx=(lmr+1):(lmr+Am); % Rowindex of the rows from A
      end
      
    else
      Fm=0;           %ML 970910
      if computeH,
	rowidx=[];
      end
    end
    
    % - ... and from the stack ...
    
    for s=1:numstk,
            udm=STKm(spnm);
            udn=STKn(spnm);
%      udm=size(STK(spnm).A,1);
      if udm>0,
        col=colflag(STK(spnm).cols);
        FNT((Fm+1):(Fm+udm),col)=STK(spnm).A;
        if computeB,
	  BFNT((Fm+1):(Fm+udm),:)=STK(spnm).b;
        end
        if computeH,
	  rowidx=[rowidx STK(spnm).rows];
        end
	Fm=Fm+udm;
      end
%      STK(spnm)=[]; clear stack memory
      spnm=spnm-1;
    end
  else
    spnm=spnm-numstk;
  end
  numorg=min(numorg,FNTm); % Rank deficient problem may have FNTm < numorg.
  
  % -- Factorization
  
  RFm=min(FNTm,FNTn);
  RFn=FNTn;
  rdim=min(RFm,RFn);

%  global AFNT
%  AFNT=FNT;
  
  if computeH,
    [Q,RR]=qr(FNT);
    RBF=triu(RR);
    H(iass).H=sparse(Q);
    H(iass).p=Pr(rowidx);
%    RR=qr(FNT);
%    RBF=triu(RR);	
%global band
%if (band==0)
%  H(iass).H=sparse(tril(RR,0));
%else
% H(iass).H=sparse(bandqr(FNT));
%end
    if computeB
      Ctemp=Q'*BFNT
      C((rp+1):(rp+numorg),:)=Ctemp(1:numorg);
    end

%subplot(131),spy(FNT),subplot(132),spy(sparse(bandqr(FNT))),subplot(133),spy(sparse(tril(qr(FNT)))),pause    
%subplot(121),spy(sparse(FNT)),subplot(122),spy(sparse(tril(qr(FNT)))),pause    

  else
    if computeB,     % A right-hand side matrix is passed to the routine.
      
      % -    Factorize the frontal matrix and dispose computed elements of Q'B.
      
      RBF=triu(qr([FNT BFNT]));
      C((rp+1):(rp+numorg),:)=RBF(1:numorg,(RFn+1):(RFn+Bn)); 
      
    else
      % -    Factorize the frontal matrix.
      RBF=triu(qr(FNT));	
    end
  end
  
% - Move the first rows of frontals R to the final R.
  
  R(FNTcols,(rp+1):(rp+numorg))=RBF(1:numorg,1:RFn)';
  if computeH,
    rowperm((rp+1):(rp+numorg))=rowidx(1:numorg);
  end
  rp=rp+numorg;
  
% --- Move the remaining rows to the stack.
  
  spnm=spnm+1;
  if rdim > numorg,
    STKm(spnm)=rdim-numorg;
    STKn(spnm)=RFn-numorg;

    STK(spnm).cols=FNTcols((numorg+1):RFn);
    STK(spnm).A=RBF((numorg+1):rdim,(numorg+1):RFn);
%    dispRspy(RBF,STK(spnm).A);
%    dispRspy(R,FNT);

    if computeB,
      STK(spnm).b=RBF((numorg+1):rdim,(RFn+1):(RFn+Bn));
    end
    if computeH,
      STK(spnm).rows=rowidx((numorg+1):rdim);
    end
  else
    STK(spnm).A=[]; 
    STK(spnm).cols=[];
    STKm(spnm)=0;
    STKn(spnm)=0;
  end
  
  if computeH,
    if FNTm>rdim,
      rowperm(last-(FNTm-rdim)+1:last)=rowidx(rdim+1:FNTm);
      last=last-(FNTm-rdim);
    end
  end
  
  lmr=lmr+Am;
  
end

% Up to now the transpose of R has been stored ...

R=R';

if ~computeH & computeB
  H=C;
  clear C
end
