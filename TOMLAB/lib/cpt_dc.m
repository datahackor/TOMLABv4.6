% cpt_c.m 
%
% function [dc,nn]=cpt_dc(x, Prob)
%
% cpt_dc is called by CONOPT to evaluate the Jacobian of the 
% nonlinear constraints and the objective function at the point x.  
%
% Automatically handles the extended constraint Jacobian when
% both upper and lower bounded constraints are present.  
%
% Anders Goran, Tomlab Optimization Inc, E-mail: anders@tomlab.biz
% Copyright (c) 2003-2004 by Tomlab Optimization Inc., $Release: 4.2.0$
% Written June 2, 2003.    Last modified June 3, 2003.

function [dc,nn] = cpt_dc(x,Prob)

cexidx = Prob.cexidx;
m2     = Prob.m2;

g  = feval('nlp_g',x,Prob);
dc = feval('nlp_dc',x,Prob);

dc = sparse( [ dc ; dc(cexidx(m2+1:end),:) ; g' ] )';

% Number of nonzeros is also returned
nn = nnz(dc);

%full(dc)
%nn
