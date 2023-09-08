% Xpress-MP MEX-interface internal callback routine
%
% Creates a global structure variable xpProblemAttrib
% where the fields corresponds to the Xpress-MP problem attributes
% as given in Xpress-Optimizer Reference Manual Release 13, section 8.
% 
% The problem attribute names always start with "XPRS_" .
% These 5 characters are NOT included in the field names!
%
% function xp2problem(xpipv,xpdpv,xpcpv1,xpcpv2,xpcpv3,xpcpv4,xpcpv5)
%
% INPUT:  
%  xpipv  Vector of doubles with Xpress-MP Integer Problem Attributes
%  xpdpv  Vector of doubles with Xpress-MP Double  Problem Attributes
%  xpcpv1 String with 1st        Xpress-MP String  Problem Attribute
%  xpcpv2 String with 2nd        Xpress-MP String  Problem Attribute
%  xpcpv3 String with 3rd        Xpress-MP String  Problem Attribute
%  xpcpv4 String with 4th        Xpress-MP String  Problem Attribute
%  xpcpv5 String with 5th        Xpress-MP String  Problem Attribute
%
% OUTPUT:  
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 2004.0.0 $
% Written Oct 22, 1999.   Last modified Dec 11, 2001.
%

function xp2problem(xpipv,xpdpv,xpcpv1,xpcpv2,xpcpv3,xpcpv4,xpcpv5)

global xpProblemAttrib

%disp('xpProblemAttrib')
%pause

PV.ACTIVENODES   = xpipv(5);
PV.BARAASIZE     = xpipv(11);
PV.BARCROSSOVER  = xpipv(7);
PV.BARDENSECOL   = xpipv(8);
PV.BARDUALINF   = xpdpv(1);
PV.BARDUALOBJ   = xpdpv(2);
PV.BARITER       = xpipv(9);
PV.BARLSIZE      = xpipv(11);
PV.BARPRIMALINF = xpdpv(4);
PV.BARPRIMALOBJ = xpdpv(5);
PV.BARSTOP      = xpdpv(6);
PV.BESTBOUND    = xpdpv(3);
PV.BOUNDNAME  =xpcpv1;
PV.COLS          = xpipv(12);
PV.CUTS          = xpipv(14);
PV.DUALINFEAS    = xpipv(21);
PV.ELEMS         = xpipv(16);
PV.ERRORCODE     = xpipv(1);
PV.IIS           = xpipv(31);
PV.LPOBJVAL     = xpdpv(10);
PV.LPSTATUS      = xpipv(33);
PV.MATRIXNAME =xpcpv2;
PV.MIPENTS       = xpipv(34);
PV.MIPINFEAS     = xpipv(3);
PV.MIPOBJVAL    = xpdpv(7);
PV.MIPSOLNODE    = xpipv(22);
PV.MIPSOLS       = xpipv(19);
PV.MIPSTATUS     = xpipv(2);
PV.NAMELENGTH    = xpipv(6);
PV.NODEDEPTH     = xpipv(15);
PV.NODES         = xpipv(23);
PV.OBJFIXED     = xpdpv(8);
PV.OBJNAME    =xpcpv3;
PV.OBJRHS       = xpdpv(9);
PV.OBJSENSE     = xpdpv(12);
PV.PARENTNODE    = xpipv(24);
PV.PRESOLVESTATE = xpipv(32);
PV.PRIMALINFEAS  = xpipv(18);
PV.QELEMS        = xpipv(26);
PV.RANGENAME  =xpcpv5;
PV.RECCONVERGE   = xpipv(25);
PV.RHSNAME    =xpcpv4;
PV.ROWS          = xpipv(27);
PV.SIMPLEXITER   = xpipv(4);
PV.SETMEMBERS    = xpipv(30);
PV.SETS          = xpipv(29);
PV.SPARECOLS     = xpipv(13);
PV.SPAREELEMS    = xpipv(20);
PV.SPAREMIPENTS  = xpipv(17);
PV.SPAREROWS     = xpipv(28);
PV.SUMPRIMALINF = xpdpv(11);
% PV.VERSION       = xpipv(35);


xpProblemAttrib=PV;

% MODIFICATION LOG:
%
% 991022 hkh Written
% 010710 hkh Revised for Release 12
% 011209 hkh Revised for Release 13
% 011210 hkh Changed name from xpProblemVariables to xpProblemAttrib
% 011212 hkh Changed to alphabetic order. 

