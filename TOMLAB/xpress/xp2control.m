% Xpress-MP MEX-interface internal callback routine
%
% Creates a global structure variable xpControlVariables
% where the fields corresponds to the Xpress-MP control parameter names
% as given in Xpress-Optimizer Reference Manual Release 13, section 7.
%
% The control parameter names always start with "XPRS_" .
% These 5 characters are NOT included in the field names!
% 
% function xp2control(xpicv,xpdcv,xpccv1,xpccv2,xpccv3,xpccv4,xpccv5,xpccv6)
%
% INPUT:  
%  xpicv  Vector of doubles with Xpress-MP Integer Control Variables
%  xpdcv  Vector of doubles with Xpress-MP Double  Control Variables
%  xpccv1 String with 1st        Xpress-MP String  Control Variable
%  xpccv2 String with 2nd        Xpress-MP String  Control Variable
%  xpccv3 String with 3rd        Xpress-MP String  Control Variable
%  xpccv4 String with 4th        Xpress-MP String  Control Variable
%  xpccv5 String with 5th        Xpress-MP String  Control Variable
%  xpccv6 String with 6th        Xpress-MP String  Control Variable
%
% OUTPUT:  
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomlab.biz
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 2004.0.0 $
% Written Oct 25, 1999.   Last modified Dec 11, 2001.
%

function xp2control(xpicv,xpdcv,xpccv1,xpccv2,xpccv3,xpccv4,xpccv5,xpccv6)

global xpControlVariables

%disp('xp2control')
%pause

PV.AUTOPERTURB   = xpicv(27);
PV.BACKTRACK     = xpicv(46);
PV.BARDUALSTOP   = xpdcv(16);
PV.BARGAPSTOP    = xpdcv(15);
PV.BARITERLIMIT  = xpicv(34);
PV.BARMEMORY     = xpicv(40);
PV.BARORDER      = xpicv(1);
PV.BAROUTPUT     = xpicv(58);
PV.BARPRIMALSTOP = xpdcv(19);
PV.BARTHREADS    = xpicv(62);
PV.BARSTEPSTOP   = xpdcv(17);
PV.BIGM          = xpdcv(2);
PV.BIGMMETHOD    = xpicv(14);
PV.BREADTHFIRST  = xpicv(35);
PV.CACHESIZE     = xpicv(2);
PV.CHOLESKYALG   = xpicv(61);
PV.CHOLESKYTOL   = xpdcv(18);
PV.COVERCUTS     = xpicv(50);
PV.CPKEEPALLCUTS = xpicv(19);
PV.CPMAXCUTS     = xpicv(42);
PV.CPMAXELEMS    = xpicv(41);
PV.CPUTIME       = xpicv(3);
PV.CRASH         = xpicv(13);
PV.CROSSOVER     = xpicv(4);
PV.CSTYLE        = xpicv(15);
PV.CUTDEPTH      = xpicv(5);
PV.CUTFREQ       = xpicv(6);
PV.CUTSTRATEGY   = xpicv(7);
PV.DEFAULTALG    = xpicv(12);
PV.DEGRADEFACTOR = xpdcv(11);
PV.DENSECOLLIMIT = xpicv(8);
PV.ELIMTOL       = xpdcv(21);
% NOT DECLARED!!! (But in the manual) PV.ERRIGNORE     = xpicv(67);
PV.ETATOL        = xpdcv(22);
PV.EXTRACOLS     = xpicv(44);
PV.EXTRAELEMS    = xpicv(53);
PV.EXTRAMIPENTS  = xpicv(48);
PV.EXTRAPRESOLVE = xpicv(47);
PV.EXTRAROWS     = xpicv(54);
PV.FEASTOL       = xpdcv(31);
PV.GOMCUTS       = xpicv(51);
PV.INVERTFREQ    = xpicv(25);
PV.INVERTMIN     = xpicv(26);
PV.KEEPBASIS     = xpicv(18);
PV.KEEPMIPSOL    = xpicv(29);
PV.KEEPNROWS     = xpicv(20);
PV.LPITERLIMIT   = xpicv(31);
PV.LPLOG         = xpicv(32);
PV.MARKOWITZTOL  = xpdcv(28);
PV.MATRIXTOL     = xpdcv(23);
PV.MAXIIS        = xpicv(49);
PV.MAXMIPSOL     = xpicv(38);
PV.MAXNODE       = xpicv(36);
PV.MAXPAGELINES  = xpicv(52);
PV.MAXSLAVE      = xpicv(37);
PV.MAXTIME       = xpicv(39);
PV.MIPABSCUTOFF  = xpdcv(3);
PV.MIPABSSTOP    = xpdcv(6);
PV.MIPADDCUTOFF  = xpdcv(1);
PV.MIPLOG        = xpicv(11);
PV.MIPPRESOLVE   = xpicv(16);
PV.MIPRELCUTOFF  = xpdcv(9);
PV.MIPPRELSTOP   = xpdcv(7);
PV.MIPTARGET     = xpdcv(5);
PV.MIPTOL        = xpdcv(24);
PV.MPSBOUNDNAME = xpccv2;
PV.MPSECHO       = xpicv(59);
PV.MPSERRIGNORE  = xpicv(24);
PV.MPSFORMAT     = xpicv(10);
PV.MPSNAMELENGTH = xpicv(43);
PV.MPSOBJNAME   = xpccv4;
PV.MPSRANGENAME = xpccv6;
PV.MPSRHSNAME   = xpccv5;
PV.NODESELECTION = xpicv(45);
PV.OMNIDATANAME = xpccv3;
PV.OMNIFORMAT    = xpicv(57);
PV.OPTIMALITYTOL = xpdcv(20);
PV.OUTPUTLOG     = xpicv(60);
PV.OUTPUTMASK   = xpccv1;
PV.OUTPUTTOL     = xpdcv(26);
PV.PENALTY       = xpdcv(8);
PV.PERTURB       = xpdcv(10);
PV.PIVOTTOL      = xpdcv(25);
PV.PPFACTOR      = xpdcv(12);
PV.PRESOLVE      = xpicv(22);
PV.PRICINGALG    = xpicv(28);
PV.PSEUDOCOST    = xpdcv(13);
PV.RECEXPAND     = xpdcv(4);
PV.RECMAXPASSES  = xpicv(63);
PV.RECSHRINK     = xpdcv(14);
PV.RECSTEPLENGTH = xpdcv(32);
PV.RECSTOP       = xpdcv(27);
PV.REFACTOR      = xpicv(9);
PV.REL10STYLE    = xpicv(21);
PV.RELPIVOTTOL   = xpdcv(29);
PV.SBBEST        = xpicv(64);
PV.SBITERLIMIT   = xpicv(65);
PV.SCALING       = xpicv(23);
PV.SOLUTIONFILE  = xpicv(17);
PV.SOSREFTOL     = xpdcv(30);
PV.TRACE         = xpicv(30);
PV.TREECOVERCUTS = xpicv(55);
PV.TREEGOMCUTS   = xpicv(56);
PV.VARSELECTION  = xpicv(33);
PV.VERSION       = xpicv(66);


xpControlVariables=PV;

% MODIFICATION LOG:
%
% 991022 hkh  Written
% 000709 hkh  Modified for Release 12
% 011209 hkh  Modified for Release 13
% 011210 hkh  Improving comments
% 011212 hkh  Alphabetic order. Add four new controls.

