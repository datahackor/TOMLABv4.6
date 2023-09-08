% TOMLAB - Large-Scale Optimization
% Version 4.6 (R4.6.0) 13-Jan-2005
% 
% -----------------------------------
% Files in the main directory, tomlab
% -----------------------------------
% contents.m    This file
% startup.m     Using Matlab, go to Tomlab directory, run startup.
% install.m     Help in installing Tomlab
% changes.m     Log of changes in Tomlab.
% tomlablic     License file (installed separately)
%
% ----------------------------------
% Driver routines, in directory base
% ----------------------------------
% tomRun        General driver routine calling any Tomlab solver
% tomSolve      General subproblem driver routine calling any Tomlab solver
%               No check is made on input problem structure.
%               Global variables are saved and restored.
% -------------------------------------------
% Files in the BASE MODULE, in directory base
% -------------------------------------------
%
% tomlabInit    Tomlab global variable initialization (not necessary to run)
% tomlabGUI     Graphical User Interface (GUI) for Optimization
% tomGUI        Short name for tomlabGUI
% tomMenu       General Tomlab menu program, similar to GUI
% tomRemote     Menu program for remote Unix use
% tomHelp       Help on the nonGUI use of Tomlab
% tomlab.bib    Bibtex file with Tomlab references
% checkdll      Checks the licenses
%
% -----------------------------------------------
% The TOM solvers, in directory base, mainly written in Matlab code.
% -----------------------------------------------
% clsSolve      Constrained nonlinear least squares solver:
%               Active set strategy. As search step algorithm using:
%               Gauss-Newton with subspace minimization, Fletcher-Xu hybrid
%               method, Al-Baali-Fletcher hybrid method and Huschens method
% conSolve      SQP algorithm. Schittkowski with Augmented Lagrangian
%               Han-Powell (Quasi-Newton update)
% cutplane      Gomorov's cutting plane algorithm for Mixed-Integer Programs
% DualSolve     Dual simplex algorithm with three selection rules. 
% expSolve      Solve exponential fitting problems
% glbFast       glbSolve DIRECT algorithm in faster Fortran version
% glbSolve      Global Optimization algorithm DIRECT by Don Jones et.al.
% glcCluster    Hybrid of glcFast / cluster alg / local solver
% glcFast       glcSolve constrained DIRECT algorithm in faster Fortran version
% glcSolve      Constrained Mixed-Integer Global Optimization, Constr.DIRECT
% goalSolve     Multi-objective Goal Attainment
% infSolve      Sparse constrained minimax solver
% L1LinSolve    Finds a linearly constrained L1 solution
% L1Solve       Sparse constrained L1 solver
% lpSolve       Simplex algorithm for general LP, structure input. 
%               Solves Phase 1 and Phase 2 problems
% mipSolve      Branch & Bound for Mixed-Integer Programs (MIP)
% nlpSolve      Filter SQP algorithm by Fletcher-Leyffer. Calls qpPhase1.
% qpSolve       Solve QP using active set method.
% slsSolve      Sparse Least Squares solver (constrained)
% sTrustr       A Structural Trust Region algorithm for unconstrained 
%               optimization (Conn, Gould, Sartenaer, Toint). Calls itrr
% stepexp       Stepwise solve exponential fitting problems
% ucSolve       Unconstrained optimization solver, handling simple bounds
%               on the parameters. Algorithms: Newton, BFGS, Inverse BFGS,
%
% -------------------------------------------------------------------
% The CGO solvers, in directory cgo, mainly written in Matlab code.
% -------------------------------------------------------------------
% rbfSolve      Costly Global Optimization (in /CGO), RBF interpolation
% ego           Efficient Global Optimization (EGO) ( in /CGO)
%
% -----------------------------------------------
% SOLVERS compiled in PC MEX DLLs or Unix MEX libs
% -----------------------------------------------
% BARNLP        Sparse barrier solver                  (/BARNLP) 
% bqpd          Dense or sparse (Prob.LargeScale == 1) QP/LP solver (/MINLP)
% conopt        Dense or sparse NLP                    (/CONOPT)
% CPLEX         Sparse LP, QP, MILP, MIQP solver       (/CPLEX)
% dfnlp         Nonlinear data fitting                 (/NLPQL)
% filterSQP     Dense or sparse (Prob.LargeScale == 1) NLP solver (/MINLP)
% LGO           Global optimization package            (/LGO)
% lsei          Linear least squares, with equality and inequality constraints
% lpopt         Dense LP solver                        (/MINOS)
% lssol         Dense QP/LP or least squares solver    (/SOL, /NPSOL)
% minlpbb       Dense or sparse (Prob.LargeScale == 1) MINLP solver (/MINLP)
% minos         Sparse LP, nonconvex QP or NLP solver  (/MINOS)
% miqpbb        Dense or sparse (Prob.LargeScale == 1) MIQP/MILP solver (/MINLP)
% nlpql         Dense SQP solver                       (/NLPQL)
% nlpjob        Multi criteria optimization            (/NLPQL)
% nlssol        Dense constrained least squares solver (/SOL, /NPSOL)
% npsol         Dense NLP solver                       (/SOL, /NPSOL)
% oqnlp         Dense or sparse NLP/MINLP              (/OQNLP) 
% pensdp        Sparse semidefinite programming solver (/PENSDP)
% pdco          Primal Dual Convex Optimization, using Tlsqr (Base Module)
% pdsco         Primal Dual Separable Convex Optimization (Base Module)
% qld           Dense convex QP/LP solver   
% qpopt         Dense nonconvex QP/LP solver           (/MINOS)
% snopt         Sparse NLP solver                      (/SOL, /SNOPT)
% SPRNLP        Sparse SQP solver                      (/SPRNLP)
% sqopt         Sparse convex QP solver                (/SOL, /SNOPT)
% tomsol        Code for speedup of computations, and some subalgorithms
% xpress-mp     Sparse LP, QP, MILP, MIQP solver       (/Xpress)
% conopt        Sparse and dense NLP                   (/CONOPT)
% Tfmin         Tomlab fmin: Finds a minimum of f(x) in [xL,xU], x 1-dim.
%               Same algorithm as Matlab fmin, but faster, using MEX-interface
% Tfzero        Finding a zero in an interval
% Tknitro       Dense or sparse interior point NLP     (/KNITRO)
% Tlsqr         Large, sparse unsymmetric equations or linear least squares
% Tnnls         Tomlab NNLS, Solves nonnegative least squares, 
%               also with linear equality constraints
% XA            Sparse LP, MILP, QP                    (/XA)
%
%
% -------------------------------------------------------------------
% Utilities to help defining problems in the Tomlab (TQ) format
% -------------------------------------------------------------------
% probAssign    Setup a Prob structure for a problem of certain problem type
% lpAssign      Define a Linear Programming problem 
% qpAssign      Define a Quadratic Programming problem 
% conAssign     Define a NonLinear Programming problem (constrained or not) 
% mipAssign     Define a mixed-integer programming problem 
% miqpAssign    Define a mixed-integer quadratic programming problem
% clsAssign     Define a nonlinear least squares problem (constrained or not)
% llsAssign     Define a linear least squares problem (constrained or not)
% glcAssign     Define a global optimization problem (constrained or not)
% sdpAssign     Define a semidefinite program
% bmiAssign     Define a bilinear semidefinite program
% minlpAssign   Define a mixed-integer nonlinear (MINLP) program
% simAssign     Both the function and the constraints are computed
% tomFiles      Set names of the m-files to be used into structure Prob
%               Only needed if using the general probAssign routine
%
% -----------------------------------------------------------------------
% Utilities to help defining problems in the Tomlab Init File (IF) format
% -----------------------------------------------------------------------
% newInitFile   Initialize the creation of a new Init File.
% addProb       Add one problem structure to the new Init File
% makeInitFile  Finish the creation of a new Tomlab Init File. A single call to
%               makeInitFile possible if 1-5 problem structures used.
%
% --------------------------------------------------------------------------
% Utilities to add/delete problem Init Files (IF) for use with GUI and Menu
% --------------------------------------------------------------------------
%
% CreateTomProb     Creates a file TomlabProblem.mat with the standard 
%                   predefined test problems in Tomlab
% AddProblemFile    Add a problem definition file to list of predefined files
% DeleteProblemFile Delete problem definition file from list of predefined files
%
% -----------------
% Testing utilities
% -----------------
% pretest       Test of presolve analysis, calls preSolve.m
% runtest       Run test of solver for sequence of problems
% systest       Run test of many solvers for sets of problems
% testtom       Script to run Tomlab system tests
% ---------------
% Other utilities
% ---------------
% penfeas_lmi   Test feasibility for semidefinite program with lmi constraints
% penfeas_bmi   Test feasibility for semidefinite program with bmi constraints
% preSolve      Presolve analysis on linear constraints, algorithm by Gondzio
% GetSolver     Returns the default solver for all types of optimization
% lp2tp         Convert LP to transportation problem
% tp2lp         Convert transportation problem to LP problem
% tp2np         Convert transportation problem to network problem
%
% ---------------------------
% END OF TOMLAB MAIN ROUTINES
% ---------------------------
%
% ---------------------------
% Sub directories in Tomlab
% ---------------------------
%
% admat         Dummys for ADMAT TB, automatic differentiation
% ampl          Interface to problems formulated in AMPL
% base          Files in the Base Module, and other utilities
% cgo           Files for Costly Global Optimization (CGO)
% cplex         The CPLEX solver package
% examples      Test routines and many example files
% ldo           Linear and discrete utility algorithms
% ldodemo       Demonstration routines for the LDO algorithms
% lib           Library functions
% mad           Dummys for MAD, MATLAB Automatic Differentiation
% matlab5.1     For users of old Matlab versions, 5.1 and 5.0
% mex           Other solvers called using MEX-file interfaces.
% optim         Tomlab versions of the Optimization TB 2.0 solvers
% optim1.x      Code to run optimization solvers in MathWorks Optimization TB
% quickguide    Guide for getting started with TOMLAB
% shared        TOMLAB shared files location
% splines       Dummys for SPLINE TB, the spline toolbox
% sqr2          Sparse multifrontal QR factorization in Matlab. Authors: Mikael 
%               Adlers, Linkoping University, Sweden and Pontus Matstoms, VTI 
%               Linkoping, Sweden. http://www.mai.liu.se/~milun/sls/
% testprob      A collection of predefined test problems, and test routines
% tsplib        A collection of test problems for Travelling Salesman 
% unix5.1       For users of old Matlab versions on Unix
% usersguide    A directory with the test files shown in the User's Guide
% xa            The XA solver package
% xpress        The Xpress solver package
%
