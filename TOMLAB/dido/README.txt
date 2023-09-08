**************************************************************************
DIDO READ ME
**************************************************************************

*----------------------------*
Last Revised 25 JANUARY 2004
Copyright (c) by I. M. Ross
*----------------------------*

=========================================
This is Beta version of DIDO Version PR.1
=========================================

DIDO Users are automatically bound to the terms of usage as described in this
file and the User's Manual.


*****************************************************************************
********** WARNING *** WARNING *** WARNING *** WARNING *** WARNING **********

DIDO is not for distribution without the explicit permission of 
I. M. Ross, Monterey, CA 93940.

You may not alter DIDO or reverse-engineer it in any form WHATSOVER!  You may
be subject to severe criminal penalties in addition to a fine of upto $500,000 
for violating this agreement.  

Reference: THE DIGITAL MILLENNIUM COPYRIGHT ACT OF 1998

********** WARNING *** WARNING *** WARNING *** WARNING *** WARNING **********
*****************************************************************************




=========================
INSTALLATION INSTRUCTIONS
=========================



STEP 1:
------- 

If you have TOMLAB installed, move to STEP 2; otherwise, install TOMLAB by 
following the instructions that came with it. Verify TOMLAB is working
properly.


STEP 2:
-------

If you have a compressed version of DIDO, uncompress it.  It will automatically
uncompress with the right folder structure.  If you have an uncompressed file, 
install DIDO by simply copying the folder DIDO_PR.x (and all subfolders) to 
your hard drive.  


STEP 3:
-------

Right click on the DIDO ICON and click on PROPERTIES.
Click on the FIND TARGET button and select MATLAB.EXE.  The file MATLAB.EXE
with its full pathname must show up on the TARGET BOX.  You may also type the
full pathname here.

Next,
Under the box START IN, select the folder DIDO_PR.x by typing the full path
name.  Ex: C:\ ....\DIDO_PR.x


STEP 4:
-------

You may now move the DIDO ICON anywhere (e.g. your desktop).  Double-Click the
DIDO ICON.  This should start MATLAB and initalize DIDO if you correctly performed 
STEP 3.


STEP 5:
-------

At the MATLAB prompt, >>, type,

>> TestDIDO


STEP 6:
------

A successful installation should produce the final output on the screen:

CONGRATULATIONS! DIDO TEST WAS SUCCESSFUL.

* If you do not get this message, go back to STEP 3.

* If you get the error message: 
MATLAB Segmentation Violation detected, go back to STEP 1

********************************************

STEP 7 (optional):
-------

You may perform a further test by typing at the MATLAB prompt, >>

>> LanderMain

A successful installation should produce a pretty graph.

============================================

DO NOT TAMPER WITH THE FOLDER DIDO or DIDO_Complete.  

This includes reading or writing files to it. 

See Warnings at the top of this file.

You may freely modify the ForUser Folder 

DO NOT use the keyword DIDO anywhere in your code
=============================================
*********************************************

-------------------
NOTES AND CREDITS:
-------------------

Queen Dido was the first person to solve an isoperimetric problem of the 
calculus of variations circa 800 BC long before the invention of calculus!

DIDO is not an acronym for anything.

Future versions of the code will enable solving VASTLY more complex problems.
And, faster too! 

----------------------------------------------------------------------------
This theoretical foundations for the methods used in this code is partially
based on the publications of I. M. Ross and F. Fahroo.  
----------------------------------------------------------------------------

________________________
FOR FURTHER INFORMATION:
------------------------
Contact 

* TBD


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) by I. M. Ross.  All rights reserved.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EOF
*************************************************************************
