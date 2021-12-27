# NIRS_dyadic_analyses
Script for dyadic correlation of conversational data

Make sure data is preprocessed before running (see:https://github.com/abinnquist/fNIRSpreProcessing). Data structure should be as follows:
dataFolder/dyads/subjects/scans/fNIRS_scandata.

Analyses only require two helper scripts (compileNIRSdata & fdr_bky.mat) and two .mat files stored in each study folder (mni_coordinates & dyad numbers).

Visualizing/imaging the analyses will require spm12 (see: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and xjview (see: https://www.alivelearn.net/xjview/). 
I am currently using Surf_Ice (see: https://www.nitrc.org/projects/surfice/) for cortex rendering of the hdr file created by the above two programs.

More notes to come.
