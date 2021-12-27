# NIRS_dyadic_analyses
Currently the script looks at two types of correlation: matched channels (i.e., channel 2 in subject 1 and subject 2) and area specific (i.e., mPFC in subject 1 and subject 2).

Make sure data is preprocessed before running (see:https://github.com/abinnquist/fNIRSpreProcessing). Data structure should be as follows:
dataFolder/dyads/subjects/scans/fNIRS_scandata.

Analyses require three helper scripts provided in the folder (compileNIRSdata, imageNIRS & fdr_bky.mat) and two .mat files stored in each study folder (mni_coordinates & dyad numbers).

Visualizing/imaging the analyses will require spm12 (see: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and xjview (see: https://www.alivelearn.net/xjview/). 
I am currently using Surf_Ice (see: https://www.nitrc.org/projects/surfice/) for cortex rendering of the hdr file created by the above two programs.

More notes to come.
