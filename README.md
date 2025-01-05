# Comparison_1D_TDEM_Forwards
Codes required to perform the analysis. Please note that if using AarhusInv the programme can be downloaded from https://hgg.au.dk/software/aarhusinv, and the installation should be made to the folder "Comparison1DForwards/AarhusInv/" such that it contains the dll's and AarhusInv64.exe & AarhusInvLic.exe files. Please use the "Aarhusinv.con" file already in the folder.

Steps:
1) Open the MATLAB script "ForwardsAnalysisComparison.m"
2) In the block of code from line 47 to 52 please choose which forward functions to compare.
3) In the block of code from line 53 to 57 please choose which analyses to perform (precision analysis 2 requires at least 2 forward functions in the previous block)
4) In the block of code from line 58 to 96 it is possible to choose the settings for the analysis, such as the number of points in the curves and the tested values for the parameters in each analysis.

