#! /bin/bash

maxInsertSize=1600
tumorMinDepth=100
tumorMinObsRate=0.02
tumorMaxObsRate=0.07
tumorMinObsNum=3
normalMinDepth=100
normalMaxObsRate=0.01
normalMaxObsNum=1
heteroSNPMinDepth=30
heteroSNPMinObsRate=0.3
heteroSNPMinObsNum=10

triAlleleMinObsRate=0.03
triAlleleMinObsNum=2
tumorMinAvgBaseQuality=25
normalMinAvgBaseQuality=25
heteroSNPMinAvgBaseQuality=25

# If you can use single reads only, comment out this line.
# isSingle=--singleReads
isSingle=

ohvar_mutGammaF=5.0,1.0
ohvar_mutGammaH=2.5,2.5,1.0
ohvar_mutAlphaL=0.1,10.0
ohvar_mutAlphaH=0.1,10.0
ohvar_mutGammaEH=1.0,1.0
ohvar_mutAlphaS=0.1,10.0

ohvar_mutAlphaB_E=10.0,10.0
ohvar_mutAlphaB_W=1.0,1.0

ohvar_errAlphaL_E=1.0,10.0
ohvar_errAlphaH_E=1.0,10.0
ohvar_errGammaEH_E=1.0,1.0
ohvar_errAlphaS_E=1.0,10.0
ohvar_errAlphaB_E=0.05,0.05

ohvar_errAlphaL_W=1.0,10.0
ohvar_errAlphaH_W=1.0,10.0
ohvar_errGammaEH_W=1.0,1.0
ohvar_errAlphaS_W=1.0,10.0
ohvar_errAlphaB_W=0.5,0.5

