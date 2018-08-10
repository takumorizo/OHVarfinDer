#! /bin/bash

maxInsertSize=1600
tumorMinDepth=12
tumorMinObsRate=0.05
tumorMaxObsRate=1.0
tumorMinObsNum=4

normalMinDepth=12
normalMaxObsRate=0.1
normalMaxObsNum=1000000

heteroSNPMinDepth=15
heteroSNPMinObsRate=0.3
heteroSNPMinObsNum=10

triAlleleMinObsRate=0.03
triAlleleMinObsNum=5

averageMapPhredQualThreshold=0
softClipPosessionFreqThreshold=1.0

tumorMinAvgBaseQuality=0
normalMinAvgBaseQuality=0
heteroSNPMinAvgBaseQuality=0

# If you can use single reads only, comment out this line.
# isSingle=--singleReads
isSingle=

pileUpBufferSize=4000000

ohvar_mutGammaF=10.0,1.0
ohvar_mutGammaH=5.0,5.0,1.0
ohvar_mutAlphaL=1.0,100.0
ohvar_mutAlphaH=1.0,100.0
ohvar_mutGammaEH=5.0,5.0
ohvar_mutAlphaS=1.0,100.0

ohvar_mutAlphaB_E=10.0,10.0
ohvar_mutAlphaB_W=1.0,1.0

ohvar_errAlphaL_E=2.0,30.0
ohvar_errAlphaH_E=2.0,30.0
ohvar_errGammaEH_E=5.0,5.0
ohvar_errAlphaS_E=2.0,30.0
ohvar_errAlphaB_E=0.05,0.05

ohvar_errAlphaL_W=2.0,30.0
ohvar_errAlphaH_W=2.0,30.0
ohvar_errGammaEH_W=5.0,5.0
ohvar_errAlphaS_W=2.0,30.0
ohvar_errAlphaB_W=0.5,0.5

