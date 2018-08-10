#! /bin/bash

maxInsertSize=800
tumorMinDepth=10
tumorMinObsRate=0.05
tumorMaxObsRate=1.0
tumorMinObsNum=4

normalMinDepth=10
normalMaxObsRate=0.1
normalMaxObsNum=1

heteroSNPMinDepth=15
heteroSNPMinObsRate=0.3
heteroSNPMinObsNum=10

triAlleleMinObsRate=0.03
triAlleleMinObsNum=2

averageMapPhredQualThreshold=15
softClipPosessionFreqThreshold=0.25

tumorMinAvgBaseQuality=25
normalMinAvgBaseQuality=25
heteroSNPMinAvgBaseQuality=25

# If you can use single reads only, comment out this line.
# isSingle=--singleReads
isSingle=

pileUpBufferSize=4000000

ohvar_mutGammaF=5.0,1.0
ohvar_mutGammaH=2.5,2.5,1.0
ohvar_mutAlphaL=0.1,10.0
ohvar_mutAlphaH=0.1,10.0
ohvar_mutGammaEH=1.0,1.0
ohvar_mutAlphaS=0.1,10.0

ohvar_mutAlphaB_E=10.0,10.0
ohvar_mutAlphaB_W=5.0,5.0

ohvar_errAlphaL_E=1.0,10.0
ohvar_errAlphaH_E=1.0,10.0
ohvar_errGammaEH_E=1.0,1.0
ohvar_errAlphaS_E=1.0,10.0
ohvar_errAlphaB_E=0.05,0.05

ohvar_errAlphaL_W=1.0,10.0
ohvar_errAlphaH_W=1.0,10.0
ohvar_errGammaEH_W=1.0,1.0
ohvar_errAlphaS_W=1.0,10.0
ohvar_errAlphaB_W=0.2,0.2

