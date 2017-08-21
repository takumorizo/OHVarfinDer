#! /bin/bash

#
# Filter condition for candidate mutation, or heterozygous SNPs
#
tumorMinDepth=10
tumorMinObsRate=0.05
tumorMaxObsRate=1.0
tumorMinObsNum=4
tumorMinAvgBaseQuality=25

normalMinDepth=10
normalMaxObsRate=0.1
normalMaxObsNum=1
normalMinAvgBaseQuality=25

heteroSNPMinDepth=15
heteroSNPMinObsRate=0.3
heteroSNPMinObsNum=10
heteroSNPMinAvgBaseQuality=25

triAlleleMinObsRate=0.03
triAlleleMinObsNum=2


# If you can use single reads only, comment out this line.
# isSingle=--singleReads
isSingle=

maxInsertSize=1600

pileUpBufferSize=4000000

#
# Hyperparameters @ OHVarfinDer model
#
ohvar2_mutGammaF=5.0,1.0
ohvar2_mutGammaH=2.5,2.5,1.0
ohvar2_mutAlphaL=0.1,10.0
ohvar2_mutAlphaH=0.1,10.0
ohvar2_mutGammaEH=1.0,1.0
ohvar2_mutAlphaS=0.1,10.0

ohvar2_mutAlphaB_E=10.0,10.0
ohvar2_mutAlphaB_W=1.0,1.0

ohvar2_errAlphaL_E=1.0,10.0
ohvar2_errAlphaH_E=1.0,10.0
ohvar2_errGammaEH_E=1.0,1.0
ohvar2_errAlphaS_E=1.0,10.0
ohvar2_errAlphaB_E=0.05,0.05

ohvar2_errAlphaL_W=1.0,10.0
ohvar2_errAlphaH_W=1.0,10.0
ohvar2_errGammaEH_W=1.0,1.0
ohvar2_errAlphaS_W=1.0,10.0
ohvar2_errAlphaB_W=0.5,0.5
