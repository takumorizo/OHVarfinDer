#! /bin/bash
#$ -S /bin/sh
#$ -cwd
: <<'#__CO__'
#__CO__

OHVarTumorMinDepth=12
OHVarTumorMinRate=0.05
OHVarTumorMinObsNum=4

OHVarNormalMinDepth=12
OHVarNormalMaxRate=0.1

OHVarSNPMinDepth=15
OHVarSNPMinRate=0.3
OHVarSNPMinObsNum=10

# If you can use single reads only, comment out this line.
# isSingle=--singleReads
isSingle=

ohvar2_mutGammaF=10.0,1.0
ohvar2_mutGammaH=5.0,5.0,1.0
ohvar2_mutAlphaL=1.0,100.0
ohvar2_mutAlphaH=1.0,100.0
ohvar2_mutGammaEH=5.0,5.0
ohvar2_mutAlphaS=1.0,100.0

ohvar2_mutAlphaB_E=10.0,10.0
ohvar2_mutAlphaB_W=1.0,1.0

ohvar2_errAlphaL_E=2.0,30.0
ohvar2_errAlphaH_E=2.0,30.0
ohvar2_errGammaEH_E=5.0,5.0
ohvar2_errAlphaS_E=2.0,30.0
ohvar2_errAlphaB_E=0.05,0.05

ohvar2_errAlphaL_W=2.0,30.0
ohvar2_errAlphaH_W=2.0,30.0
ohvar2_errGammaEH_W=5.0,5.0
ohvar2_errAlphaS_W=2.0,30.0
ohvar2_errAlphaB_W=0.5,0.5

