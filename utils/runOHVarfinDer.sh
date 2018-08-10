#! /bin/bash

readonly DIR=`dirname ${0}`

readonly REF=$1
readonly TUMOR=$2
readonly NORMAL=$3
readonly OUTPUTDIR=$4
readonly REGION=${5-""}
readonly PARAM=${6-${DIR}/params/params.sh}
readonly FILTER=${7-${DIR}/params/filterList.sh}
source ${PARAM}
source ${FILTER}

echo "
Reference  : ${REF}
Tumor bam  : ${TUMOR}
Normal bam : ${NORMAL}
OUTPUTDIR  : ${OUTPUTDIR}
parameter  : ${PARAM}
filter     : ${FILTER}
region     : ${REGION}
"


check_mkdir()
{
if [ -d $1 ]
then
    echo "$1 exists."
else
    echo "$1 does not exits."
    mkdir -p $1
fi
}

check_error()
{
if [ $1 -ne 0 ]; then
    echo "FATAL ERROR: pipeline script"
    echo "ERROR CODE: $1"
    exit $1
fi
}

check_file_exists()
{
if [ -f $1 ]; then
  echo "$1 exists."
  return
fi
echo "$1 does not exists."
exit 1
}

echo "check_mkdir ${OUTPUTDIR}"
check_mkdir ${OUTPUTDIR}

echo "${DIR}/../bin/ohvarfinder \
--algorithm OHVarfinDer \
-f ${REF}    \
-a ${TUMOR}  \
-b ${NORMAL} \
-o ${OUTPUTDIR}/output \
--maxInsertSize=${maxInsertSize} \
--tumorMinDepth=${tumorMinDepth} \
--tumorMinObsRate=${tumorMinObsRate} \
--tumorMaxObsRate=${tumorMaxObsRate} \
--tumorMinObsNum=${tumorMinObsNum} \
--normalMinDepth=${normalMinDepth} \
--normalMaxObsRate=${normalMaxObsRate} \
--normalMaxObsNum=${normalMaxObsNum} \
--heteroSNPMinDepth=${heteroSNPMinDepth} \
--heteroSNPMinObsRate=${heteroSNPMinObsRate} \
--heteroSNPMinObsNum=${heteroSNPMinObsNum} \
--ohvar_mutGammaF=${ohvar_mutGammaF} \
--ohvar_mutGammaH=${ohvar_mutGammaH} \
--ohvar_mutAlphaL=${ohvar_mutAlphaL} \
--ohvar_mutAlphaH=${ohvar_mutAlphaH} \
--ohvar_mutGammaEH=${ohvar_mutGammaEH} \
--ohvar_mutAlphaS=${ohvar_mutAlphaS} \
--ohvar_mutAlphaB_E=${ohvar_mutAlphaB_E} \
--ohvar_mutAlphaB_W=${ohvar_mutAlphaB_W} \
--ohvar_errAlphaL_E=${ohvar_errAlphaL_E} \
--ohvar_errAlphaH_E=${ohvar_errAlphaH_E} \
--ohvar_errGammaEH_E=${ohvar_errGammaEH_E} \
--ohvar_errAlphaS_E=${ohvar_errAlphaS_E} \
--ohvar_errAlphaB_E=${ohvar_errAlphaB_E} \
--ohvar_errAlphaL_W=${ohvar_errAlphaL_W} \
--ohvar_errAlphaH_W=${ohvar_errAlphaH_W} \
--ohvar_errGammaEH_W=${ohvar_errGammaEH_W} \
--ohvar_errAlphaS_W=${ohvar_errAlphaS_W} \
--ohvar_errAlphaB_W=${ohvar_errAlphaB_W} \
--tumorMinAvgBaseQuality=${tumorMinAvgBaseQuality} \
--normalMinAvgBaseQuality=${normalMinAvgBaseQuality} \
--heteroSNPMinAvgBaseQuality=${heteroSNPMinAvgBaseQuality} \
--triAlleleMinObsRate=${triAlleleMinObsRate} \
--triAlleleMinObsNum=${triAlleleMinObsNum} \
--averageMapPhredQualThreshold=${averageMapPhredQualThreshold} \
--softClipPosessionFreqThreshold=${softClipPosessionFreqThreshold} \
--pileUpBufferSize=${pileUpBufferSize} \
${isSingle} \
-R ${REGION}"
${DIR}/../bin/ohvarfinder \
--algorithm OHVarfinDer \
-f ${REF}    \
-a ${TUMOR}  \
-b ${NORMAL} \
-o ${OUTPUTDIR}/output \
--maxInsertSize=${maxInsertSize} \
--tumorMinDepth=${tumorMinDepth} \
--tumorMinObsRate=${tumorMinObsRate} \
--tumorMaxObsRate=${tumorMaxObsRate} \
--tumorMinObsNum=${tumorMinObsNum} \
--normalMinDepth=${normalMinDepth} \
--normalMaxObsRate=${normalMaxObsRate} \
--normalMaxObsNum=${normalMaxObsNum} \
--heteroSNPMinDepth=${heteroSNPMinDepth} \
--heteroSNPMinObsRate=${heteroSNPMinObsRate} \
--heteroSNPMinObsNum=${heteroSNPMinObsNum} \
--ohvar_mutGammaF=${ohvar_mutGammaF} \
--ohvar_mutGammaH=${ohvar_mutGammaH} \
--ohvar_mutAlphaL=${ohvar_mutAlphaL} \
--ohvar_mutAlphaH=${ohvar_mutAlphaH} \
--ohvar_mutGammaEH=${ohvar_mutGammaEH} \
--ohvar_mutAlphaS=${ohvar_mutAlphaS} \
--ohvar_mutAlphaB_E=${ohvar_mutAlphaB_E} \
--ohvar_mutAlphaB_W=${ohvar_mutAlphaB_W} \
--ohvar_errAlphaL_E=${ohvar_errAlphaL_E} \
--ohvar_errAlphaH_E=${ohvar_errAlphaH_E} \
--ohvar_errGammaEH_E=${ohvar_errGammaEH_E} \
--ohvar_errAlphaS_E=${ohvar_errAlphaS_E} \
--ohvar_errAlphaB_E=${ohvar_errAlphaB_E} \
--ohvar_errAlphaL_W=${ohvar_errAlphaL_W} \
--ohvar_errAlphaH_W=${ohvar_errAlphaH_W} \
--ohvar_errGammaEH_W=${ohvar_errGammaEH_W} \
--ohvar_errAlphaS_W=${ohvar_errAlphaS_W} \
--ohvar_errAlphaB_W=${ohvar_errAlphaB_W} \
--tumorMinAvgBaseQuality=${tumorMinAvgBaseQuality} \
--normalMinAvgBaseQuality=${normalMinAvgBaseQuality} \
--heteroSNPMinAvgBaseQuality=${heteroSNPMinAvgBaseQuality} \
--triAlleleMinObsRate=${triAlleleMinObsRate} \
--triAlleleMinObsNum=${triAlleleMinObsNum} \
--averageMapPhredQualThreshold=${averageMapPhredQualThreshold} \
--softClipPosessionFreqThreshold=${softClipPosessionFreqThreshold} \
--pileUpBufferSize=${pileUpBufferSize} \
${isSingle} \
-R ${REGION}

echo "cat ${OUTPUTDIR}/output.calls.txt | grep -E -v ${filterExpression} | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.filt.variant"
cat ${OUTPUTDIR}/output.calls.txt | grep -E -v ${filterExpression} | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.filt.variant

echo "cat ${OUTPUTDIR}/output.calls.txt | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.variant"
cat ${OUTPUTDIR}/output.calls.txt | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.variant
