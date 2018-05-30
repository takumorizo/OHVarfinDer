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
--algorithm OHVarfinDer2 \
-f ${REF}    \
-a ${TUMOR}  \
-b ${NORMAL} \
-o ${OUTPUTDIR} \
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
--ohvar2_mutGammaF=${ohvar2_mutGammaF} \
--ohvar2_mutGammaH=${ohvar2_mutGammaH} \
--ohvar2_mutAlphaL=${ohvar2_mutAlphaL} \
--ohvar2_mutAlphaH=${ohvar2_mutAlphaH} \
--ohvar2_mutGammaEH=${ohvar2_mutGammaEH} \
--ohvar2_mutAlphaS=${ohvar2_mutAlphaS} \
--ohvar2_mutAlphaB_E=${ohvar2_mutAlphaB_E} \
--ohvar2_mutAlphaB_W=${ohvar2_mutAlphaB_W} \
--ohvar2_errAlphaL_E=${ohvar2_errAlphaL_E} \
--ohvar2_errAlphaH_E=${ohvar2_errAlphaH_E} \
--ohvar2_errGammaEH_E=${ohvar2_errGammaEH_E} \
--ohvar2_errAlphaS_E=${ohvar2_errAlphaS_E} \
--ohvar2_errAlphaB_E=${ohvar2_errAlphaB_E} \
--ohvar2_errAlphaL_W=${ohvar2_errAlphaL_W} \
--ohvar2_errAlphaH_W=${ohvar2_errAlphaH_W} \
--ohvar2_errGammaEH_W=${ohvar2_errGammaEH_W} \
--ohvar2_errAlphaS_W=${ohvar2_errAlphaS_W} \
--ohvar2_errAlphaB_W=${ohvar2_errAlphaB_W} \
--tumorMinAvgBaseQuality=${tumorMinAvgBaseQuality} \
--normalMinAvgBaseQuality=${normalMinAvgBaseQuality} \
--heteroSNPMinAvgBaseQuality=${heteroSNPMinAvgBaseQuality} \
--triAlleleMinObsRate=${triAlleleMinObsRate} \
--triAlleleMinObsNum=${triAlleleMinObsNum} \
--pileUpBufferSize=${pileUpBufferSize} \
${isSingle} \
-R ${REGION}"
${DIR}/../bin/ohvarfinder \
--algorithm OHVarfinDer2 \
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
--ohvar2_mutGammaF=${ohvar2_mutGammaF} \
--ohvar2_mutGammaH=${ohvar2_mutGammaH} \
--ohvar2_mutAlphaL=${ohvar2_mutAlphaL} \
--ohvar2_mutAlphaH=${ohvar2_mutAlphaH} \
--ohvar2_mutGammaEH=${ohvar2_mutGammaEH} \
--ohvar2_mutAlphaS=${ohvar2_mutAlphaS} \
--ohvar2_mutAlphaB_E=${ohvar2_mutAlphaB_E} \
--ohvar2_mutAlphaB_W=${ohvar2_mutAlphaB_W} \
--ohvar2_errAlphaL_E=${ohvar2_errAlphaL_E} \
--ohvar2_errAlphaH_E=${ohvar2_errAlphaH_E} \
--ohvar2_errGammaEH_E=${ohvar2_errGammaEH_E} \
--ohvar2_errAlphaS_E=${ohvar2_errAlphaS_E} \
--ohvar2_errAlphaB_E=${ohvar2_errAlphaB_E} \
--ohvar2_errAlphaL_W=${ohvar2_errAlphaL_W} \
--ohvar2_errAlphaH_W=${ohvar2_errAlphaH_W} \
--ohvar2_errGammaEH_W=${ohvar2_errGammaEH_W} \
--ohvar2_errAlphaS_W=${ohvar2_errAlphaS_W} \
--ohvar2_errAlphaB_W=${ohvar2_errAlphaB_W} \
--tumorMinAvgBaseQuality=${tumorMinAvgBaseQuality} \
--normalMinAvgBaseQuality=${normalMinAvgBaseQuality} \
--heteroSNPMinAvgBaseQuality=${heteroSNPMinAvgBaseQuality} \
--triAlleleMinObsRate=${triAlleleMinObsRate} \
--triAlleleMinObsNum=${triAlleleMinObsNum} \
--pileUpBufferSize=${pileUpBufferSize} \
${isSingle} \
-R ${REGION}

echo "cat ${OUTPUTDIR}/output.calls.txt | grep -E -v ${filterExpression} | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.filt.variant"
cat ${OUTPUTDIR}/output.calls.txt | grep -E -v ${filterExpression} | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.filt.variant

echo "cat ${OUTPUTDIR}/output.calls.txt | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.variant"
cat ${OUTPUTDIR}/output.calls.txt | python ${DIR}/rmSNP.py > ${OUTPUTDIR}/output.variant
