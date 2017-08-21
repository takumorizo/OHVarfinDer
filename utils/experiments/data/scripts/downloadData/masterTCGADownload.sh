#! /bin/sh
#$ -S /bin/sh
#$ -cwd

: <<'#__CO__'
# HapMuC real data download commands
bash ./masterTCGADownload.sh \
${outputDir} \
../../TCGA/TCGA_mutation_calling_benchmark_files/gdc_manifest_tcga_mutation_calling_benchmark_files.txt
#__CO__


readonly DIR=`dirname ${0}`

config=${DIR}/config.sh
utility=${DIR}/utility.sh
source ${config}
source ${utility}


outputDir=$1
manifest=$2

check_mkdir ${outputDir}

CURLOGDIR=${LOGDIR}
check_mkdir ${CURLOGDIR}
LOGSTR=-e\ ${CURLOGDIR}\ -o\ ${CURLOGDIR}

__sum_up_manifest_filesize() {
  local line total=0
  while read ID FILENAME MD5 FILESIZE STATE; do
	if [[ ${line} -ge 1 ]]; then
	    total=$(( total + FILESIZE ))
		export line=$((line + 1))
	elif [[ ${line} -eq 0  ]]; then
		export line=$((line + 1))
	fi
  done

  echo $total
  return 0
}

total=`cat $manifest | __sum_up_manifest_filesize `
echo "total=${total}"

TOTALSIZETB=$((total / 1000000000000))
TOTALSIZEGB=$((total / 1000000000))
TOTALSIZEMB=$((total / 1000000))

echo "${TOTALSIZETB} (Tb) ${TOTALSIZEGB} (Gb) ${TOTALSIZEMB} (Mb) ( ${total} (byte) ) files will be down loaded. Do you exec downloads? [y/n]"
read ANSWER

case $ANSWER in
    "" | "Y" | "y" | "yes" | "Yes" | "YES" ) echo "YES!! Start downloading.";;
    * )
		echo "NO!! Stop downloading."
		exit 0
		;;
esac

LINE=0
cat ${manifest} | while read ID FILENAME MD5 FILESIZE STATE ; do
	if [[ ${LINE} -ge 1 ]]; then
		echo "qsub -v CONFIG=${config} -v UTIL=${utility} -v DIR=${DIR} \
		-l s_vmem=4G,mem_req=4G ${LOGSTR} \
		${DIR}/eachTCGADownload.sh ${outputDir} ${FILENAME} ${BASEURL} ${ID} ${MD5}"
		qsub -v CONFIG=${config} -v UTIL=${utility} -v DIR=${DIR} \
		-l s_vmem=4G,mem_req=4G ${LOGSTR} \
		${DIR}/eachTCGADownload.sh ${outputDir} ${FILENAME} ${BASEURL} ${ID} ${MD5}

		LINE=$((LINE + 1))
	elif [[ ${LINE} -eq 0  ]]; then
		LINE=$((LINE + 1))
	fi
done

