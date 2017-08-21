#! /bin/sh
#$ -S /bin/sh
#$ -cwd

outputDir=$1
FILENAME=$2
BASEURL=$3
ID=$4
MD5=$5

echo "outputDir : ${outputDir}"
echo "FILENAME : ${FILENAME}"
echo "BASEURL : ${BASEURL}"
echo "ID : ${ID}"
echo "MD5 : ${MD5}"


source ${UTIL}

check_mkdir ${outputDir}
cd ${outputDir}

again=true

while ${again} ; do
	echo "wget -O ${outputDir}/${FILENAME} ${BASEURL}/${ID}"
	wget -O ${outputDir}/${FILENAME} ${BASEURL}/${ID}

	if [ -e ${outputDir}/${FILENAME} ]; then
		echo "${outputDir}/${FILENAME}"
	else
		echo "${outputDir}/${FILENAME} not exists"
		echo "${outputDir}/${FILENAME} not exists" >&2
		exit 1
	fi

	MD5RESULT=$( md5sum ${outputDir}/${FILENAME} | awk '{print $1}')
	echo "${MD5RESULT} ; result md5sum"
	echo "${MD5} ; expected md5sum"

	if [[ "${MD5}" = "${MD5RESULT}" ]]; then
		echo "MD5 sum is ok"
		echo "finish downloading file ${outputDir}/${FILENAME}"
		again=false
	else
		if [ -e ${outputDir}/${FILENAME} ]; then
			echo "rm ${outputDir}/${FILENAME}"
			rm ${outputDir}/${FILENAME}
		fi
		again=true
	fi
done



