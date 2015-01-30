#!/bin/bash

if [[ $# != 1 ]]
then
    echo "Usage: getReferenceInput.sh <output_dir>"
    exit 1
fi

echo "*************************"
echo "   GET REFERENCE INPUT"
echo "*************************"

FTP_ENSEMBLE=ftp://ftp.ensembl.org/pub
REL=release-78
SEQ_FORMAT=fasta
REL_NAME=GRCh38
ORGANISM=homo_sapiens
SEQ_TYPE=dna
FILE_PREFIX_NAME=Homo_sapiens.GRCh38.dna.chromosome
echo "Get input DNA sequences of human chromosomes in FASTA format:"
echo " - server: ${FTP_ENSEMBLE}"
echo " - release: ${REL} - (${REL_NAME})"
echo " - file format: ${SEQ_FORMAT}"
echo " - organism: ${ORGANISM}"
echo " - seq. type: ${SEQ_TYPE}"
echo " - file prefix names: ${FILE_PREFIX_NAME}"

URL=${FTP_ENSEMBLE}/${REL}/${SEQ_FORMAT}/${ORGANISM}/${SEQ_TYPE}
OUT_DIR=${1}/${REL_NAME}_${SEQ_FORMAT}_chr

mkdir -p ${OUT_DIR}

rm -f ${OUT_DIR}/CHECKSUMS
wget -P ${OUT_DIR} ${URL}/CHECKSUMS

for chr in {1..22} X Y
do
	echo "Downloading chr${chr}"
	FILE=${FILE_PREFIX_NAME}.${chr}.fa.gz
	SUM=`sum ${OUT_DIR}/${FILE}`
	if [ `grep -c "${SUM}" ${OUT_DIR}/CHECKSUMS` == 1 ]
	then
		echo "File ${OUT_DIR}/${FILE} already downloaded"
	else
		rm -f ${OUT_DIR}/${FILE}
		wget -P ${OUT_DIR} --limit-rate 3M ${URL}/${FILE}
	fi
done
echo "Download complete."


SEQ_FORMAT=gtf
SEQ_TYPE=genes
FILE_PREFIX_NAME=Homo_sapiens.GRCh38.78
echo "Get input annotations of human genes in GFF format:"
echo " - server: ${FTP_ENSEMBLE}"
echo " - release: ${REL} - (${REL_NAME})"
echo " - file format: ${SEQ_FORMAT}"
echo " - organism: ${ORGANISM}"
echo " - seq. type: ${SEQ_TYPE}"
echo " - file prefix names: ${FILE_PREFIX_NAME}"

URL=${FTP_ENSEMBLE}/${REL}/${SEQ_FORMAT}/${ORGANISM}
OUT_DIR=${1}/${REL_NAME}_${SEQ_FORMAT}_genes

mkdir -p ${OUT_DIR}

rm -f ${OUT_DIR}/CHECKSUMS
wget -P ${OUT_DIR} ${URL}/CHECKSUMS

echo "Downloading gene annotations"
FILE=${FILE_PREFIX_NAME}.gtf.gz
SUM=`sum ${OUT_DIR}/${FILE}`
if [ `grep -c "${SUM}" ${OUT_DIR}/CHECKSUMS` == 1 ]
then
	echo "File ${OUT_DIR}/${FILE} already downloaded"
else
	rm -f ${OUT_DIR}/${FILE}
	wget -P ${OUT_DIR} ${URL}/${FILE}
fi
echo "Download complete."

echo "Started parsing gene annotation file"
zcat ${OUT_DIR}/${FILE} | \
awk '\
BEGIN {FS="\t"} ($3 == "gene") && ($1 ~ /^[1-9XY]/) && ($2 ~ /havana/) \
	{split($9,f, " "); \
	gsub(";", "", f[2]); \
	gsub(";", "", f[4]); \
	gsub(";", "", f[6]); \
	print $1"\t"$4"\t"$5"\t"$7"\t"f[1]"="f[2]"\t"f[3]"="f[4]"\t"f[5]"="f[6] \
}' > ${OUT_DIR}/${FILE_PREFIX_NAME}.csv
echo "Finished parsing gene annotation file"
