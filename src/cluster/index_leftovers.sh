#!/bin/bash

if [ $# -ne 8 ]; then
	echo "$0 mapping_files_list contig_breaks output_dir read_len mean_insert_size stdev_insert_size min_mps_per_cluster contig_name_file"
fi 

# The file containing a list of the mapping files, one file per line
RMAPPING_FILES=$1
# The directory where the clusters and link edges are stored
OUTPUT_DIR=$3
CONTIG_NAME_FILE=$8
CLUSTER_FILE=${OUTPUT_DIR}/clusters
LINK_FILE=${OUTPUT_DIR}/links



mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/indices
rm -rf ${OUTPUT_DIR}/indexmaps_script
rm -rf ${OUTPUT_DIR}/indexclus_script

MMAPPING_FILES=
for RMAPPINGS in `cat ${RMAPPING_FILES}`
do
	MMAPPINGS=${OUTPUT_DIR}/`basename ${RMAPPINGS}`.disc 
	MMAPPING_FILES="$MMAPPING_FILES ${MMAPPINGS}"
done

#build map indices
for RMAPPINGS in `cat ${RMAPPING_FILES}`
do
	while read CONTIG
	do
		MMAPPINGS=${OUTPUT_DIR}/`basename ${RMAPPINGS}`.disc.${CONTIG}
		INDEX=${OUTPUT_DIR}/indices/`basename ${RMAPPINGS}`.disc.${CONTIG}.idx
		echo  "${CNVER_FOLDER}/cluster/idx_build ${INDEX} ${MMAPPINGS} " >> ${OUTPUT_DIR}/indexmaps_script
	done < ${CONTIG_NAME_FILE}
done

#build cluster indices
while read CONTIG; do
	for TYPE in 0 1 2 3; do
		INDEX=${OUTPUT_DIR}/clusters.${CONTIG}.t${TYPE}.idx
		DATA=${CLUSTER_FILE}.$CONTIG.t$TYPE
		echo "${CNVER_FOLDER}/cluster/idx_build $INDEX $DATA delim HEADER1" >> ${OUTPUT_DIR}/indexclus_script
	done
done < ${CONTIG_NAME_FILE}
