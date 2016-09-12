#!/bin/bash

function listFiles {
    DELIMITER=$1
    shift
    FILES=$@

    OUTPUT=$(ls -1 $FILES | tr '\n' "${DELIMITER}")
    echo ${OUTPUT%$DELIMITER}
}

function map_reads {
    SAMPLE=$1
    INDEX_DIR=$2
    NUM_THREADS=$3
    READ_1_FILES=$4
    READ_2_FILES=$5
    OUTPUT_DIR=$6

    star_tmp=${SAMPLE}.tmp
    mkdir ${star_tmp}

    STAR --runThreadN ${NUM_THREADS} --genomeDir ${INDEX_DIR} --readFilesIn ${READ_1_FILES} ${READ_2_FILES} --outFileNamePrefix ${star_tmp}/star --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

    mv ${star_tmp}/starAligned.sortedByCoord.out.bam ${OUTPUT_DIR}/${SAMPLE}.bam
    mv ${star_tmp}/starLog.final.out ${OUTPUT_DIR}/${SAMPLE}.log.out

    rm -rf ${star_tmp}
}

function count_reads_for_features {
    NUM_THREADS=$1
    FEATURES_GTF=$2
    BAM_FILE=$3
    COUNTS_OUTPUT_FILE=$4

    counts_tmp=.counts_tmp

    featureCounts -T ${NUM_THREADS} -p -a ${FEATURES_GTF} -o ${counts_tmp} -s 2 ${BAM_FILE}
    tail -n +3 ${counts_tmp} | cut -f 1,7 > ${COUNTS_OUTPUT_FILE}

    rm ${counts_tmp}
    mv ${counts_tmp}.summary ${COUNTS_OUTPUT_FILE}.summary
}
