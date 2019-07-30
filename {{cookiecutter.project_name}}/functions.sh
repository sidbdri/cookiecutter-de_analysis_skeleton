#!/bin/bash

STAR=STAR{{cookiecutter.star_version}}
FEATURE_COUNTS=featureCounts{{cookiecutter.featurecounts_version}}
PICARD=/opt/picard-tools-{{cookiecutter.picard_version}}/picard.jar

function listFiles {
    local DELIMITER=$1
    shift
    local FILES=$@

    local OUTPUT=$(ls -1 $FILES | tr '\n' "${DELIMITER}")
    echo ${OUTPUT%$DELIMITER}
}

function listFilesNoNewLine {
    local DELIMITER=$1
    shift
    local FILES=$@

    local OUTPUT=$(ls -1 $FILES | tr '\n' "${DELIMITER}")
    echo -n ${OUTPUT%$DELIMITER}
}

function map_reads {
    local SAMPLE=$1
    local INDEX_DIR=$2
    local NUM_THREADS=$3
    local READ_1_FILES=$4
    local READ_2_FILES=$5
    local MAPPING_DIR=$6


    local star_tmp=${SAMPLE}.tmp
    mkdir ${star_tmp}

    if [ -z "$READ_2_FILES" ]
    then
        #if read_2 files is empty, the input is single end read
        read_files_opt="--readFilesIn ${READ_1_FILES}"
    else
        read_files_opt="--readFilesIn ${READ_1_FILES} ${READ_2_FILES}"
    fi

    ${STAR} --runThreadN ${NUM_THREADS} --genomeDir ${INDEX_DIR} --outFileNamePrefix ${star_tmp}/star --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate Unsorted --readFilesCommand zcat ${read_files_opt}

    mv ${star_tmp}/starAligned.out.bam ${MAPPING_DIR}/${SAMPLE}.bam
    mv ${star_tmp}/starAligned.sortedByCoord.out.bam ${MAPPING_DIR}/${SAMPLE}.sorted.bam
    mv ${star_tmp}/starLog.final.out ${MAPPING_DIR}/${SAMPLE}.log.out

    rm -rf ${star_tmp}
}

function picard_rnaseq_metrics {
    local SAMPLE=$1
    local INPUT_BAM=$2
    local OUTPUT_DIR=$3
    local REF_FLAT=$4
    local RIBOSOMAL_DIR=$5
    local STRANDNESS_FLAG=${6:-2}

    local S=''

    case ${STRANDNESS_FLAG} in
    "0") S=NONE ;;
    "1") S=FIRST_READ_TRANSCRIPTION_STRAND ;;
    "2") S=SECOND_READ_TRANSCRIPTION_STRAND ;;
    *) echo "Unrecognized strandness: ${STRANDNESS_FLAG}. Can only be {0,1,2}"; exit 1 ;;
    esac

    sambamba view -H ${INPUT_BAM} > ${RIBOSOMAL_DIR}/${SAMPLE}_header.txt
    cat ${RIBOSOMAL_DIR}/${SAMPLE}_header.txt ${RIBOSOMAL_DIR}/intervalListBody.txt > ${RIBOSOMAL_DIR}/${SAMPLE}.txt

    java -jar ${PICARD} CollectRnaSeqMetrics I=${INPUT_BAM} O=${OUTPUT_DIR}/${SAMPLE}.txt REF_FLAT=${REF_FLAT} STRAND=${S} RIBOSOMAL_INTERVALS=${RIBOSOMAL_DIR}/${SAMPLE}.txt
}

function count_reads_for_features {
    local NUM_THREADS=$1
    local FEATURES_GTF=$2
    local BAM_FILE=$3
    local COUNTS_OUTPUT_FILE=$4
    local STRANDNESS_FLAG=${5:-2}


    local counts_tmp=.$(basename "${BAM_FILE}").counts_tmp

    ${FEATURE_COUNTS} -T ${NUM_THREADS} -p -a ${FEATURES_GTF} -o ${counts_tmp} -s ${STRANDNESS_FLAG} ${BAM_FILE}
    tail -n +3 ${counts_tmp} | cut -f 1,7 > ${COUNTS_OUTPUT_FILE}

    rm ${counts_tmp}
    mv ${counts_tmp}.summary ${COUNTS_OUTPUT_FILE}.summary
}


function detect_stranness {
    ## If 1 is close to 2, then it is 0
    ## Otherwise, it is the larger one among 1 and 2.
    local DIR=$1
    local SAMPLE_NAME=${2:-''}

    # we need to find out the correct species file for this sample, that is whichever has
    # the largest amount of assigned reads from the [0]s. See https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/85
    local species=`grep Assigned ${DIR}/*${SAMPLE_NAME}*[0].testsummary* | \
    awk '{print $2"\t"$1}' | sort -nr | head -1 | cut -f 2 | awk -F"." '{print $4}'`

    local zero="$(grep Assigned ${DIR}/*${SAMPLE_NAME}.${species}.counts.0.testsummary* | awk '{print $2}' ) "
    local one="$(grep Assigned ${DIR}/*${SAMPLE_NAME}.${species}.counts.1.testsummary* | awk '{print $2}'  )"
    local two="$(grep Assigned ${DIR}/*${SAMPLE_NAME}.${species}.counts.2.testsummary* | awk '{print $2}'  )"
#    echo -n ${SAMPLE_NAME}: ${zero} ${one} ${two} ""
#    echo ${one} ${two} | awk '{print ($1-$2)/($1+$2)}'

    echo ${one} ${two} | awk 'function abs(v) {return v < 0 ? -v : v}
                                { d=($1-$2)/($1+$2)
                                if( abs(d)<0.75) {print 0}
                                else if(d>=0.75) {print 1}
                                else if(d<=-0.75) {print 2}
                                else {print -1}
                                }'
    }


function count_reads_for_features_strand_test {
    local NUM_THREADS=$1
    local FEATURES_GTF=$2
    local BAM_FILE=$3
    local COUNTS_OUTPUT_FILE=$4

    local counts_tmp=.$(basename "${BAM_FILE}").counts_tmp


    for i in 0 1 2; do
        ${FEATURE_COUNTS} -T ${NUM_THREADS} -p -a ${FEATURES_GTF} -o ${counts_tmp} -s $i ${BAM_FILE}
        tail -n +3 ${counts_tmp} | cut -f 1,7 > ${COUNTS_OUTPUT_FILE}.$i

        rm ${counts_tmp}
        mv ${counts_tmp}.summary ${COUNTS_OUTPUT_FILE}.$i.testsummary
    done
}

function clean_de_results {
    local DE_RESULTS_FILE=$1

    sed -i "s/-Inf/'-Inf/g;s/,NA,/,,/g;s/,NA,/,,/g;s/,NA$/,/g;s/,NaN,/,,/g;s/,NaN,/,,/g;s/,NaN$/,/g" ${DE_RESULTS_FILE}
}

function addSample2tsv {
    local SAMPLE_TSV=$1 && shift
    local BASE_DIR=$1 && shift
    local READ1_IDENTIFIER=$1 && shift
    local READ2_IDENTIFIER=$1 && shift
    local FASTQ_SUFFIX=$1 && shift
    local PAIRED_READ=$1 && shift
    if [ -z "$READ2_IDENTIFIER" ];then
        local PAIRED_READ=0
    fi
    echo -ne '' > ${SAMPLE_TSV}
    local SAMPLE=$@
    for sample in ${SAMPLE}; do
        echo -ne ${sample}" " >> ${SAMPLE_TSV}
        echo -n $(listFilesNoNewLine "," ${BASE_DIR}/${sample}/*${READ1_IDENTIFIER}.${FASTQ_SUFFIX}) >> ${SAMPLE_TSV}
        if [ "$PAIRED_READ"=="yes" ]; then
            echo -n " "  >> ${SAMPLE_TSV}
            echo -n $(listFilesNoNewLine "," ${BASE_DIR}/${sample}/*${READ2_IDENTIFIER}.${FASTQ_SUFFIX}) >> ${SAMPLE_TSV}
        fi
        echo "" >> ${SAMPLE_TSV}
    done
}
