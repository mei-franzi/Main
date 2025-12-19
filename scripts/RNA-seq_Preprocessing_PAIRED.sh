#!/bin/bash

### Sequencing R-ODAF, Omics Data Analysis Framework for Regulatory application ###

# Load environment
# conda activate rnaseq # Ensure this is active before running or uncomment if supported

project="rat-aging"
OUTPUT_DIR="/data/fm346/ngs/${project}/data/output"
RAW_SAMPLE_DIR="/data/fm346/ngs/${project}/data/fastq/"
SUFFIX_INPUTFILES='.fastq.gz'
SEQTYPE="RNASeq"
SEQMODE="paired"
PAIRED_END_SUFFIX_FORWARD="_1"
PAIRED_END_SUFFIX_REVERSE="_2"

ORGANISM_GENOME_ID="mRatBN7.2"
GENOME_FILES_DIR="/data/fm346/ngs/rat-aging/databases/rat_STAR_index"
GENOME_FILE_NAME="Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
GTF_FILE_NAME="Rattus_norvegicus.mRatBN7.2.110.gtf"

GENOME_INDEX_ALREADY_AVAILABLE="No" 
RSEM_INDEX_ALREADY_AVAILABLE="No"
LARGE_GENOME="Yes" # Set to Yes for Rat (mammalian genomes are 'Large' for STAR)

# System parameters
CPU_FOR_ALIGNMENT=10
CPU_FOR_OTHER=6

#######################################################################
### Variables Setup
#######################################################################

declare BASEDIR=${OUTPUT_DIR}
declare SUFFIX_IN=${SUFFIX_INPUTFILES}
declare PAIR1=${PAIRED_END_SUFFIX_FORWARD}
declare PAIR2=${PAIRED_END_SUFFIX_REVERSE}

declare GENOMEDIR=${GENOME_FILES_DIR}
declare GENOME="${GENOMEDIR}/${GENOME_FILE_NAME}"
declare GTF="${GENOMEDIR}/${GTF_FILE_NAME}"
declare GenomeID=${ORGANISM_GENOME_ID}

declare TRIMM_DIR="${BASEDIR}/Trimmed_reads/"
declare QC_DIR_fastp="${TRIMM_DIR}/fastpQCoutput/"
declare QC_DIR_multiQC="${TRIMM_DIR}/MultiQC/"
declare align_DIR="${TRIMM_DIR}/STAR/"
declare Quant_DIR="${TRIMM_DIR}/RSEM/"
declare RSEM_GENOMEDIR="${GENOMEDIR}/RSEM/"

declare SUFFIX_out="_trimmed${SUFFIX_IN}"

# Create directories
mkdir -p ${TRIMM_DIR} ${QC_DIR_fastp} ${QC_DIR_multiQC} ${align_DIR} ${Quant_DIR} ${RSEM_GENOMEDIR}

# Logging
mydate="$(date +'%d.%m.%Y.%H-%M')"
exec > >(tee -a "${OUTPUT_DIR}/log_${mydate}.out") 2>&1

echo "--- Starting Preprocessing for Project: ${project} ---"

##################################
### Trimming raw reads : Fastp ###
##################################

if [ ${SEQMODE} == "paired" ]; then
    echo "Checking for raw files in ${RAW_SAMPLE_DIR}..."
    for READ1 in ${RAW_SAMPLE_DIR}*${PAIR1}${SUFFIX_IN}; do
        [ -e "$READ1" ] || continue
        FILENAME=$(basename "$READ1")
        READ2="${RAW_SAMPLE_DIR}${FILENAME/${PAIR1}${SUFFIX_IN}/${PAIR2}${SUFFIX_IN}}"
        SAMPLE_ID="${FILENAME%${PAIR1}${SUFFIX_IN}}"
        
        OUT1="${TRIMM_DIR}${SAMPLE_ID}${PAIR1}${SUFFIX_out}"
        OUT2="${TRIMM_DIR}${SAMPLE_ID}${PAIR2}${SUFFIX_out}"

        if [ -f "$OUT1" ] && [ -f "$OUT2" ]; then
            echo "[SKIP] Trimmed files for ${SAMPLE_ID} already exist."
        else
            echo "[TRIMMING] fastp: ${SAMPLE_ID}"
            fastp --in1 "${READ1}" --in2 "${READ2}" --out1 "$OUT1" --out2 "$OUT2" \
                --json "${QC_DIR_fastp}${SAMPLE_ID}PE_fastp.json" \
                --html "${QC_DIR_fastp}${SAMPLE_ID}PE_fastp.html" \
                --cut_front --cut_tail --cut_right --length_required 36 --thread ${CPU_FOR_OTHER}
        fi
    done
fi

####################################################
### Alignment : STAR ###
####################################################



if [ ${SEQTYPE} == "RNASeq" ]; then
    # Indexing
    if [ "${GENOME_INDEX_ALREADY_AVAILABLE}" == "No" ]; then
        echo "[INDEXING] STAR Genome Indexing starting..."
        STAR --runMode genomeGenerate --genomeDir "${GENOMEDIR}" --genomeFastaFiles "${GENOME}" \
             --sjdbGTFfile "${GTF}" --sjdbOverhang 99 --runThreadN ${CPU_FOR_OTHER} \
             $( [ "$LARGE_GENOME" == "Yes" ] && echo "--genomeSAsparseD 1 --genomeChrBinNbits 15" )
    fi

    # Alignment
    for R1_PATH in ${TRIMM_DIR}*${PAIR1}${SUFFIX_out}; do
        [ -e "$R1_PATH" ] || continue
        FILENAME=$(basename "$R1_PATH")
        SAMPLE_ID="${FILENAME%${PAIR1}${SUFFIX_out}}"
        R2_PATH="${TRIMM_DIR}${SAMPLE_ID}${PAIR2}${SUFFIX_out}"
        
        if [ -f "${align_DIR}${SAMPLE_ID}Log.final.out" ]; then
            echo "[SKIP] Alignment for ${SAMPLE_ID} already exists."
        else
            echo "[ALIGNING] STAR: ${SAMPLE_ID}"
            STAR --runThreadN ${CPU_FOR_ALIGNMENT} --genomeDir "${GENOMEDIR}" \
                --readFilesIn "${R1_PATH}" "${R2_PATH}" \
                --quantMode TranscriptomeSAM --readFilesCommand zcat \
                --outFileNamePrefix "${align_DIR}${SAMPLE_ID}"
        fi
    done
fi

#######################
### Quantification : RSEM ###
#######################

# Indexing RSEM
if [ "${RSEM_INDEX_ALREADY_AVAILABLE}" == "No" ]; then
    echo "[INDEXING] RSEM Reference starting..."
    rsem-prepare-reference --gtf "${GTF}" "${GENOME}" "${RSEM_GENOMEDIR}${GenomeID}"
fi

# Quantifying
for R1_PATH in ${TRIMM_DIR}*${PAIR1}${SUFFIX_out}; do
    [ -e "$R1_PATH" ] || continue
    FILENAME=$(basename "$R1_PATH")
    SAMPLE_ID="${FILENAME%${PAIR1}${SUFFIX_out}}"
    
    if [ -f "${Quant_DIR}${SAMPLE_ID}.genes.results" ]; then
        echo "[SKIP] RSEM for ${SAMPLE_ID} already exists."
    else
        echo "[QUANTIFYING] RSEM: ${SAMPLE_ID}"
        rsem-calculate-expression -p ${CPU_FOR_OTHER} --paired-end \
            --bam "${align_DIR}${SAMPLE_ID}Aligned.toTranscriptome.out.bam" \
            --no-bam-output "${RSEM_GENOMEDIR}${GenomeID}" "${Quant_DIR}${SAMPLE_ID}"
    fi
done

# Summarize Results
echo "Generating Data Matrix..."
cd "${Quant_DIR}"
rsem-generate-data-matrix *.genes.results > genes.data.tsv
rsem-generate-data-matrix *.isoforms.results > isoforms.data.tsv

#######################
### Quality Control ###
#######################

echo "Running MultiQC..."
multiqc --cl-config "extra_fn_clean_exts: { '_fastp.json' }" "${BASEDIR}" \
        --filename MultiQC_Report.html --outdir "${QC_DIR_multiQC}"

echo "Preprocessing complete."
