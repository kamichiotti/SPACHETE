#!/bin/bash -eu

## The main wrapper script to generate reports containing read count and statistical score for each circular or linear junction
# detected in the data. This works on our specific deployment of a SLURM job scheduler and is intended as a template to modify to
# work with the specific syntax required by your scheduler. Here, the task array feature of SLURM is used to run tasks in parallel
# on individual files where possible, waiting to proceed until all files successfully processed at a specific step when necessary.
# The main steps are:
#  1) align to fastq files to genome, transcriptome, linear junction, circular junction, and ribosomal indices using Bowtie2 
#  2) parse sam files and analyze alignments using "naive" statistical method
#  3) run GLM to report posterior probability and p-value per junction
#  4) generate alignment statistics per sample (useful for evaluating efficiency of ribo depletion, global circular:linear ratio, etc)
#     and identify unaligned reads 
#  5) identify candidate de novo junctions from unaligned reads and generate a Bowtie2 de novo index from these candidates
#  6) align unaligned reads to de novo index and analyze using "naive" statistical method (requires 2nd call to findCircularRNA.sh using "unaligned" mode)

# Usage (see README for descriptions of parameters):
#    sh findCircularRNA.sh read_directory read_id_style alignment_parent_directory dataset_name junction_overlap [mode] [report_directory_name] [ntrim] [denovoCircMode] [junction_id_suffix] 
#    And then run again with same parameters as before but append "_unaligned" to the mode parameter

cd /share/PI/horence/circularRNApipeline_SLURM/

source analysis/depends_python.sh  # use correct version of python

# very basic error checking
if [ $# -eq 11 ]
then
  echo "Read 1 threshold was specified without Read 2 threshold. Please include Read 2 threshold or check your arguments if you did not intend to specify Read 1 threshold ${11}"
  exit 2
fi

CODE_DIR=`pwd`
READ_DIR=${1}
READ_STYLE=${2}
ALIGN_PARDIR=${3}
DATASET_NAME=${4}
OVERLAP=${5}
if [ $# -ge 6 ]
then
  MODE=${6}
else
  MODE=sam
fi

if [ $# -ge 7 ]
then
  REPORTDIR_NAME=${7}
else
  REPORTDIR_NAME=circReads
fi

# ntrim for denovo
if [ $# -ge 8 ]
then
  NTRIM=${8}
else
  NTRIM=50
fi

# should denovo contains only circles (1 means only circles)
if [ $# -ge 9 ]
then
  DENOVOCIRC=${9}
else
  DENOVOCIRC=1
fi

JUNCTION_DIR_SUFFIX=${10}
RD1_THRESH=${11}
RD2_THRESH=${12}

if [ $# -ge 13 ]
then
  JUNCTION_MIDPOINT=${13}
else
  JUNCTION_MIDPOINT=150
fi

# change resource allocations based on job size
if [[ "$MODE" = *large* ]]
then
  JUNC_VMEM="80000"
  GENOME_VMEM="80000"
  TRANSC_VMEM="80000"
  REG_VMEM="80000"
  RIBO_VMEM="80000"
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_VMEM="200000"
  FILTER_VMEM="200000"
  PREPROCESS_MAX_RT="23:0:0"
  FILTER_MAX_RT="23:0:0"
  PFA_VMEM="200000"
elif [[ "$MODE" = *bigmem* ]]
then
  JUNC_VMEM="80000"
  GENOME_VMEM="80000"
  TRANSC_VMEM="80000"
  REG_VMEM="80000"
  RIBO_VMEM="80000"
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_VMEM="80000"
  FILTER_VMEM="200000"  # this is used for the naive reports
  PREPROCESS_MAX_RT="16:0:0"  # this is a limit put on bigmem
  FILTER_MAX_RT="16:0:0"  # this is a limit put on bigmem
  PFA_VMEM="200000"  # this is used for glm reports
else
  JUNC_VMEM="20000"
  GENOME_VMEM="6000"
  TRANSC_VMEM="6000"
  REG_VMEM="20000"
  RIBO_VMEM="5000"
  ALIGN_MAX_RT="12:0:00"
  PREPROCESS_VMEM="15000"
  FILTER_VMEM="15000"
  PREPROCESS_MAX_RT="4:0:0"
  FILTER_MAX_RT="12:0:0"
  PFA_VMEM="40000"
fi

if [[ "$MODE" = *bam* ]]
then
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_MAX_RT="23:0:0"
  FILTER_MAX_RT="23:0:0"
fi

# we can ask specifically to run on our own machine
if [[ "$MODE" = *horence* ]]
then
  RESOURCE_FLAG="-p horence"
elif [[ "$MODE" = *owners* ]]
then
  RESOURCE_FLAG="-p owners"
fi

# sometimes we need more memory for analysis part only
if [[ "$MODE" = *bigmem* ]]
then
  ANALYSIS_RESOURCE_FLAG="-p bigmem --qos=bigmem"
else
  ANALYSIS_RESOURCE_FLAG=$RESOURCE_FLAG
fi


# a few special flags for the unaligned mode
if [[ $MODE = *unaligned* ]]
then
  TASK_FILE_READ_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned
  UFLAG="-u"
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt
else
  TASK_FILE_READ_DIR=${READ_DIR}
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
fi  

# set up info about this dataset, create directory structures
if [[ $MODE != *analysis* ]]
then
  # set up
  python analysis/writeTaskIdFiles.py -r ${TASK_FILE_READ_DIR} -a ${ALIGN_PARDIR} -d ${DATASET_NAME} ${UFLAG}
  
  # select correct prefix name to use for bowtie index files
  
  if [[ $MODE = *grch38* ]]
  then
    bt_prefix="grch38"
  elif [[ $MODE = *mouse* ]]
  then
    bt_prefix="mm10"
  elif [[ $MODE = *rat* ]]
  then
    bt_prefix="rn5"
  elif [[ $MODE = *elegans* ]]
  then
    bt_prefix="celegans"
  elif [[ $MODE = *fly* ]]
  then
    bt_prefix="dm3"
  elif [[ $MODE = *pombe* ]]
  then
    bt_prefix="ASM294v2_23"
  elif [[ $MODE = *crypto* ]]
  then
    bt_prefix="cryptococcus_neoformans_grubii_h99"
  elif [[ $MODE = *cerevisiae* ]]
  then
    bt_prefix="Scer"
  elif [[ $MODE = *mikatae* ]]
  then
    bt_prefix="Smik"
  elif [[ $MODE = *bayanus* ]]
  then
    bt_prefix="Sban"
  elif [[ $MODE = *HSV* ]]
  then
    bt_prefix="KOS"
  elif [[ $MODE = *capsas* ]]
  then
    bt_prefix="capsaspora_atcc_30864_2"
  elif [[ $MODE = *rosetta* ]]
  then
    bt_prefix="salpingoeca_rosetta_1"
  else
    bt_prefix="hg19"
  fi
fi

# break out of alignment-only loop since this is needed for both alignment and analysis
NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

# create directory and subdirectories for log output
LOG_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/logs
mkdir -p ${LOG_DIR}
mkdir -p ${LOG_DIR}/align
mkdir -p ${LOG_DIR}/analysis
mkdir -p ${LOG_DIR}/glm
mkdir -p ${LOG_DIR}/sampleStats
mkdir -p ${LOG_DIR}/denovo_script_out
mkdir -p ${LOG_DIR}/denovo_index

# alignments
if [[ $MODE != *analysis* ]]
then
  if [[ $MODE = *unaligned* ]]
  then
    j_id=`sbatch -J DeNovoAlign${DATASET_NAME}${DENOVOCIRC} ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${ALIGN_MAX_RT} --mem=6000 -D ${CODE_DIR}/index -o ${LOG_DIR}/align/${DATASET_NAME}AlignDeNovo${DENOVOCIRC}_%A_%a.out -e ${LOG_DIR}/align/${DATASET_NAME}AlignDeNovo${DENOVOCIRC}_%A_%a.err analysis/align.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE denovo_${DATASET_NAME}_${DENOVOCIRC} denovo denovo_${DATASET_NAME}_onlycircles${DENOVOCIRC}.fa | awk '{print $4}'`
    depend_str="--depend=afterok:"${j_id}
  else
    g_id=`sbatch -J GenomeAlign${DATASET_NAME} ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${ALIGN_MAX_RT} --mem=${GENOME_VMEM} -D ${CODE_DIR}/index -o ${LOG_DIR}/align/${DATASET_NAME}AlignGenome_%A_%a.out -e ${LOG_DIR}/align/${DATASET_NAME}AlignGenome_%A_%a.err analysis/align.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_genome genome ${bt_prefix}_genome.fa | awk '{print $4}'`
    j_id=`sbatch -J JunctionAlign${DATASET_NAME} ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${ALIGN_MAX_RT} --mem=${JUNC_VMEM} -D ${CODE_DIR}/index -o ${LOG_DIR}/align/${DATASET_NAME}AlignJunction_%A_%a.out -e ${LOG_DIR}/align/${DATASET_NAME}AlignJunction_%A_%a.err analysis/align.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_junctions_scrambled junction ${bt_prefix}_junctions_scrambled.fa | awk '{print $4}'`
    r_id=`sbatch -J RiboAlign${DATASET_NAME} ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${ALIGN_MAX_RT} --mem=${RIBO_VMEM} -D ${CODE_DIR}/index -o ${LOG_DIR}/align/${DATASET_NAME}AlignRibo_%A_%a.out -e ${LOG_DIR}/align/${DATASET_NAME}AlignRibo_%A_%a.err analysis/align.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_ribosomal ribo ${bt_prefix}_ribosomal.fa | awk '{print $4}'`
    t_id=`sbatch -J TranscAlign${DATASET_NAME} ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${ALIGN_MAX_RT} --mem=${TRANSC_VMEM} -D ${CODE_DIR}/index -o ${LOG_DIR}/align/${DATASET_NAME}AlignTranscriptome_%A_%a.out -e ${LOG_DIR}/align/${DATASET_NAME}AlignTranscriptome_%A_%a.err analysis/align.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_transcriptome transcriptome ${bt_prefix}_transcriptome.fa | awk '{print $4}'`
    reg_id=`sbatch -J RegAlign${DATASET_NAME} ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${ALIGN_MAX_RT} --mem=${REG_VMEM} -D ${CODE_DIR}/index -o ${LOG_DIR}/align/${DATASET_NAME}AlignReg_%A_%a.out -e ${LOG_DIR}/align/${DATASET_NAME}AlignReg_%A_%a.err analysis/align.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $MODE ${bt_prefix}_junctions_reg reg ${bt_prefix}_junctions_reg.fa | awk '{print $4}'`  
    depend_str="--depend=afterok:${g_id}:${j_id}:${r_id}:${t_id}:${reg_id}"
  fi
fi



# preprocessing
p_id=`sbatch -J ${DATASET_NAME}Preprocess ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${PREPROCESS_MAX_RT} --mem=${PREPROCESS_VMEM} -D ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/${DATASET_NAME}Preprocess_%A_%a.out -e ${LOG_DIR}/analysis/${DATASET_NAME}Preprocess_%A_%a.err ${depend_str} analysis/preprocessAlignedReads.sh ${TASK_DATA_FILE} ${ALIGN_PARDIR} ${DATASET_NAME} ${JUNCTION_MIDPOINT} ${OVERLAP} ${JUNCTION_DIR_SUFFIX} | awk '{print $4}'`

# analysis
a_id=`sbatch -J ${DATASET_NAME}Analysis ${ANALYSIS_RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=${FILTER_VMEM} --time=${FILTER_MAX_RT} -D ${CODE_DIR}/analysis -o ${LOG_DIR}/analysis/${DATASET_NAME}Analysis_%A_%a.out -e ${LOG_DIR}/analysis/${DATASET_NAME}Analysis_%A_%a.err --depend=afterok:${p_id} analysis/filterFDR.sh ${MODE} ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${READ_STYLE} ${OVERLAP} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} | awk '{print $4}'`

# quality stats, glm, and denovo only if this is not unaligned mode 
if [[ ${MODE} != *unaligned* ]]
then
  
  if [[ ${MODE} != *skipGLM* ]]
  then
    TEMP_OUT_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/orig/temp # to store the output files with data parsed from SAM
    mkdir -p ${TEMP_OUT_DIR}
    pfa_id=`sbatch -J ${DATASET_NAME}ParseForAnalysis ${RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${PREPROCESS_MAX_RT} --mem=${PREPROCESS_VMEM} -D ${CODE_DIR}/analysis -o ${LOG_DIR}/glm/${DATASET_NAME}ParseForAnalysis_%A_%a.out -e ${LOG_DIR}/glm/${DATASET_NAME}ParseForAnalysis_%A_%a.err --depend=afterok:${a_id} analysis/parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} | awk '{print $4}'`
    sbatch -J PFA${DATASET_NAME} ${ANALYSIS_RESOURCE_FLAG} --array=1-${NUM_FILES} --time=${PREPROCESS_MAX_RT} --mem=${PFA_VMEM} -D ${CODE_DIR}/analysis -o ${LOG_DIR}/glm/${DATASET_NAME}PredictJunctions_%A_%a.out -e ${LOG_DIR}/glm/${DATASET_NAME}PredictJunctions_%A_%a.err --depend=afterok:${pfa_id} analysis/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}
  fi
  
  # make the output directories
  OUT_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/sampleStats
  mkdir -p ${OUT_DIR}
  mkdir -p ${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned
  mkdir -p ${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned/forDenovoIndex
    
  # for each sample, generate sample stats in the output directory
  q_id=`sbatch -J Quality${DATASET_NAME} ${RESOURCE_FLAG} -a 1-${NUM_FILES} --time=48:0:0 --mem=44000 -D ${CODE_DIR}/qualityStats -o ${LOG_DIR}/sampleStats/${DATASET_NAME}Quality_%A_%a.out -e ${LOG_DIR}/sampleStats/${DATASET_NAME}Quality_%A_%a.err --depend=afterok:${a_id} qualityStats/qualityStatsSingleSample.sh $TASK_DATA_FILE $ALIGN_PARDIR $DATASET_NAME $REPORTDIR_NAME ${OUT_DIR} ${NTRIM} ${JUNCTION_DIR_SUFFIX} | awk '{print $4}'`
    
  # cat all of those files then delete the original separate files
  sbatch -J Cat${DATASET_NAME} ${RESOURCE_FLAG} --mem=5000 --time=1:0:0 -D ${CODE_DIR}/qualityStats -o ${LOG_DIR}/sampleStats/${DATASET_NAME}Cat_%A_%a.out -e ${LOG_DIR}/sampleStats/${DATASET_NAME}Cat_%A_%a.err --depend=afterok:${q_id} qualityStats/qualityStatsCat.sh ${OUT_DIR}

  if [[ ${MODE} != *skipDenovo* ]]    # create the index for this dataset
  then
    sbatch -J Denovo_${DENOVOCIRC}_${DATASET_NAME} ${RESOURCE_FLAG} --mem=55000 --time=23:0:0 -D ${CODE_DIR}/denovo_scripts -o ${LOG_DIR}/denovo_script_out/${DATASET_NAME}_${DENOVOCIRC}DenovoScript.out -e ${LOG_DIR}/denovo_script_out/${DATASET_NAME}_${DENOVOCIRC}DenovoScript.err --depend=afterok:${q_id} denovo_scripts/denovo_scripts.batch ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${NTRIM} ${DENOVOCIRC}
  fi
fi

  # delete files
  if [[ "$MODE" = *delete* ]]
  then
    sbatch -J ${TRIMMED_DATASET_NAME}DelFiles ${RESOURCE_FLAG} -o ${ALIGN_PARDIR}/${TRIMMED_DATASET_NAME}/${TRIMMED_DATASET_N
AME}DelDirs.out -e ${ALIGN_PARDIR}/${TRIMMED_DATASET_NAME}/${TRIMMED_DATASET_NAME}DelDirs.err --depend=afterany:${ua_id} dele
teDirs.sh ${RAW_READ_DIR} ${TRIMMED_READ_DIR} ${ALIGN_PARDIR}/${TRIMMED_DATASET_NAME}
  fi








