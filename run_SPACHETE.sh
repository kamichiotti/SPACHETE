#!/bin/bash
#SPACHETE call wrapper
#
#First please set the EXONS_DIR, INDEL_INDICES, and CIRC_REF which will be the same as MACHETE
#
# To run spachete you need to specify:
#(1) The path to the KNIFE named directory of interest
#(2) The path to where you want to save the output.
#    will be built if it doesn't exist
#
#Can optionally provide a list of keys to specify which files in the circpipe
#directory get run. Only samples that contain at least one of the keys get run.
#If this is not specified that every file in the circpipe directory will get run
#(3) This is the STEM_INCLUDE_LIST=("key") optional variable

#####################################
#    Example Input and outputs      #
#####################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/peter/A011_and_A012/"               #<-- Input path to the KNIFE output
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/peter_ovarian"                  #<-- Output path, will get made if does not exist
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("A012-CTTGTA_S3_L006") #<-- Optional, if not set all files will be run
#Choices for KNIFE dir
#A002-CGATGT_S1_L006/ A011_and_A012/       A016_and_A018/       HA08_and_HA09/
#A003_and_A007/       A014_and_A015/       A019_and_HA07/

#####################################
#             Fetal Test            #
#####################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/fetal_lung/fetal_lung"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/fetal_lung_test"
#MODE="hg19"

#####################################
#          CML Test                 #
#####################################
KNIFE_DIR="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562"
OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CML_uconn_test_3_1_17"
MODE="hg19"
STEM_INCLUDE_ONLY_LIST=("SRR3192411")


#####################################
#       References Paths            #
# (should be the same across runs)  #
#####################################
#NOTE!!
#Have to change these paths to to same ones that MACHETE points to
EXONS_DIR="/scratch/PI/horence/gillian/HG19exons"
INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/toyIndelIndices/"
CIRC_REF="/share/PI/horence/circularRNApipeline_Cluster/index"


#####################################
#      Dont have to change          #
#####################################
#Get the absolute path to this wrapper script
ABS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #<-- very cool line gets loc of this dir
TIME_STAMP=`date +"%Y_%m_%d-%H_%M"`
#OUT_DIR=$OUT_DIR"_"$TIME_STAMP #<-- uncomment to avoid write-over

#Build up the OPTIONS string
OPTIONS=""
OPTIONS="$OPTIONS --root-dir $ABS_PATH"
OPTIONS="$OPTIONS --circpipe-dir $KNIFE_DIR"
OPTIONS="$OPTIONS --mode $MODE"
OPTIONS="$OPTIONS --output-dir $OUT_DIR"
OPTIONS="$OPTIONS --hg19Exons $EXONS_DIR"
OPTIONS="$OPTIONS --reg-indel-indices $INDEL_INDICES"
OPTIONS="$OPTIONS --circref-dir $CIRC_REF"
OPTIONS="$OPTIONS --SLURM True" #<-- comment this line out if not on SLURM
if [ "$STEM_INCLUDE_ONLY_LIST" -a ${#STEM_INCLUDE_ONLY_LIST[@]} -gt 0 ]
then
    INCLUDE_LIST=$(printf ",%s" "${STEM_INCLUDE_ONLY_LIST[@]}")
    INCLUDE_LIST=${INCLUDE_LIST:1}
    OPTIONS="$OPTIONS --stem-include-list $INCLUDE_LIST"
fi

LOG_STEM="$ABS_PATH/logs/spachete_feeder_$TIME_STAMP" 
python "$ABS_PATH/wrappers/spachete_feeder.py" $OPTIONS 1> "$LOG_STEM.out" 2> "$LOG_STEM.err"


