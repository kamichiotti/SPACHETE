#!/bin/bash
#SPACHETE call wrapper
#
#To run spachete you need to specify:
#(1) The path to the KNIFE named directory of interest
#(2) The path to where you want to save the output.
#    will be built if it doesn't exist
#
#Can optionally provide a list of keys to specify which files in the circpipe
#directory get run. Only samples that contain at least one of the keys get run.
#If this is not specified that every file in the circpipe directory will get run
#(3) This is the STEM_INCLUDE_LIST=("key") optional variable

#############################
#    Input and outputs      #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/peter/A011_and_A012/"               #<-- Input path to the KNIFE output
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/peter_ovarian"                  #<-- Output path, will get made if does not exist
#STEM_INCLUDE_ONLY_LIST=("A012-CTTGTA_S3_L006") #<-- Optional, if not set all files will be run
#Choices for KNIFE dir
#A002-CGATGT_S1_L006/ A011_and_A012/       A016_and_A018/       HA08_and_HA09/
#A003_and_A007/       A014_and_A015/       A019_and_HA07/

######################################################################################
###                                                                                ###
###                                                                                ###
###                           Gilian Samples Machete                               ###
###                                                                                ###
###                                                                                ###
######################################################################################

#############################
#    CML UConn samples      #
#############################
KNIFE_DIR="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562"
OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CML_uconn_test"
MODE="hg19"
STEM_INCLUDE_ONLY_LIST=("SRR3192409")

#############################
#     CML test samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_test/aligned/CML"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/cml_test_9_15"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("ENCFF000HOC2")

#################################
#    Normal breast samples      #
#################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_breast/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/normal_breast_9_9"
#STEM_INCLUDE_ONLY_LIST=("SRR1027188")
#MODE="hg19"

######################################
#    Simulated engstrom samples      #
###############################3######
#KNIFE_DIR="/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/engstrom_9_9"
#MODE="hg19"

################################
#    Normal fetal samples      #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_fetal/circpipe_fetal"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/normal_fetal_500K_11_19"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("Fetal_Lung_361_CAGATC_L005"  "Fetal_Stomach_361_GCCAAT_L008" "Fetal_Stomach_403b_AGTTCC_L006")

################################
#        Bladder samples       #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/bladder/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/bladder_8_27"
#MODE="hg19"

#############################
#       Ewing samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/Ewing/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/ewing_9_14"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("SRR1594020"  "SRR1594021"  "SRR1594022"  "SRR1594023"  "SRR1594024"  "SRR1594025")
#STEM_INCLUDE_ONLY_LIST=("SRR1594020" "SRR1594021")

#############################
#      RNaseR samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/ov_RNaseR_Qatar/ovCircPipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/ovcar3_9_15"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("SRR1772257" "SRR1772957" "SRR1777309" "SRR1777310")

#############################
#    Prostate samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/prostate/aligned/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/gillian/prostate"
#MODE="hg19"

#############################
#        SEQC study         #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/SEQC_study_set/circpipe_SEQC"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/gillian/seqc"
#MODE="hg19"




######################################################################################
###                                                                                ###
###                                                                                ###
###                              Personal Samples                                  ###
###                                                                                ###
###                                                                                ###
######################################################################################

################################
#    Quake Mouse samples       #
################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/quake/quake5"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/quake_12_5"
#MODE="mm10"
#STEM_INCLUDE_ONLY_LIST=("SRR3137883")

################################
#  Salzman Ovarian samples     #
################################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/OvarianCancer2014_cutAdapt"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/gillian/salzman_ov"
#MODE="hg19"

#############################
#    EWS-Fli  Test run      #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/EWS_FLI_bigmem/EWS_bigmem"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/EWS_FLI_test"
#MODE="hg19"

#############################
#    NSC-Diff  Test run     #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/NSC_diff/NSC_diff_test"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/NSC_diff_test"
#MODE="hg19"

##############################
#  Normal Genome Samples     #
##############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/normal_human_genome/ERR1549500"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/normal_genome_collapsed"
#MODE="hg19"

##############################################################
#  Just running on cell free DNA from a transplant study:    #
#            http://www.ncbi.nlm.nih.gov/pubmed/24267896     #
##############################################################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/cfDNA_SRR1024346/single_end"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/cfDNA_SRR1024346_8_25"
#MODE="hg19"


######################################################################################
###                                                                                ###
###                                                                                ###
###                      Paths that are the same as MACHETE                        ###
###                                                                                ###
###                                                                                ###
######################################################################################
# Don't have to change any of this code regardless of sample
EXONS_DIR="/scratch/PI/horence/gillian/HG19exons"
INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/toyIndelIndices/"
CIRC_REF="/share/PI/horence/circularRNApipeline_Cluster/index"
#EXONS_DIR=""
#INDEL_INDICES=""
#CIRC_ReF=""

#Build up the OPTIONS string
OPTIONS=""
OPTIONS="$OPTIONS --circpipe-dir $KNIFE_DIR"
OPTIONS="$OPTIONS --mode $MODE"
OPTIONS="$OPTIONS --output-dir $OUT_DIR"
OPTIONS="$OPTIONS --hg19Exons $EXONS_DIR"
OPTIONS="$OPTIONS --reg-indel-indices $INDEL_INDICES"
OPTIONS="$OPTIONS --circref-dir $CIRC_REF"
if [ "$STEM_INCLUDE_ONLY_LIST" -a ${#STEM_INCLUDE_ONLY_LIST[@]} -gt 0 ]
then
    INCLUDE_LIST=$(printf ",%s" "${STEM_INCLUDE_ONLY_LIST[@]}")
    INCLUDE_LIST=${INCLUDE_LIST:1}
    OPTIONS="$OPTIONS --stem-include-list $INCLUDE_LIST"
fi

python /scratch/PI/horence/rob/SPACHETE_dirs/SPACHETE/spachete_feeder.py $OPTIONS 1> out_spachete_wrapper.txt 2> out_spachete_wrapper.err


