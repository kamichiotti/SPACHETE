#!/bin/bash
#MACHETE call scheme with arguments
#description = "Required args: --circpipe-dir, --output-dir, --hg19Exons, --reg-indel-indices, and --circref-dir."

#Optional arguments and description:
#parser.add_argument("--user-bp-dist",dest="USERBPDIST",type=int,default=1000,help="Default is %(default)s.")
#parser.add_argument("--num-junc-bases",dest="NUMBASESAROUNDJUNC",type=int,default=8,help="Default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70.")
#parser.add_argument("--numIndels",dest="NumIndels",type=int,default=5,help="Default is %(default)s.")


#Example call: (pretending everything is on the same line)
#run.py --circpipe-dir /scratch/PI/horence/alignments/EWS_FLI_bigmem/
#       --output-dir /scratch/PI/horence/alignments/EWS_FLI_bigmem/FarJunc/
#       --hg19Exons /scratch/PI/horence/gillian/HG19exons/
#       --reg-indel-indices /scratch/PI/horence/gillian/HG19_reg_indels/
#       --circref-dir /share/PI/horence/circularRNApipeline_Cluster/index/ 


#############################
#    Input and outputs      #
#############################
#KNIFE_DIR=""
#OUT_DIR=""

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
OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/CML_new_gtf_test/"
STEM_INCLUDE_ONLY_LIST=("SRR3192409")

#############################
#     CML test samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_test/aligned/CML"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/cml_test"

#################################
#    Normal breast samples      #
#################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_breast/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_breast_8_15"

######################################
#    Simulated engstrom samples      #
###############################3######
#KNIFE_DIR="/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/engstrom"

################################
#    Normal fetal samples      #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_fetal/circpipe_fetal"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_fetal_8_15"
#STEM_INCLUDE_ONLY_LIST=("Fetal_Adrenal_360_CTTGTA_L006"
#                        "Fetal_Adrenal_403b_GTCCGC_L008"
#                        "Fetal_Adrenal_408_CGATGT_L008"
#                        "Fetal_Adrenal_419_TGACCA_L008"
#                        "Fetal_Brain_408_AGTCAA_L006"
#                        "Fetal_Heart_361_TGACCA_L007"
#                        "Fetal_Heart_365_ACAGTG_L007"
#                        "Fetal_Intestine_395_CGATGT_L005"
#                        "Fetal_Lung_384_CTTGTA_L007"
#                        "Fetal_Stomach_360_ACAGTG_L008"
#                        "Fetal_Stomach_401_GTCCGC_L006")

################################
#        Bladder samples       #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/bladder/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/bladder"

#############################
#       Ewing samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/Ewing/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/ewing"

#############################
#      RNaseR samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/ov_RNaseR_Qatar/ovCircPipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/RNaseR"

#############################
#    Prostate samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/prostate/aligned/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/prostate"

#############################
#        SEQC study         #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/SEQC_study_set/circpipe_SEQC"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/seqc"




######################################################################################
###                                                                                ###
###                                                                                ###
###                              Personal Samples                                  ###
###                                                                                ###
###                                                                                ###
######################################################################################

################################
#  Salzman Ovarian samples     #
################################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/OvarianCancer2014_cutAdapt"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/salzman_ov"

#############################
#    EWS-Fli  Test run      #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/EWS_FLI_bigmem/EWS_bigmem"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/EWS_FLI_test"

#############################
#    NSC-Diff  Test run     #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/NSC_diff/NSC_diff_test"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/NSC_diff_test"

##############################
#  Normal Genome Samples     #
##############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/normal_human_genome/ERR1549500"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_genome_collapsed"



######################################################################################
###                                                                                ###
###                                                                                ###
###                      Paths that are the same as MACHETE                        ###
###                                                                                ###
###                                                                                ###
######################################################################################

EXONS_DIR="/scratch/PI/horence/gillian/HG19exons"
INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/toyIndelIndices/"
CIRC_REF="/share/PI/horence/circularRNApipeline_Cluster/index"
#EXONS_DIR=""
#INDEL_INDICES=""
#CIRC_ReF=""

#Build up the OPTIONS string
OPTIONS=""
OPTIONS="$OPTIONS --circpipe-dir $KNIFE_DIR"
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

python spachete_feeder.py $OPTIONS 1> out_spachete_feeder.txt 2> out_spachete_feeder.err


