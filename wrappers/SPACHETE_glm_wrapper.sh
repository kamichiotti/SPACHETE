#!/bin/bash

#This small wrapper script restarts a SPACHETE (or MACHETE) job
#starting near the very end for the GLM so we can more easily tinker w/
#the GLM and then get rapid results

#############################
#    Input and outputs      #
#############################
#KNIFE_DIR="" #<- This is the 'input' dir
#OUT_DIR=""   #<- This is the 'output' dir

#############################
#       Ewing samples       #
#############################
KNIFE_DIR="/scratch/PI/horence/gillian/Ewing/circpipe"
OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/ewing_9_1"

#############################
#    CML UConn samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/CML_uconn_8_31/"
#STEM_INCLUDE_ONLY_LIST=("SRR3192410")

#############################
#     CML test samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_test/aligned/CML"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/cml_test_8_25"

#################################
#    Normal breast samples      #
#################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_breast/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_breast_9_1_36-cutoff-glm"
#STEM_INCLUDE_ONLY_LIST=("SRR1027188")

######################################
#    Simulated engstrom samples      #
###############################3######
#KNIFE_DIR="/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/gillian/engstrom"

################################
#    Normal fetal samples      #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_fetal/circpipe_fetal"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_fetal_8_31"
#STEM_INCLUDE_ONLY_LIST=("Fetal_Heart_397_CAGATC_L004")

################################
#        Bladder samples       #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/bladder/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/bladder_8_27"

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
###                      Paths that are the same as MACHETE                        ###
###                                                                                ###
###                                                                                ###
######################################################################################

#Build up the OPTIONS string
OPTIONS=""
OPTIONS="$OPTIONS --circpipe-dir $KNIFE_DIR"
OPTIONS="$OPTIONS --output-dir $OUT_DIR"
if [ "$STEM_INCLUDE_ONLY_LIST" -a ${#STEM_INCLUDE_ONLY_LIST[@]} -gt 0 ]
then
    INCLUDE_LIST=$(printf ",%s" "${STEM_INCLUDE_ONLY_LIST[@]}")
    INCLUDE_LIST=${INCLUDE_LIST:1}
    OPTIONS="$OPTIONS --stem-include-list $INCLUDE_LIST"
fi

echo "spachete_glm_feeder.py $OPTIONS 1> out_spachete_feeder.txt 2> out_spachete_feeder.err"
python spachete_glm_feeder.py $OPTIONS 1> out_spachete_feeder.txt 2> out_spachete_feeder.err


