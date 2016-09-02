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


