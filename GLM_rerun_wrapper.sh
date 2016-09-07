
########################################
#              Main paths              #
########################################
#Likely WILL want to change these
OUT_PATH="/scratch/PI/horence/rob/SPACHETE/glm_reruns/"
GLM_PATH="/scratch/PI/horence/rob/SPACHETE/GLM_script_UseIndel.r"

########################################
#      Function to call GLM rerun      #
########################################
glm_rerun(){

    echo "In glm_rerun fnct: CIRC_PIPE = $1"
    echo "In glm_rerun fnct: FJDIR     = $2"
    echo "In glm_rerun fnct: OUT_DIR   = $3"
    echo "In glm_rerun fnct: GLM_PATH  = $4"

    cmd="/scratch/PI/horence/rob/SPACHETE/All_Samples_GLMandReports.sh "
    cmd+="$1 "
    cmd+="$2 "
    cmd+="$3 "
    cmd+="$4"

    #echo "sh $cmd"
    sh $cmd
}

########################################
#            Ewing samples             #
########################################
#CIRC_PIPE="/scratch/PI/horence/gillian/Ewing/circpipe/"
#FJDIR_STEM="/scratch/PI/horence/rob/spachete_outputs/ewing_9_1/"
#OUT_DIR="ewing_9_1/"
#STEMS=("SRR1594020")

########################################
#           Uconn CML samples          #
########################################
CIRC_PIPE="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562/"
FJDIR_STEM="/scratch/PI/horence/rob/spachete_outputs/CML_uconn_8_31/"
OUT_DIR="cml_uconn_8_31/"
STEMS=("SRR3192409" "SRR3192410" "SRR3192411")


########################################
#  Loop through calling the GLM rerun  #
########################################
#Shouldn't have to change this
for STEM in ${STEMS[@]}
do
    FJDIR="$FJDIR_STEM""$STEM/"
    OUT="$OUT_PATH""$OUT_DIR""$STEM/""reports/"
    glm_rerun $CIRC_PIPE $FJDIR $OUT $GLM_PATH
done




