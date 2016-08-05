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
#    EWS-Fli  Test run      #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/EWS_FLI_bigmem/EWS_bigmem"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/EWS_FLI_test"

#############################
#    NSC-Diff  Test run     #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/NSC_diff/NSC_diff_test"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/NSC_diff_test"

#############################
#    CML UConn samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/CML_filtered"

#################################
#    Normal breast samples      #
#################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_breast/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_breast"

######################################
#    Simulated engstrom samples      #
###############################3######
#KNIFE_DIR="/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/engstrom"

################################
#    Normal fetal samples      #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_fetal/circpipe_fetal"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_fetal"

################################
#  Salzman Ovarian samples     #
################################
#KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/OvarianCancer2014_cutAdapt"
#OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/salzman_ov"

##############################
#  Normal Genome Samples     #
##############################
KNIFE_DIR="/scratch/PI/horence/rob/parent_dirs/normal_human_genome/ERR1549500"
OUT_DIR="/scratch/PI/horence/rob/spachete_outputs/normal_genome"









############################################
#PATHS THAT DONT NEED TO CHANGE ON SHERLOCK#
############################################
#INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices"
EXONS_DIR="/scratch/PI/horence/gillian/HG19exons"
INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/toyIndelIndices/"
CIRC_REF="/share/PI/horence/circularRNApipeline_Cluster/index"


#Build up the OPTIONS string
OPTIONS=""
OPTIONS="$OPTIONS --circpipe-dir $KNIFE_DIR"
OPTIONS="$OPTIONS --output-dir $OUT_DIR"
OPTIONS="$OPTIONS --hg19Exons $EXONS_DIR"
OPTIONS="$OPTIONS --reg-indel-indices $INDEL_INDICES"
OPTIONS="$OPTIONS --circref-dir $CIRC_REF"

python spachete_feeder.py $OPTIONS 1> o.txt 2> o.err


