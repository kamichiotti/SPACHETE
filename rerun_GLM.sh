#!/bin/sh

#  rerun_GLM.sh
#  
#
#  Created by Gillian Hsieh on 3/3/16.
#
GLM_PATH=${1}
CIRC_PIPE=${2} #directory that contains circReads, orig, logs, etc
FJDIR=${3}
INSTALLDIR=${4}

if [ $# -ge 5 ]
then
OUTPUTDIR=${5}
else
OUTPUTDIR=${FJDIR}reports/glmReports/
fi
echo $OUTPUTDIR

STEMFILE=${FJDIR}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
#STEM=SRR1594021

REG_INPUTDIR=${CIRC_PIPE}circReads/ids/
FJ_INPUTDIR=${FJDIR}GLM_classInput/

mkdir -p ${OUTPUTDIR}

for file in ${REG_INPUTDIR}*${STEM}*output.txt
do
reg_class_input=${file}
done

for file in ${FJ_INPUTDIR}*${STEM}*output_FJ.txt
do
FJ_input=${file}
done

for file in ${REG_INPUTDIR}*${STEM}*output_RegIndel.txt
do
RegIndel_input=${file}
done

for file in ${FJ_INPUTDIR}*${STEM}*output_FJIndels.txt
do
FJIndel_input=${file}
done

ml load R/3.0.2
echo "Rscript ${GLM_PATH} ${FJ_input} ${reg_class_input} ${STEM} ${OUTPUTDIR} ${RegIndel_input} ${FJIndel_input}"
Rscript ${GLM_PATH} ${FJ_input} ${reg_class_input} ${STEM} ${OUTPUTDIR} ${RegIndel_input} ${FJIndel_input}

echo "rerun_GLM.sh complete for ${STEM} -- check ${OUTPUTDIR}/*${STEM}* for Far Junction GLM outputs"
