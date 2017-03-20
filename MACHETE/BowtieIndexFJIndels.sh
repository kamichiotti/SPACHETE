#!/bin/bash -eu

#  BowtieIndexFJIndels.sh
#  
#
#  Created by Gillian Hsieh on 2/9/16.
#
## The shell script BowtieIndexFJIndels.sh calls bowtie to index the indels_N.fa files that were created in the previous step.  The indels are output to the directory FJDir/BowtieIndels/<STEM>/<STEM>_Indels_N where N is the number of indels in that index.

FastaDir=${1}
IndelNumber=${2}
BowtieIndexDir=${3}
FJDir=${4}
TASK_ID=${5}

STEMFILE=${4}/StemList.txt
STEM=`awk 'FNR == '${TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${3}/${STEM}

bowtie2-build ${1}/${STEM}/${STEM}_FJ_Indels_${2}.fa ${3}/${STEM}/${STEM}_Indels_${2}

echo "BowtieIndexFJIndels.sh complete. Check for the bowtie index ${4}/BowtieIndels/${STEM}/${STEM}_Indels_${2}" >> ${4}/MasterError.txt
