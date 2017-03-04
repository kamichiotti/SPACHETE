#!/bin/bash -eu

#  linkfastafiles.sh
#  
#
#  Created by Gillian Hsieh on 1/31/16.
#

##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.


FJDir=${1} #MACHETE output directory
STEMFILE=${1}/StemList.txt
TASK_ID=${2}
SPORK_fasta=${3}
INSTALLDIR=${4}

STEM=`awk 'FNR == '${TASK_ID}' {print $1}' ${STEMFILE}`

BigFastaDir=${1}/fasta/
ChrFastaDir=${1}/fasta/${STEM}/
BigFastaFile=${1}/$SPORK_fasta
BowtieIndex=${1}/BowtieIndex/${STEM}/${STEM}_FJ_Index
mkdir -p ${1}/BowtieIndex/${STEM}/
CopyFastaName=${BigFastaDir}${STEM}_FarJunctions.fa
cp ${BigFastaFile} ${CopyFastaName}

## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>_FJ_Index
bowtie2-build ${CopyFastaName} ${BowtieIndex}

#Making individual chromosome fastas just to follow gillian's methods
echo "python ${INSTALLDIR}/split_large_fasta.py ${CopyFastaName}"
python ${INSTALLDIR}/split_large_fasta.py ${CopyFastaName}

echo "linkfastafiles and Bowtie Index built for Far Junctions for sample ${STEM}" >> ${1}/MasterError.txt
