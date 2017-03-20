#!/bin/sh

#  All_Samples_GLMandReports.sh
#  
#
#  Created by Gillian Hsieh on 4/23/16.
#
## This shell will loop through the following data sets on Sherlock to regenerate the GLM reports and appended naive reports
## 1. Ewing Sarcoma KD
## 2. Normal Breast
## 3. bladder cancer
## 4. CML_test (CML from caltech)
## 5. prostate cancer
## 6. normal_fetal
## 7. Engstrom simulated data set
## 8. SEQC - normal human and brain sets
## 9. Ov_RNaseR_Qatar
## 10. ovarian data from salzman lab !!!!! -RErun FJ indels class input?
## 11. CML_UConn

## This is the output directory for ALL above samples.  Directories for each sample will be created.
CIRC_PIPE=${1}
FJDIR=${2}
OUTPUT_DIR=${3}
GLM_PATH=${4}

if [ $# -ge 5 ]
then
MODE=${5}
else
MODE="horence"
fi
if [[ "$MODE" = *owners* ]]
then
RESOURCE_FLAG="-p owners"
fi
if [[ "$MODE" = *horence* ]]
then
RESOURCE_FLAG="-p horence"
fi


INSTALLDIR="/scratch/PI/horence/rob/SPACHETE/"

mkdir -p ${OUTPUT_DIR}



#######EWING SARCOMA KD###################
#CIRC_PIPE=/scratch/PI/horence/gillian/Ewing/circpipe/
#FJDIR=/scratch/PI/horence/gillian/Ewing/FarJunc/

Output=${OUTPUT_DIR}
mkdir -p ${Output}
mkdir -p ${Output}err/
mkdir -p ${Output}glmReports/
mkdir -p ${Output}AppendedReports/

echo "In All-rerun: CIRC_PIPE=$CIRC_PIPE"
echo "In All-rerun: FJDIR=$FJDIR"
echo "In All-rerun: OUTPUT_DIR=$OUTPUT_DIR"
echo "In All-rerun: GLM_PATH=$GLM_PATH"
NUM_FILES=$((`more ${FJDIR}StemList.txt | wc -l`))

#Build up the rerunGLM command
glm_cmd="sbatch -J GLM.r ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --time=5:0:0"
glm_cmd+=" -o ${Output}err/out_GLM_r.txt"
glm_cmd+=" -e ${Output}err/err_GLM_r.txt"
glm_cmd+=" ${INSTALLDIR}rerun_GLM.sh"
glm_cmd+=" ${GLM_PATH}"
glm_cmd+=" ${CIRC_PIPE}"
glm_cmd+=" ${FJDIR}"
glm_cmd+=" ${INSTALLDIR}"
glm_cmd+=" ${Output}glmReports/"

#Call the rerunGLM command
j1_id=`$glm_cmd | awk '{print $4}'`
#j1_id=`sbatch -J GLM.r ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --time=5:0:0 -o ${Output}err/out_GLM_r.txt -e ${Output}err/err_GLM_r.txt ${INSTALLDIR}rerun_GLM.sh ${GLM_PATH} ${CIRC_PIPE} ${FJDIR} ${INSTALLDIR} ${Output}glmReports/ | awk '{print $4}'`


#Work on building up the AppendedReports cmd
app_cmd="sbatch -J AppendNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --time=2:0:0"
app_cmd+=" -o ${Output}err/out_AppendNaive.txt"
app_cmd+=" -e ${Output}err/err_AppendNaive.txt"
app_cmd+=" --depend=afterok:${j1_id}"
app_cmd+=" ${INSTALLDIR}rerun_AppendNaiveRept.sh"
app_cmd+=" ${FJDIR}"
app_cmd+=" ${CIRC_PIPE}circReads/glmReports/"
app_cmd+=" ${INSTALLDIR}"
app_cmd+=" ${Output}glmReports/"
app_cmd+=" ${Output}AppendedReports/"

#Call the appendedReports cmd
j2_id=`$app_cmd | awk '{print $4}'`
#j2_id=`sbatch -J AppendNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --time=2:0:0 -o ${Output}err/out_AppendNaive.txt -e ${Output}err/err_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}rerun_AppendNaiveRept.sh ${FJDIR} ${CIRC_PIPE}circReads/glmReports/ ${INSTALLDIR} ${Output}glmReports/ ${Output}AppendedReports/ | awk '{print $4}'`



