from argparse import ArgumentParser
import subprocess
import shutil
import glob
import sys
import pdb
import re
import os

############################
#Get self path to this file#
############################
MACHETE = os.path.dirname(os.path.abspath(__file__))

#######################
#  Parse Arguments    #
#######################
description = "Required args: --circpipe-dir, --output-dir, --hg19Exons, --reg-indel-indices, and --circref-dir."
parser = ArgumentParser(description=description)
parser.add_argument("--circpipe-dir",required=True,dest="CIRCPIPE_DIR",help="Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)")
parser.add_argument("--output-dir",required=True,dest="OUTPUT_DIR",help="Output directory for resuts. Directory path will be created recursively if it doesn't already exist. If it exists already, the directory will be deleted then created again.")
parser.add_argument("--stem-include-list",required=False,dest="STEM_INCLUDE_LIST",help="If only want to run some of the stems then define this, otherwise will run all the stems")


#################
#  Build Paths  #
#################
args = parser.parse_args()
CIRCPIPE_DIR = args.CIRCPIPE_DIR
OUTPUT_DIR = args.OUTPUT_DIR

ORIG_DIR = os.path.join(CIRCPIPE_DIR,"orig")
FullStemFile = os.path.join(OUTPUT_DIR,"StemList.txt")


##########################
#  Specify stems to run  #
##########################
include_stems = args.STEM_INCLUDE_LIST
if include_stems:
    include_stems = include_stems.split(",")


##########################################################
# Running and Building output directories for each Stem  #
##########################################################
processes = []
out_files = []
err_files = []
with open(FullStemFile,"r") as stem_file:
    stems = stem_file.readlines()
    for stem in sorted(stems):
        stem = stem.strip()
        print stem
        #If include_stems is not None then
        #skip all the stems not in the include stems list
        if include_stems and stem not in include_stems:
            print stem+" not in include stems, skipping "+stem
            continue

        OUTPUT_DIR = args.OUTPUT_DIR
        OUTPUT_DIR = os.path.join(OUTPUT_DIR,stem)


        OPTIONS = ""
        OPTIONS += " --circpipe-dir "+CIRCPIPE_DIR
        OPTIONS += " --output-dir "+OUTPUT_DIR


        machete_out = open(os.path.join(OUTPUT_DIR,"spachete_"+stem+".out"),"w")
        machete_err = open(os.path.join(OUTPUT_DIR,"spachete_"+stem+".err"),"w")
        out_files.append(machete_out)
        err_files.append(machete_err)

        #Sending jobs to the SLURM scheduler (comment out if not on SLURM)
        job_name = stem[-5:] #annoying that the full job name doesn't show
        slurm_out_name = os.path.join(OUTPUT_DIR,"slurm_glm_"+stem+".out")
        slurm_err_name = os.path.join(OUTPUT_DIR,"slurm_glm_"+stem+".err")

        sub_SLURM_name = os.path.join(OUTPUT_DIR,stem+"_job.sh")
        with open(sub_SLURM_name,"w") as sub_SLURM:
            sub_SLURM.write("#!/bin/bash\n")
            sub_SLURM.write("python "+os.path.join(MACHETE,"spachete_glm_run.py")+OPTIONS+"\n")

        subprocess.call(["sbatch","-p","horence","--mem=30000","--time=8:00:00",
                         "-o",slurm_out_name,"-e",slurm_err_name,"-J",job_name,sub_SLURM_name],
                         stdout=machete_out,stderr=machete_err)
        """

        #This is for multithreading w/out submitting a SLURM job (comment out if on SLURM)
        processes.append(subprocess.Popen(["python",os.path.join(MACHETE,"spachete_glm_run.py"),
                                           "--circpipe-dir",CIRCPIPE_DIR,
                                           "--output-dir",OUTPUT_DIR,
                                           "--hg19Exons",EXONS,
                                           "--reg-indel-indices",REG_INDEL_INDICES,
                                           "--circref-dir",CIRCREF],
                                           stdout=machete_out,stderr=machete_err))

        """

#Force all the jobs to complete before exiting
for p_ind in range(len(processes)):
    processes[p_ind].communicate()
    out_files[p_ind].close()
    err_files[p_ind].close()

sys.stdout.write("Finished spachete feeder\n")


