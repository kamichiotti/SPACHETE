from argparse import ArgumentParser
import subprocess
import shutil
import glob
import sys
import pdb
import re
import os

# NOTE THAT AS OF JULY 2016, THIS VERSION ASSUMES THAT THERE IS ONLY ON
#   PAIR OF FASTQ FILES. SMALL CHANGES COULD FIX THIS PRESUMABLY
#   BY ADDING LOOPS OVER FILES TO 
#   THE CHANGES (MARKED BELOW) MADE IN JULY 2016.

#MACHETE = os.path.dirname(__file__)
MACHETE = os.path.dirname(os.path.abspath(__file__))
print "path:",MACHETE
# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for

def checkProcesses(popenDict):
    """
    Function :
    Args     : popenDict - A dict. whose keys are subprocess.Popen instances. The value of a key is also a dict. that has the following
                                keys: 'stdout', 'stderr', and 'cmd'.
    """
    for popen in popenDict:
        popen.communicate() #hangs until job finishes
    for popen in popenDict:
        retcode = popen.returncode
        cmd = popenDict[popen]['cmd']
        stdout = popenDict[popen]['stdout']
        stdout.close()
        stdout_text = open(stdout.name,'r').read()
        stderr = popenDict[popen]['stderr']
        stderr.close()
        stderr_text = open(stderr.name,'r').read()
        if retcode:
            raise Exception("Command '{cmd}' failed with return code {retcode}. Log files are {stdout} and {stderr}. {stdout} contains {stdout_text}. {stderr} contains {stderr_text}.".format(cmd=cmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name,stdout_text=stdout_text,stderr_text=stderr_text))
            # raise Exception("Command '{cmd}' failed with return code {retcode}. Log files are {stdout} and {stderr}.".format(cmd=cmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

description = "Required args: --circpipe-dir, --output-dir, --hg19Exons, --reg-indel-indices, and --circref-dir."

parser = ArgumentParser(description=description)
parser.add_argument("--circpipe-dir",required=True,dest="CIRCPIPE_DIR",help="Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)")
parser.add_argument("--output-dir",required=True,dest="OUTPUT_DIR",help="Output directory for resuts. Directory path will be created recursively if it doesn't already exist. If it exists already, the directory will be deleted then created again.")
parser.add_argument("--user-bp-dist",dest="USERBPDIST",type=int,default=1000,help="Default is %(default)s.")
#parser.add_argument("REFGENOME",help="HG19 vs HG38;could upgrade to HG38.")
parser.add_argument("--num-junc-bases",dest="NUMBASESAROUNDJUNC",type=int,default=8,help="Default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70.")
parser.add_argument("--numIndels",dest="NumIndels",type=int,default=5,help="Default is %(default)s.")
parser.add_argument("--hg19Exons",required=True,dest="EXONS",help="Path to HG19Exons. Formerly called PICKLEDIR.")
parser.add_argument("--reg-indel-indices",required=True,dest="REG_INDEL_INDICES",help="Path to files with names like hg19_junctions_reg_indels_1.1.bt2l,hg19_junctions_reg_indels_2.rev.1.bt2l ... These are sometimes in a folder called IndelIndices.")
parser.add_argument("--circref-dir",required=True,dest="CIRCREF",help="Path to reference libraries output by KNIFE - directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices.")

args = parser.parse_args()

CIRCPIPE_DIR = args.CIRCPIPE_DIR
OUTPUT_DIR = args.OUTPUT_DIR
if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
#else:
#   shutil.rmtree(OUTPUT_DIR)
#   os.makedirs(OUTPUT_DIR)
USERBPDIST = args.USERBPDIST
REFGENOME = "HG19" #args.REFGENOME
NUMBASESAROUNDJUNC = args.NUMBASESAROUNDJUNC
NumIndels = args.NumIndels
EXONS = args.EXONS
REG_INDEL_INDICES = args.REG_INDEL_INDICES
CIRCREF = args.CIRCREF
#end arg parsing


ORIG_DIR = os.path.join(CIRCPIPE_DIR,"orig")
FullStemFile = os.path.join(OUTPUT_DIR,"StemList.txt")

if os.path.isfile(FullStemFile):
## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
    print("STATUS:using existing FullStemList.txt")
else:
    print("STATUS:generating FullStemList.txt from KNIFE output directory filenames")
    subprocess.check_call("python {MACHETE}/writeStemIDFiles.py -o {ORIG_DIR} -f {OUTPUT_DIR}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR),shell=True)


include_stems = ["Human_simulated_reads1"]

#Building output directories for each Stem
processes = []
with open(FullStemFile,"r") as stem_file:
    stems = stem_file.readlines()
    for stem in sorted(stems):
        stem = stem.strip()
        print stem
        #Uncomment to skip some of the stems as defined in the include list above
        if stem not in include_stems:
            print stem+" not in include stems, skipping "+stem
            continue

        OUTPUT_DIR = args.OUTPUT_DIR
        OUTPUT_DIR += "/"+stem
        if not os.path.isdir(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        UNALIGNEDDIR = os.path.join(ORIG_DIR,"unaligned")
        GLM_DIR = os.path.join(CIRCPIPE_DIR,"circReads/glmReports")
        DistantPEDir = os.path.join(OUTPUT_DIR,"DistantPEFiles")
        if not os.path.exists(DistantPEDir):
            os.mkdir(DistantPEDir)
        FASTADIR = os.path.join(OUTPUT_DIR,"fasta")
        if not os.path.exists(FASTADIR):
            os.mkdir(FASTADIR)
        BOWTIE_DIR = os.path.join(OUTPUT_DIR,"BowtieIndex")
        if not os.path.exists(BOWTIE_DIR):
            os.mkdir(BOWTIE_DIR)
        FARJUNCDIR = os.path.join(OUTPUT_DIR,"FarJunctionAlignments")
        if not os.path.exists(FARJUNCDIR):
            os.mkdir(FARJUNCDIR)
        SECONDFARJUNCDIR = os.path.join(OUTPUT_DIR,"FarJuncSecondary")
        if not os.path.exists(SECONDFARJUNCDIR):
            os.mkdir(SECONDFARJUNCDIR)
        BadFJDir = os.path.join(OUTPUT_DIR,"BadFJ")
        if not os.path.exists(BadFJDir):
            os.mkdir(BadFJDir)
        StemFile = os.path.join(OUTPUT_DIR,"StemList.txt")
        with open(StemFile,"w") as sub_stem_file:
            sub_stem_file.write(stem+"\n")

        reportsDir = os.path.join(OUTPUT_DIR,"reports")
        if not os.path.exists(reportsDir):
            os.mkdir(reportsDir)
        LOG_DIR = os.path.join(OUTPUT_DIR,"err_and_out")
        if not os.path.exists(LOG_DIR):
            os.mkdir(LOG_DIR)

        BowtieIndelsDir = os.path.join(OUTPUT_DIR,"BowtieIndels")
        if not os.path.exists(BowtieIndelsDir):
            os.mkdir(BowtieIndelsDir)
        FarJuncIndelsDir = os.path.join(OUTPUT_DIR,"FarJuncIndels")
        if not os.path.exists(FarJuncIndelsDir):
            os.mkdir(FarJuncIndelsDir)
        AlignedIndelsDir = os.path.join(SECONDFARJUNCDIR,"AlignedIndels")
        if not os.path.exists(AlignedIndelsDir):
            os.mkdir(AlignedIndelsDir)
        IndelsHistogramDir = os.path.join(OUTPUT_DIR,"IndelsHistogram")
        if not os.path.exists(IndelsHistogramDir):
            os.mkdir(IndelsHistogramDir)
        AppendedReportsDir = os.path.join(OUTPUT_DIR,"reports/AppendedReports")
        if not os.path.exists(AppendedReportsDir):
            os.mkdir(AppendedReportsDir)
        appendGlmDir = os.path.join(GLM_DIR,"AppendGLM")
        if not os.path.isdir(appendGlmDir):
            os.mkdir(os.path.join(GLM_DIR,"AppendGLM"))
        GLM_classInputDir = os.path.join(OUTPUT_DIR,"GLM_classInput")
        if not os.path.exists(GLM_classInputDir):
            os.mkdir(GLM_classInputDir)

        OPTIONS = ""
        OPTIONS += " --circpipe-dir "+CIRCPIPE_DIR
        OPTIONS += " --output-dir "+OUTPUT_DIR
        OPTIONS += " --hg19Exons "+EXONS
        OPTIONS += " --reg-indel-indices "+REG_INDEL_INDICES
        OPTIONS += " --circref-dir "+CIRCREF

        sub_SLURM_name = os.path.join(OUTPUT_DIR,stem+"_job.sh")
        with open(sub_SLURM_name,"w") as sub_SLURM:
            sub_SLURM.write("#!/bin/bash\n")
            sub_SLURM.write("python clean_spachete_run.py"+OPTIONS+"\n")

        slurm_out_name = os.path.join(OUTPUT_DIR,"slurm_"+stem+".out")
        slurm_err_name = os.path.join(OUTPUT_DIR,"slurm_"+stem+".err")
        job_name = stem[-5:] #annoying that the full job name doesn't show
        machete_out = open(os.path.join(OUTPUT_DIR,"mach_"+stem+".out"),"w")
        machete_err = open(os.path.join(OUTPUT_DIR,"mach_"+stem+".err"),"w")

        #NOTE I'm giving 30GBs to each job, which might not be enough
        subprocess.call(["sbatch","-p","horence","--mem=30000","--time=24:00:00",
                         "-o",slurm_out_name,"-e",slurm_err_name,"-J",job_name,sub_SLURM_name],
                         stdout=machete_out,stderr=machete_err)


        """
        #This is for multithreading w/out submitting a new job
        processes.append(subprocess.Popen(["python","spachete_run.py",
                                           "--circpipe-dir",CIRCPIPE_DIR,
                                           "--output-dir",OUTPUT_DIR,
                                           "--hg19Exons",EXONS,
                                           "--reg-indel-indices",REG_INDEL_INDICES,
                                           "--circref-dir",CIRCREF]))
        """
#Force all the jobs to complete
#for p in processes:
#    p.communicate()
#print "Done"


