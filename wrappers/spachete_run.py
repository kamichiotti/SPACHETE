from argparse import ArgumentParser
from bisect import bisect_left
import subprocess
import shutil
import time
import glob
import sys
import pdb
import re
import os

# NOTE THAT AS OF JULY 2016, THIS VERSION ASSUMES THAT THERE IS ONLY ON
#   PAIR OF FASTQ FILES. SMALL CHANGES COULD FIX THIS PRESUMABLY
#   BY ADDING LOOPS OVER FILES TO 
#   THE CHANGES (MARKED BELOW) MADE IN JULY 2016.

# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for

def binary_search(a, x, lo=0, hi=None):   # can't use a to specify default for hi
    hi = hi if hi is not None else len(a) # hi defaults to len(a)   
    pos = bisect_left(a,x,lo,hi)          # find insertion position
    return (pos if pos != hi and a[pos] == x else -1) # don't walk off the end

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

full_time = time.time()
description = "Required args: --root-dir, --circpipe-dir, --output-dir, --hg19Exons, --reg-indel-indices, and --circref-dir."

parser = ArgumentParser(description=description)
parser.add_argument("--root-dir",required=True,dest="ROOT_DIR",help="Dir that run_MACHETE.sh is in")
parser.add_argument("--circpipe-dir",required=True,dest="CIRCPIPE_DIR",help="Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)")
parser.add_argument("--mode",required=True,dest="MODE",help="How to run SPACHETE/MACHETE for human --> 'hg19', for mouse --> 'mm10'")
parser.add_argument("--output-dir",required=True,dest="OUTPUT_DIR",help="Output directory for resuts. Directory path will be created recursively if it doesn't already exist. If it exists already, the directory will be deleted then created again.")
parser.add_argument("--user-bp-dist",dest="USERBPDIST",type=int,default=1000,help="Default is %(default)s.")
parser.add_argument("--num-junc-bases",dest="NUMBASESAROUNDJUNC",type=int,default=8,help="Default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70.")
parser.add_argument("--numIndels",dest="NumIndels",type=int,default=5,help="Default is %(default)s.")
#parser.add_argument("--hg19Exons",required=True,dest="EXONS",help="Path to HG19Exons. Formerly called PICKLEDIR.")
parser.add_argument("--reg-indel-indices",required=True,dest="REG_INDEL_INDICES",help="Path to files with names like hg19_junctions_reg_indels_1.1.bt2l,hg19_junctions_reg_indels_2.rev.1.bt2l ... These are sometimes in a folder called IndelIndices.")
parser.add_argument("--circref-dir",required=True,dest="CIRCREF",help="Path to reference libraries output by KNIFE - directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices.")

args = parser.parse_args()
ROOT_DIR = args.ROOT_DIR
MACHETE_DIR = os.path.join(ROOT_DIR,'MACHETE')
SPORK_DIR = os.path.join(ROOT_DIR,'SPORK')
WRAPPERS_DIR = os.path.join(ROOT_DIR,'wrappers')

CIRCPIPE_DIR = args.CIRCPIPE_DIR
MODE = args.MODE
OUTPUT_DIR = args.OUTPUT_DIR
sys.stdout.write(OUTPUT_DIR+"\n")
if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

USERBPDIST = args.USERBPDIST
NUMBASESAROUNDJUNC = args.NUMBASESAROUNDJUNC
NumIndels = args.NumIndels
#EXONS = args.EXONS #RB 3/13/17 Isn't being used, commenting out
REG_INDEL_INDICES = args.REG_INDEL_INDICES
CIRCREF = args.CIRCREF

#end arg parsing


ORIG_DIR = os.path.join(CIRCPIPE_DIR,"orig")
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

#RB adding timer functionality
timer_file_path = LOG_DIR+"/timed_events.txt"
def write_time(message,start_time,timer_file_path,append=True,uniform_len=70):
    """
    Goal: write a timed event out to a file
    Arguments:
        message is the text of the job to put next to the time (i.e. 'Time to align reads')
        start_time is the start of the job gotten by time.time()
        time_file_path is the full path to the timer store file
        append is a boolean saying whether or not to append to the file vs. overwriting (default append)
        uniform_len is the max message length to allow aligning the timer output file

    Returns:
        nothing (just writes to the timer file)
    """
    seconds_duration = float(time.time()-start_time)
    minutes_duration = seconds_duration/60
    hours_duration = seconds_duration/3600
    seconds_str= ("{:3.2f}".format(seconds_duration)).rjust(5," ")
    minutes_str= ("{:3.2f}".format(minutes_duration)).rjust(5," ")
    hours_str= ("{:3.2f}".format(hours_duration)).rjust(5," ")

    if len(message) < uniform_len:
        message += " "*(uniform_len-len(message))
    else:
        message = message[:uniform_len]
    time_out_str  = message+":    "
    time_out_str += seconds_str+":seconds    "
    time_out_str += minutes_str+":minutes    "
    time_out_str += hours_str+":hours\n"
    timer_file = open(timer_file_path,"a") if append else open(timer_file_path,"w")
    timer_file.write(time_out_str)
    timer_file.close()

#RB filter out used ids from glm class files
def filter_glm_class_file(glm_input_name,used_ids_name,timer_file_path):
    start_time = time.time()
    filtered_glm_outupt_name = glm_input_name+".tmp"

    sys.stdout.write("Starting filtering: "+glm_input_name+"\n")
    glm_output = open(filtered_glm_outupt_name,"w")
    with open(used_ids_name,"r") as used_ids_file:
        used_ids = []
        for used_id in used_ids_file:
            clean_used_id = used_id.split(" ")[0][1:].strip()
            used_ids.append(clean_used_id)
        used_ids.sort()

        for glm_line in open(glm_input_name,"r"):
            glm_read_id = glm_line.split("\t")[0].strip()
            #If the glm input is NOT in the used ids names
            #then write to the glm output
            res_ind = binary_search(used_ids,glm_read_id)
            if res_ind == -1:
                glm_output.write(glm_line.strip()+"\n")

    glm_output.flush()
    glm_output.close()

    prefilter_name = glm_input_name+".prefilter"
    os.rename(glm_input_name,prefilter_name)
    os.rename(filtered_glm_outupt_name,glm_input_name)
    write_time("Filter out reads used in: "+glm_input_name,start_time,timer_file_path)



subprocess.check_call("rm -f {LOG_DIR}/*".format(LOG_DIR=LOG_DIR),shell=True)
subprocess.check_call("rm -f {LOG_DIR}/MasterError.txt".format(LOG_DIR=LOG_DIR),shell=True)

start_time = time.time()
## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
if os.path.isfile(StemFile):
    print("STATUS:using existing StemList.txt")
    write_time("Use existing unique file name pairs",start_time,timer_file_path,False) #false overwrites the timer file
    for stem in open(StemFile,"r"):
        write_time("--> Found "+stem.strip(),time.time(),timer_file_path)
else:
    print("STATUS:generating StemList.txt from KNIFE output directory filenames")
    subprocess.check_call("python {MACHETE_DIR}/writeStemIDFiles.py -o {ORIG_DIR} -f {OUTPUT_DIR}".format(MACHETE_DIR=MACHETE_DIR,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR),shell=True)
    write_time("Detect all unique file name pairs",start_time,timer_file_path,False) #false overwrites the timer file

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
NUM_FILES = len(open(StemFile,"r").readlines())
print("NUM_FILES: "+str(NUM_FILES))

# Run SPORK
#ALIGN_PARDIR = "/".join(CIRCPIPE_DIR.split("/")[:-1])+"/"
#DATASET_NAME = CIRCPIPE_DIR.split("/")[-1]
#MODE = "complete"
NTRIM = "13"
#DENOVOCIRC = "1"
#NUM_FLANKING = "150"
#This should be blocking
start_time = time.time()
SPORK_STEM_NAME = open(StemFile,"r").readline().strip()
SPORK_OPTIONS = ""
SPORK_OPTIONS += "--input-dir "+CIRCPIPE_DIR+" "
SPORK_OPTIONS += "--output-dir "+OUTPUT_DIR+" "
SPORK_OPTIONS += "--stem-name "+SPORK_STEM_NAME+" "
spork_process = subprocess.Popen(["python",os.path.join(SPORK_DIR,"SPORK_main.py"),
                                    "--root-dir",ROOT_DIR,
                                    "--input-dir",CIRCPIPE_DIR,
                                    "--output-dir",OUTPUT_DIR,
                                    "--ref-dir",CIRCREF,
                                    "--mode",MODE,
                                    "--stem-name",SPORK_STEM_NAME],
                                    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
spork_out,spork_err = spork_process.communicate()
spork_ret_code = spork_process.returncode
sys.stdout.write(spork_out)
sys.stdout.flush()
sys.stderr.write(spork_err)
sys.stderr.flush()
write_time("Run spork",start_time,timer_file_path)

#RB NOTE it looks like it is possible for SPORK to yield an empty "fusion-fasta" file to machete, which breaks it
#RB I'm first going to see if the file is there, then I'm going to check if the file is empty, if so I'll just not run MACHETE
SPORK_FASTA= os.path.join("spork_out","novel_junctions.fasta")
fusion_fasta_path = os.path.join(OUTPUT_DIR,SPORK_FASTA)
if not os.path.isfile(fusion_fasta_path):
    if str(spork_ret_code) == "0":
        sys.stderr.write("SPORK: exited fine but didn't create 'fusions-fasta' file, likely problem with input file\n")
        sys.stdout.write("SPORK: exited fine but didn't create 'fusions-fasta' file, likely problem with input file\n")
        sys.exit(0)

    elif str(spork_ret_code) == "1":
        sys.stderr.write("SPORK ERROR: Some error in SPORK, didn't create 'fusions-fasta' file, exiting immediately\n")
        sys.stdout.write("SPORK ERROR: Some error in SPORK, didn't create 'fusions-fasta' file, exiting immediately\n")
        sys.exit(1)

num_fusion_lines = 0
with open(fusion_fasta_path,'r') as fusion_fasta:
    num_fusion_lines = len(fusion_fasta.readlines())

if num_fusion_lines == 0:
    if str(spork_ret_code) == "0":
        sys.stderr.write("SPORK: exited fine but didn't identify any junctions at all to pass to MACHETE, exiting immediately\n")
        sys.stdout.write("SPORK: exited fine but didn't identify any junctions at all to pass to MACHETE, exiting immediately\n")
        sys.exit(0)
    elif str(spork_ret_code) == "1":
        sys.stderr.write("SPORK ERROR: didn't identify any junctions at all to pass to MACHETE, exiting immediately\n")
        sys.stdout.write("SPORK ERROR: didn't identify any junctions at all to pass to MACHETE, exiting immediately\n")
        sys.exit(1)


#Make the bowtie index building call on the spork fasta
start_time = time.time()
print("STATUS:make FJ bowtie indices for each experiment")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_5FJIndexing.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_5FJIndexing.txt"),"w")
    #print(cmd)
    cmd = "{MACHETE_DIR}/linkfastafiles_SPORK.sh {OUTPUT_DIR} {index} {SPORK_FASTA} {INSTALLDIR}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index,SPORK_FASTA=SPORK_FASTA,INSTALLDIR=MACHETE_DIR)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Link fastas inds",start_time,timer_file_path)


### END CHANGES for SPORK JS
# Added July 5 2016
# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
#j8new_id this is B2
start_time = time.time()
print("STATUS:align unaligned reads to FJ index:  - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq")
print("this align of unaligned reads was added in july 2016; it takes place before the original one")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_7newAlignFJ.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_7newAlignFJ.txt"),"w")
    cmd = "{MACHETE_DIR}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("AlignUnalignedtoFJ inds",start_time,timer_file_path)

#
#
###make FJ naive report
## FarJuncNaiveReport.sh is a shell script that calls the python script FarJuncNaiveReport.py to generate the "Naive Reports".  Inputs include the MACHETE output directory, paths to the KNIFE alignment files, the amount a read should overlap the junction in order to be considered a "true" junctional alignment, and the MACHETE installation directory.
## see the FarJuncNaiveReport.sh for more info on FarJuncNaiveReport.py and details about how alignments are selected as "true" or "false", and how the a p value is calculated.
## The rate of true or anomaly alignments and p values are output to FJDir/reports/<STEM>_naive_report.txt.  Specific read ID's are also tracked and information on them can be found in FJDir/reports/IDs_<STEM>.txt.


#
# Added July 5 2016
#j9new_id this is C
start_time = time.time()
print("STATUS:make naive rpt, new as of jul 2016, earlier than previous one ")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_8newNaiveRpt.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_8NnewaiveRpt.txt"),"w")
    cmd = "{MACHETE_DIR}/FarJuncNaiveReport.sh {OUTPUT_DIR} {ORIG_DIR} {NUMBASESAROUNDJUNC} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("FarJuncNaiveReport inds",start_time,timer_file_path)

#
#

# Added July 5 2016
# ASSUMES THERE IS ONLY ONE PAIR OF FASTQ FILES! (ok for spachete_feeder)
start_time = time.time()
i=1
stdout = open(os.path.join(LOG_DIR,str(i) + "assumingonefile_out_getStem.txt"),"w")
stderr = open(os.path.join(LOG_DIR,str(i) + "assumingonefile_err_getStem.txt"),"w")
stemCmd = "awk 'FNR == '{i}' {{print $1}}' {StemFile}".format(i=i,StemFile=StemFile)
popen = subprocess.Popen(stemCmd,stdout=stdout,stderr=stderr,shell=True)
popen.communicate() 
stdout.close()
stderr.close()
retcode = popen.returncode
if retcode:
    raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

STEM_ASSUMING_ONE_FILE = open(stdout.name,'r').read().strip()

print("STATUS:STEM_ASSUMING_ONE_FILE is:" + STEM_ASSUMING_ONE_FILE + "\n")

fasta_stem_dir = os.path.join(FASTADIR,STEM_ASSUMING_ONE_FILE)
write_time("Filter stem file w/ awk",start_time,timer_file_path)


## Now call parse_to_remove_FJ.py

stdout = open(os.path.join(LOG_DIR,"parse_out.txt"),"w")
stderr = open(os.path.join(LOG_DIR,"parse_err.txt"),"w")

cmd="python {MACHETE_DIR}/parse_to_remove_FJ.py --stem {STEM_ASSUMING_ONE_FILE} --outputdir {OUTPUT_DIR}".format(MACHETE_DIR=MACHETE_DIR, STEM_ASSUMING_ONE_FILE=STEM_ASSUMING_ONE_FILE, OUTPUT_DIR=OUTPUT_DIR)

print(cmd)
start_time = time.time()
subprocess.call(cmd, shell=True, stdout=stdout, stderr=stderr)
write_time("Parse_to_remove_FJ",start_time,timer_file_path)

stdout.close()
stderr.close()

# Added July 5 2016
##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.
## The script linkfastafiles.sh uses linux to concatenate the <FJDir>/fasta/<STEM>/<STEM>_chr1,2,3,...,X,Y_FarJunctions.fa into a single large fasta <FJDir>/fasta/<STEM>_FarJunctions.fa.
## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>/<STEM>_FJ_Index
#j6anew_id- this is B1
start_time = time.time()
print("STATUS:make FJ bowtie indices for each experiment")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_5secondnewFJIndexing.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_5secondnewFJIndexing.txt"),"w")
    print(cmd)
    cmd = "{MACHETE_DIR}/linkfastafiles.sh {OUTPUT_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Compile FJ fasta inds",start_time,timer_file_path)

# Added July 5 2016
# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
#j8new_id this is B2
start_time = time.time()
print("STATUS:align unaligned reads to FJ index:  - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq")
print("this align of unaligned reads was added in july 2016; it takes place before the original one")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_7secondnewAlignFJ.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_7secondnewAlignFJ.txt"),"w")
    cmd = "{MACHETE_DIR}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("AlignUnalignedtoFJ fasta inds",start_time,timer_file_path)


## If there is homology between a FarJunctions fasta sequence and the genome or transcriptome or a linear junction or circular junction, then the fusion read is less likely.  Alignments of the FarJunctions fasta sequences to the KNIFE reference indices, genome, transcriptome, linear junctions (reg), and scrambled junctions (junc) are created with two different bowtie parameters.  Bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align. These are just aligning the FJ Fasta to the bad juncs with various alignment parameters. Any junctions aligning to here will eventually be tagged as "BadFJ=1" in the final reports whereas if junctions don't align, they will receive a "BadFJ=0" in the final reports.
#RB 12/5/16: Added the if elif for MODE type
genomeIndex = ""
transcriptomeIndex = ""
regIndex = ""
juncIndex = ""

if MODE.lower() == "hg19":
    genomeIndex = os.path.join(CIRCREF,"hg19_genome")
    transcriptomeIndex = os.path.join(CIRCREF,"hg19_transcriptome")
    regIndex = os.path.join(CIRCREF,"hg19_junctions_reg")
    juncIndex = os.path.join(CIRCREF,"hg19_junctions_scrambled")

elif MODE.lower() == "mm10":
    genomeIndex = os.path.join(CIRCREF,"mm10_genome")
    transcriptomeIndex = os.path.join(CIRCREF,"mm10_transcriptome")
    regIndex = os.path.join(CIRCREF,"mm10_junctions_reg")
    juncIndex = os.path.join(CIRCREF,"mm10_junctions_scrambled")

   

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --np 0 --rdg 50,50 --rfg 50,50"

for i in range(1,NUM_FILES + 1):
    start_time = time.time()
    stdout = open(os.path.join(LOG_DIR,str(i) + "_out_getStem.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(i) + "_err_getStem.txt"),"w")
    stemCmd = "awk 'FNR == '{i}' {{print $1}}' {StemFile}".format(i=i,StemFile=StemFile)
    sub_start_time = time.time()
    popen = subprocess.Popen(stemCmd,stdout=stdout,stderr=stderr,shell=True)
    popen.communicate() 
    write_time("--filter with awk "+str(i),sub_start_time,timer_file_path)
    stdout.close()
    stderr.close()
    retcode = popen.returncode
    if retcode:
        raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))
    STEM = open(stdout.name,'r').read().strip()
    FarJuncFasta = glob.glob(os.path.join(FASTADIR,STEM + "*FarJunctions.fa")) #should be a list of size 1
    FarJuncFasta = FarJuncFasta[0]
    BadFJStemDir =os.path.join(BadFJDir,STEM)
    if  os.path.isdir(BadFJStemDir):
        shutil.rmtree(BadFJStemDir)
    os.mkdir(BadFJStemDir)

    BadFJver2Dir=os.path.join(OUTPUT_DIR,"BadFJ_ver2",STEM)
    if os.path.isdir(BadFJver2Dir):
        shutil.rmtree(BadFJver2Dir)
    os.makedirs(BadFJver2Dir)

    r1file=os.path.join(BadFJver2Dir,STEM + "_FarJunctions_R1.fa")
    r2file=os.path.join(BadFJver2Dir,STEM + "_FarJunctions_R2.fa")
    #Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM
    # is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file 
    # -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the 
    # last 40 non-N reads from the sequence.
    
    #j7_id=
    print("STATUS:Identify Bad FJ's")
    processes = {}
    stdout = open(os.path.join(LOG_DIR,str(STEM) + "_out_6BadJunc.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(STEM) + "_err_6BadJunc.txt"),"w")
    cmd = "{MACHETE_DIR}/LenientBadFJ_SLURM.sh {FarJuncFasta} {BadFJver2Dir} {OUTPUT_DIR} {MACHETE_DIR}".format(MACHETE_DIR=MACHETE_DIR,FarJuncFasta=FarJuncFasta,BadFJver2Dir=BadFJver2Dir,OUTPUT_DIR=OUTPUT_DIR)
    print(cmd)
    sub_start_time = time.time()
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    stdout.close()
    stderr.close()
    retcode = popen.returncode
    write_time("--lenientBadFJ_SLURM"+str(i),sub_start_time,timer_file_path)
    if retcode:
        raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))
    BadFJtoGenomeFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoGenome.sam")
    sub_start_time = time.time()
    if os.path.exists(BadFJtoGenomeFile):
        print("{BadFJtoGenomeFile} exists. To realign, please manually delete this file first".format(BadFJtoGenomeFile=BadFJtoGenomeFile))
        write_time("--using existing BADFJtoGenomeFile "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
        fasta = os.path.join(FASTADIR,"{STEM}_FarJunctions.fa".format(STEM=STEM))
        cmd = "{MACHETE_DIR}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {genomeIndex} {fasta} {BadFJtoGenomeFile}".format(MACHETE_DIR=MACHETE_DIR,BOWTIEPARAM=BOWTIEPARAM,genomeIndex=genomeIndex,fasta=fasta,BadFJtoGenomeFile=BadFJtoGenomeFile)
        print(cmd)
        print("BadFJ to genome")
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        stdout.close()
        stderr.close()
        retcode = popen.returncode
        write_time("--Make BadFJ to genome "+str(i),sub_start_time,timer_file_path)
        if retcode:
            raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

    BadFJtotranscriptomeFile = os.path.join(BadFJStemDir,STEM + "__BadFJtotranscriptome.sam")
    sub_start_time = time.time()
    if os.path.exists(BadFJtotranscriptomeFile):
        print("{BadFJtotranscriptomeFile} exists.  To realign, please manually delete this file first".format(BadFJtotranscriptomeFile=BadFJtotranscriptomeFile))
        write_time("--BadFJtotranscriptome found "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
        cmd = "{MACHETE_DIR}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {transcriptomeIndex} {fasta} {BadFJtotranscriptomeFile}".format(MACHETE_DIR=MACHETE_DIR,BOWTIEPARAM=BOWTIEPARAM,transcriptomeIndex=transcriptomeIndex,fasta=fasta,BadFJtotranscriptomeFile=BadFJtotranscriptomeFile)
        print(cmd)
        print("BadFJ to transcriptome")
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        stdout.close()
        stderr.close()
        retcode = popen.returncode
        write_time("--Make BadFJ to transcriptome "+str(i),sub_start_time,timer_file_path)
        if retcode:
            raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

    sub_start_time = time.time()
    BadFJtoRegFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoReg.sam")
    if os.path.exists(BadFJtoRegFile):
        print("{BadFJtoRegFile} exists. To realign, please manually delete this file first.".format(BadFJtoRegFile=BadFJtoRegFile))
        write_time("--BadFJtoRegFile exists "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
        cmd = "{MACHETE_DIR}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {regIndex} {fasta} {BadFJtoRegFile}".format(MACHETE_DIR=MACHETE_DIR,BOWTIEPARAM=BOWTIEPARAM,regIndex=regIndex,fasta=fasta,BadFJtoRegFile=BadFJtoRegFile)
        print(cmd)
        print("STATUS:BadFJ to reg")
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        stdout.close()
        stderr.close()
        retcode = popen.returncode
        write_time("--Making BadFJtoRegFile "+str(i),sub_start_time,timer_file_path)
        if retcode:
            raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))

    sub_start_time = time.time()
    BadFJtoJuncFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoJunc.sam")
    if os.path.exists(BadFJtoJuncFile):
        print("{BadFJtoJuncFile} exists. To realign, please manually delete this file first".format(BadFJtoJuncFile=BadFJtoJuncFile))
        write_time("--BadFJtoJuncFile exists "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJStemDir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJStemDir,"err.txt"),"w")
        print("STATUS:BadFJ to junc: ")
        cmd = "{MACHETE_DIR}/BowtieFJAligner.batch.sh \"{BOWTIEPARAM}\" {juncIndex} {fasta} {BadFJtoJuncFile}".format(MACHETE_DIR=MACHETE_DIR,BOWTIEPARAM=BOWTIEPARAM,juncIndex=juncIndex,fasta=fasta,BadFJtoJuncFile=BadFJtoJuncFile)
        print(cmd)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        stdout.close()
        stderr.close()
        retcode = popen.returncode
        write_time("--Making BadFJtoJuncFile "+str(i),sub_start_time,timer_file_path)
        if retcode:
            raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}.".format(cmd=stemCmd,retcode=retcode,stdout=stdout.name,stderr=stderr.name))



    ## Read gaps are disallowed in the first version of BadJuncs.  A second version of BadJuncs was created to also find genome/reg/junc/transcriptome alignments with gapped alignments.
    ## For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.
    badfj2_thresh = 500000 #RB 11/18/16: Updated this from 50,000 to 500,000
    genomeBOWTIEPARAM = "--no-unal --no-mixed --no-sq -p 8 -I 0 -X {badfj2_thresh} -f --ff -x {genomeIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtoGenome.sam".format(badfj2_thresh=badfj2_thresh,genomeIndex=genomeIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)
    transcriptomeBOWTIEPARAM = "--no-unal --no-mixed  --no-sq -p 8 -I 0 -X {badfj2_thresh} -f --ff -x {transcriptomeIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtotranscriptome.sam".format(badfj2_thresh=badfj2_thresh,transcriptomeIndex=transcriptomeIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)
    regBOWTIEPARAM = "--no-unal --no-mixed --no-sq -p 8 -I 0 -X {badfj2_thresh} -f --ff -x {regIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtoReg.sam".format(badfj2_thresh=badfj2_thresh,regIndex=regIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)
    juncBOWTIEPARAM ="--no-unal --no-mixed --no-sq -p 8 -I 0 -X {badfj2_thresh} -f --ff -x {juncIndex} -1 {r1file} -2 {r2file} -S {BadFJver2Dir}/{STEM}_BadFJtoJunc.sam".format(badfj2_thresh=badfj2_thresh,juncIndex=juncIndex,r1file=r1file,r2file=r2file,BadFJver2Dir=BadFJver2Dir,STEM=STEM)

    
    #do the bowtie alignment for each of the BadFJ Ver 2 indices above
    sub_start_time = time.time()
    BadFJtoGenomeFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtoGenome.sam")
    if os.path.exists(BadFJtoGenomeFile):
        print( "{BadFJtoGenomeFile} exists. To realign, please manually delete this file first".format(BadFJtoGenomeFile=BadFJtoGenomeFile))
        write_time("--BadFJtoGenomeFile v2 exists "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
        #RB: This call is failing 3/1/17, not finding STEM_FarJunctions_R1.fa (the R1 file)
        cmd = "{MACHETE_DIR}/BowtieAligner_BadFJv2.sh \"{genomeBOWTIEPARAM}\"".format(MACHETE_DIR=MACHETE_DIR,genomeBOWTIEPARAM=genomeBOWTIEPARAM)
        print(cmd)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes = {}
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
        checkProcesses(processes)
        write_time("--Making BadFJtoGenomeFile v2 "+str(i),sub_start_time,timer_file_path)
    
    print("STATUS:BadFJ_ver2 to genome")
    BadFJtotranscriptomeFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtotranscriptome.sam")
    sub_start_time = time.time()
    if os.path.exists(BadFJtotranscriptomeFile):
        print("{BadFJtotranscriptomeFile} exists. To realign, please manually delete this file first.".format(BadFJtotranscriptomeFile=BadFJtotranscriptomeFile))
        write_time("--BadFJtoTranscriptomeFile v2 exists "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
        cmd = "{MACHETE_DIR}/BowtieAligner_BadFJv2.sh \"{transcriptomeBOWTIEPARAM}\"".format(MACHETE_DIR=MACHETE_DIR,transcriptomeBOWTIEPARAM=transcriptomeBOWTIEPARAM)
        print(cmd)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes = {}
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
        checkProcesses(processes)
        write_time("--Making BadFJtoTranscriptomeFile v2 "+str(i),sub_start_time,timer_file_path)
    print("STATUS:BadFJ_ver2 to transcriptome")

    BadFJtoRegFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtoReg.sam")
    sub_start_time = time.time()
    if os.path.exists(BadFJtoRegFile):
        print("{BadFJtoRegFile} exists. To realign, please manually delete this file first.".format(BadFJtoRegFile=BadFJtoRegFile))
        write_time("--BadFJtoRegFile v2 exists "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
        cmd = "{MACHETE_DIR}/BowtieAligner_BadFJv2.sh \"{regBOWTIEPARAM}\"".format(MACHETE_DIR=MACHETE_DIR,regBOWTIEPARAM=regBOWTIEPARAM)
        print(cmd)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes = {}
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
        checkProcesses(processes)
        write_time("--Making BadFJtoRegFile v2 "+str(i),sub_start_time,timer_file_path)
    print("STATUS:BadFJ_ver2 to reg")

    BadFJtoJuncFile = os.path.join(BadFJver2Dir,STEM + "_BadFJtoJunc.sam")
    sub_start_time = time.time()
    if os.path.exists(BadFJtoJuncFile):
        print("{BadFJtoJuncFile} exists. To realign, please manually delete this file first.".format(BadFJtoJuncFile=BadFJtoJuncFile))
        write_time("--BadFJtoJuncFile v2 exists "+str(i),sub_start_time,timer_file_path)
    else:
        stdout = open(os.path.join(BadFJver2Dir,"out.txt"),"w")
        stderr = open(os.path.join(BadFJver2Dir,"err.txt"),"w")
        cmd = "{MACHETE_DIR}/BowtieAligner_BadFJv2.sh \"{juncBOWTIEPARAM}\"".format(MACHETE_DIR=MACHETE_DIR,juncBOWTIEPARAM=juncBOWTIEPARAM)
        print(cmd)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes = {}
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
        checkProcesses(processes)
        write_time("--Making BadFJtoJuncFile v2 "+str(i),sub_start_time,timer_file_path)
    print("STATUS:BadFJ_ver2 to junc")
    write_time("Big alignment loop "+str(i),start_time,timer_file_path)



# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
#j8_id
start_time = time.time()
print("STATUS:align unaligned reads to FJ index:  - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_7AlignFJ.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_7AlignFJ.txt"),"w")
    cmd = "{MACHETE_DIR}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("AlignUnalignedtoFJ again",start_time,timer_file_path)


#
#
###make FJ naive report
## FarJuncNaiveReport.sh is a shell script that calls the python script FarJuncNaiveReport.py to generate the "Naive Reports".  Inputs include the MACHETE output directory, paths to the KNIFE alignment files, the amount a read should overlap the junction in order to be considered a "true" junctional alignment, and the MACHETE installation directory.
## see the FarJuncNaiveReport.sh for more info on FarJuncNaiveReport.py and details about how alignments are selected as "true" or "false", and how the a p value is calculated.
## The rate of true or anomaly alignments and p values are output to FJDir/reports/<STEM>_naive_report.txt.  Specific read ID's are also tracked and information on them can be found in FJDir/reports/IDs_<STEM>.txt.


#
#j9_id
start_time = time.time()
print("STATUS:make naive rpt")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_8NaiveRpt.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_8NaiveRpt.txt"),"w")
    cmd = "{MACHETE_DIR}/FarJuncNaiveReport.sh {OUTPUT_DIR} {ORIG_DIR} {NUMBASESAROUNDJUNC} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Make naive report again",start_time,timer_file_path)

##
####### GLM #######
##Make Class input files for GLM
##
###
#j15a_id
start_time = time.time()
print("STATUS:make FJ class input files")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_15FJforGLM.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_15FJforGLM.txt"),"w")
    cmd = "{MACHETE_DIR}/parse_FJ_ID_for_GLM.sh {OUTPUT_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Parse FJ ID for GLM",start_time,timer_file_path)


###### ESTIMATE LIGATION ARTIFACT
## Ligation artifact refers to the rate at which cDNA are artificially ligated, producing what appears to be a true junctional alignment.  In this step, reads that remained unaligned to the Far Junctions Bowtie index that are located in FarJuncSecondary are then aligned to new indel fasta indices under the assumption that in any local area, the rate of false ligation of a junction is similar to the rate of false ligation of any two cDNAs that do not end at exon boundaries.
#
###
### MakeIndelFiles.sh is a shell script that calls the python script AddIndelsToFasta.py.  It takes the FarJunctions fasta files as inputs (FJDir/fasta/<STEM>_FarJunctions.fa) and outputs five files called FJDir/FarJuncIndels/<STEM>/<STEM>_FJ_Indels_1,2,3,4,5.fa where the numbers 1-5 indicate the N's inserted on each side of the breakpoint or deletions on each side of the breakpoint.  For example, the FJ_indels_3 file is the same as the FarJunctions.fa file except every sequence has 3 N's on each side of the breakpoint (total of 6 N's inserted at the breakpoint), or 3 bases are deleted from each exon on each side of the breakpoint (total of 6 deletions at the breakpoint).
####

#j10_id
start_time = time.time()
print("STATUS:make indel files")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_10FJIndels.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_10FJIndels.txt"),"w")
    cmd = "{MACHETE_DIR}/MakeIndelFiles.sh {OUTPUT_DIR} {NumIndels} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,NumIndels=NumIndels,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("MakeIndelFiles",start_time,timer_file_path)

#
#
# make Bowtie Indices for Far Junc Indel files
## The shell script BowtieIndexFJIndels.sh calls bowtie to index the indels_N.fa files that were created in the previous step.  The indels are output to the directory FJDir/BowtieIndels/<STEM>/<STEM>_Indels_N where N is the number of indels in that index.
#
print("STATUS:index indels")
for i in range(1,NumIndels + 1):
    start_time = time.time()
    processes = {}
    for index in range(1,NUM_FILES + 1):
        stdout = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_out_11indexindels.txt".format(i=i,index=index)),"w")
        stderr = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_err_11indexindels.txt".format(i=i,index=index)),"w")
        cmd = "{MACHETE_DIR}/BowtieIndexFJIndels.sh {FarJuncIndelsDir} {i} {BowtieIndelsDir} {OUTPUT_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,FarJuncIndelsDir=FarJuncIndelsDir,BowtieIndelsDir=BowtieIndelsDir,OUTPUT_DIR=OUTPUT_DIR,i=i,index=index)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
    checkProcesses(processes)
    write_time("Index indels "+str(i),start_time,timer_file_path)


##Align FarJuncSecondary (unaligned to FJ index) to FJ indels

# This section calls the shell BowtieAlignFJIndels.sh to align the fq files FJdir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_1/2.fq to the Bowtie2 indices of the far junction indels created in the previous step.  Aligned indels are stored in FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam where N is the number of indels in the Bowtie index where the reads aligned.  The bowtie parameters include a max of ~4 mismatches / 100 basepairs, max #N is the read length, and prohibits gapped alignments or gaps in the read.

#
##
print("STATUS:align to indels")
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
for i in range(1,NumIndels + 1):
    start_time = time.time()
    processes = {}
    for index in range(1,NUM_FILES + 1):
        stdout = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_out_12alignindels.txt".format(i=i,index=index)),"w")
        stderr = open(os.path.join(LOG_DIR,"NumIndes{i}_{index}_err_12alignindels.txt".format(i=i,index=index)),"w")
        cmd = "{MACHETE_DIR}/BowtieAlignFJIndels.sh {OUTPUT_DIR} \"{BOWTIEPARAMETERS}\" {i} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,BOWTIEPARAMETERS=BOWTIEPARAMETERS,i=i,index=index)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
    checkProcesses(processes)
    write_time("Align indels "+str(i),start_time,timer_file_path)


## Calls FindAlignmentArtifact_SLURM.sh which is a shell that calls MakeIndelsHisto.py.  The MakeIndelsHisto.py script reads in the aligned indels from the sam files FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam.  It concatenates all the indels_1,2,3,4,5.sam files into a larger file All_<STEM>_1/2_indels.sam in the same directory as the original sam files. In the All_<STEM>_1/2_indels.sam files, any indels that did not overlap the junction by the user specified # base pairs around the breakpoint are removed. Additionally, since the FarJuncSecondary files were aligned independently to the indels_1-5 bowtie indices, the same read could align to multiple indices.  Therefore for each read, the read with the best alignment score is placed in the All_<STEM>_1/2_indels sam file, and all other alignments to other indices are discarded.
#  Then the python script checks the FJdir/FarJunctionAlignments/<STEM>/ sam files and creates an array for each junction.  The array is of length 2*#indels+1. In the case of 5 indels, the length is 11 and represents the number reads that aligned to the junction with [5Del, 4Del, 3Del, 2Del, 1Del, aligned to junction exactly, 1Ins, 2Ins, 3Ins, 4Ins, 5Ins]
## the junction name and this array are output to FJDir/IndelsHistogram/indels_<STEM>_1/2.txt.  These outputs will be used to generate the Appended reports later.
#j14_id
start_time = time.time()
print("STATUS:make indels histo")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_13filterIndels.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_13filterIndels.txt"),"w")
    cmd = "{MACHETE_DIR}/FindAlignmentArtifact_SLURM.sh {OUTPUT_DIR} {NUMBASESAROUNDJUNC} {NumIndels} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,NumIndels=NumIndels,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Make indel histogram",start_time,timer_file_path)


###
####
##### REG INDELS ################
##
# Align unaligned files to the expanded reg junctions with indels
AlignedIndels = os.path.join(CIRCPIPE_DIR,"orig/RegIndelAlignments")
if not os.path.exists(AlignedIndels):
    os.makedirs(AlignedIndels)

#  To train the GLM, indel alignments are also created for the linear junctions.
#The reference index of indels to the linear junctions is static and has already been created
#It is referenced above as "REG_INDEL_INDICES" on line 44.
#The script AlignRegIndels calls bowtie to align reads that were unaligned the the KNIFE indices (in KNIFEdir/orig/unaligned/*.fq)
#to the REG_INDEL_INDICES, with the parameters of 1) approx 4 mismatches / 100 bases, maximum number N's = readlength,
#and no gapped alignments or read gaps.
#j16_id
print("STATUS:Aligning unaligned files to linear junc indels")
BOWTIEPARAMETERS = "--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
for i in range(1,NumIndels + 1):
    start_time
    processes = {}
    for index in range(1,NUM_FILES + 1):
        stdout = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_out_15AlignRegIndels.txt".format(i=i,index=index)),"w")
        stderr = open(os.path.join(LOG_DIR,"NumIndels{i}_{index}_err_15AlignRegIndels.txt".format(i=i,index=index)),"w")
        cmd = "{MACHETE_DIR}/AlignUnalignedtoRegIndel.sh {CIRCPIPE_DIR} {i} {OUTPUT_DIR} \"{BOWTIEPARAMETERS}\" {REG_INDEL_INDICES} {index}".format(MACHETE_DIR=MACHETE_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,i=i,OUTPUT_DIR=OUTPUT_DIR,BOWTIEPARAMETERS=BOWTIEPARAMETERS,REG_INDEL_INDICES=REG_INDEL_INDICES,index=index)
        popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
        processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
    checkProcesses(processes)
    write_time("Align linear indels "+str(i),start_time,timer_file_path)

###
####
######## MAKE REG AND FJ INDELS CLASS OUTPUT FILES ###########
####
###
### reg indels class output file
#Calls RegIndelsClassID.sh shell which calls RegIndels_ClassIDFile.py.
# Inputs are the regular indel alignments located at KNIFEdir/orig/RegIndelAlignments/<STEM>/unaligned_<STEM>_1/2_indel1,2,3,4,5.sam.  These are concatenated into a single file in the same directory called All_<STEM>_1/2_Regindels.sam.  In the concatenation step, like for the far junctions indels, any reads are omitted if they fail to overlie the junction by the user specified overlap, and also if a read aligns to multiple indel indices, the one with the best alignment score is put into the concatenated file.
## Then, the partner reads of all reads from All_<STEM>_1/2_Regindels.sam are identified and labeled as "good" or "bad".  If a read partner is found in genome, it has priority over transcriptome, which has priority over reg, and finally junc.  A far junction R2 cannot be found in another dictionary, as FJ reads are generated from previously unaligned reads.  If the read partner is in genome, it must be located on the same chromosome, on the opposite reference strand from R1.  If a read partner is in reg, then the downstream exon must be upstream of the uptstream reg indel exon, or the upstream read partner exon must be downstream of the downstream reg indel exon, on the same chromosome. Reference strands must be opposite.    If the partner is in junc, then the reg indel alignment must be located inside the circle created by the scrambled junction, and on the opposite reference strand.  In this manner, class input files are generated for the reg indels, which are located at KNIFE dir/circReads/ids/<STEM>_output_RegIndel.txt
#

#j18_id
start_time = time.time()
print("STATUS:Reg Indels Class Output")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_18RegIndelsClassOutput.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_18RegIndelsClassOutput.txt"),"w")
    cmd = "{MACHETE_DIR}/RegIndelsClassID.sh {OUTPUT_DIR} {CIRCPIPE_DIR} {NUMBASESAROUNDJUNC} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("RegIndelsClassID",start_time,timer_file_path)


# FJ indels class output file
## This script calls FJIndelsClassID.sh
## This takes the FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/All_<STEM>_1/2_FJindels.sam and identifies read partners.  The same criteria to identify read partners as FarJuncNaiveReport.sh are applied (see above).
## Output files are placed into FJDir/GLM_classInput/<STEM>_output_FJIndels.txt
#j19_id
start_time = time.time()
print("STATUS:FJ Indels Class Output")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_19FJIndelsClassOutput.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_19FJIndelsClassOutput.txt"),"w")
    cmd = "{MACHETE_DIR}/FJIndelsClassID.sh {OUTPUT_DIR} {CIRCPIPE_DIR} {NUMBASESAROUNDJUNC} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("FJIndelsClassID",start_time,timer_file_path)
print("STATUS:FJ Indels Class Output")


##Filter out reads used to build the SPORK junctions before they get to the GLM.
used_ids_name = os.path.join(OUTPUT_DIR,"spork_out","used_read_ids.txt")
glm_input_name = os.path.join(OUTPUT_DIR,"GLM_classInput",SPORK_STEM_NAME+"_1__output_FJ.txt")
filter_glm_class_file(glm_input_name,used_ids_name,timer_file_path)

glm_input_name = os.path.join(OUTPUT_DIR,"GLM_classInput",SPORK_STEM_NAME+"_output_FJIndels.txt")
filter_glm_class_file(glm_input_name,used_ids_name,timer_file_path)

glm_input_name = os.path.join(OUTPUT_DIR,"GLM_classInput",SPORK_STEM_NAME+"_temp_output_FJIndels.txt")
filter_glm_class_file(glm_input_name,used_ids_name,timer_file_path)


###### RUN GLM ###########################
#
## Run GLM
##  This calls the GLM script.  Class input files from KNIFE dir/circReads/ids/ are fed into the GLM and GLM reports are generated in FJDir/reports/glmReports.  Please see GLM documentation for additional information.
#j15b_id
start_time = time.time()
print("STATUS:Run GLM")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_15GLM_r.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_15GLM_r.txt"),"w")
    cmd = "{MACHETE_DIR}/run_GLM.sh {CIRCPIPE_DIR} {OUTPUT_DIR} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Run glm "+str(index),start_time,timer_file_path)

# Commenting this out for now- June 8 2016; had an error of this type:
#    circInputFile= open( glob.glob(glmDir+"*"+args.stem+"*circJuncProbs.txt")[0], mode="rU")
#    IndexError: list index out of range
# ## Append linear junctions GLM report with anomalies, indels
# # AddIndelstolinearGLM.sh calls the python script
# # KNIFEglmReportsForMachete.py.  This script parses circular and linear
# # glmReports.  For the linear glmReports from KNIFE, the script collects
# # any junctions where 1) the two exons from the linear report are from
# # different genes or 2) the posterior probability is >0.9.  It adds on
# # the rate of anomaly reads and indels to the reports and feeds them
# # into FJDir/reports/AppendedReports.  For ciruclar reports, the script
# # collects any junctions where the posterior probability is <0.9,
# # appends the "Decoy" rate, and feeds the reports into
# # FJDir/reports/Appended reports.
# ## The purpose of this script is to place all reports in a single directory for the user.
# #j17_id
# print("Appending linearJuncs GLM report")
# processes = {}
# for index in range(1,NUM_FILES + 1):
#   stdout = open(os.path.join(LOG_DIR,str(index) + "_out_17AppendRegGLM.txt"),"w")
#   stderr = open(os.path.join(LOG_DIR,str(index) + "_err_17AppendGLM.txt"),"w")
#   cmd = "{MACHETE_DIR}/AddIndelstolinearGLM.sh {CIRCPIPE_DIR} {OUTPUT_DIR} {MACHETE_DIR} {index}".format(MACHETE_DIR=MACHETE_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
#   popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
#   processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
# checkProcesses(processes)
##
####

# The AppendNaiveRept.sh shell calls the AppendNaiveRept.py script.  This reads in the IndelsHistogram, BadFJ and BadFJ_ver2 files, and GLM report results and outputs all the results into a single file in /FJDir/reports/AppendedReports/<STEM>_naive_report_Appended.txt
# j15_id
start_time = time.time()
print("STATUS: append naive rpt")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_14AppendRpt.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_14AppendRpt.txt"),"w")
    cmd = "{MACHETE_DIR}/AppendNaiveRept.sh {OUTPUT_DIR} {GLM_DIR} {MACHETE_DIR} {OUTPUT_DIR}/reports/glmReports {index}".format(MACHETE_DIR=MACHETE_DIR,OUTPUT_DIR=OUTPUT_DIR,GLM_DIR=GLM_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
write_time("Append naive report "+str(index),start_time,timer_file_path)
print("STATUS: done")

write_time("Full run time",full_time,timer_file_path)




