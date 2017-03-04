from argparse import ArgumentParser
import subprocess
import sys
import os

#MACHETE = os.path.dirname(__file__)
MACHETE = os.path.dirname(os.path.abspath(__file__))
print "path:",MACHETE

#Check the processes for good output
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


# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for

description = "Required args: --circpipe-dir, --output-dir, --hg19Exons, --reg-indel-indices, and --circref-dir."

parser = ArgumentParser(description=description)
parser.add_argument("--circpipe-dir",required=True,dest="CIRCPIPE_DIR",help="Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)")
parser.add_argument("--output-dir",required=True,dest="OUTPUT_DIR",help="Output directory for resuts. Directory path will be created recursively if it doesn't already exist. If it exists already, the directory will be deleted then created again.")

args = parser.parse_args()

CIRCPIPE_DIR = args.CIRCPIPE_DIR
OUTPUT_DIR = args.OUTPUT_DIR
LOG_DIR = os.path.join(OUTPUT_DIR,"glm_rerun_logs")
if not os.path.isdir(LOG_DIR):
    os.makedirs(LOG_DIR)

StemFile = os.path.join(OUTPUT_DIR,"StemList.txt")
NUM_FILES = len(open(StemFile,"r").readlines())
GLM_DIR = os.path.join(CIRCPIPE_DIR,"circReads/glmReports")
#end arg parsing

###### RUN GLM ###########################
#
## Run GLM
##  This calls the GLM script.  Class input files from KNIFE dir/circReads/ids/ are fed into the GLM and GLM reports are generated in FJDir/reports/glmReports.  Please see GLM documentation for additional information.
#j15b_id
print("STATUS:Run GLM")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_run_GLM.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_run_GLM.txt"),"w")
    cmd = "{MACHETE}/run_GLM.sh {CIRCPIPE_DIR} {OUTPUT_DIR} {MACHETE} {index}".format(MACHETE=MACHETE,CIRCPIPE_DIR=CIRCPIPE_DIR,OUTPUT_DIR=OUTPUT_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#Appending to the naive report
print("STATUS: append naive rpt")
processes = {}
for index in range(1,NUM_FILES + 1):
    stdout = open(os.path.join(LOG_DIR,str(index) + "_out_AppendRpt.txt"),"w")
    stderr = open(os.path.join(LOG_DIR,str(index) + "_err_AppendRpt.txt"),"w")
    cmd = "{MACHETE}/AppendNaiveRept.sh {OUTPUT_DIR} {GLM_DIR} {MACHETE} {OUTPUT_DIR}/reports/glmReports {index}".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,GLM_DIR=GLM_DIR,index=index)
    print(cmd)
    popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
    processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
print("STATUS: done")





