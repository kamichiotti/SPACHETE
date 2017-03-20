import sys
import os

def make_machete_output_dirs(OUTPUT_DIR):

    #MACHETE output directories
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
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
    GLM_classInputDir = os.path.join(OUTPUT_DIR,"GLM_classInput")
    if not os.path.exists(GLM_classInputDir):
        os.mkdir(GLM_classInputDir)
        


if __name__ == "__main__":
    OUTPUT_DIR = "/home/rbierman/test_spachete"
    CIRCPIPE_DIR = "/home/rbierman/test_spachete"
    make_output_dir(OUTPUT_DIR+"/MACHETE")
    print "done"
