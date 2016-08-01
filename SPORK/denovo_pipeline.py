#################
# Documentation #
#################
#Last update : 7/13/16
#Script name : SPORK (Small Powerful Orthogonal Read Mapper from KNIFE)
#Supported by: Salzman Lab Stanford University <http://med.stanford.edu/salzmanlab.html>
#Algorithm by: Julia Salzman
#Written by  : Rob Bierman <rbierman@stanford.edu>
#Relies on   : KNIFE <https://github.com/lindaszabo/KNIFE>
#Python      : 2.7.5
#Wrapper     : B2_run_denovo_scripts.batch
#Example call: python denovo_pipeline.py $ALIGN_PARDIR $DATASET_NAME $MODE $NUM_FLANKING $NTRIM $DENOVOCIRC 1> "SPORK_out.log" 2> "SPORK_out.err"

###################
# General Imports #
###################
import subprocess
import time
import sys
import os

###########################
# Script specific Imports #
###########################
self_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(self_path+"/denovo_class_utils/")
from denovo_consensus_utils import *
from denovo_utils import *
from Junction import Junction
from BinPair import BinPair
from GTFEntry import GTFEntry
from SAMEntry import SAMEntry
from FastQEntry import FastQEntry

##########################
# Command line arguments #
##########################
# Arguments from denovo_scripts.bash
# sys.argv[0] is this program's name
# sys.argv[1] is the align_parent_directory
# sys.argv[2] is the dataset name
# sys.argv[3] is the mode
# sys.argv[4] is the seq len on each side of the splice site in denovo jcts
# sys.argv[5] is the ntrim
# sys.argv[6] is the circle prefs (0 is linear and circular, 1 is circular only)
# sys.argv[7] is optionally the output directory (for SPACHETE)
# sys.argv[8] is optionally the stem name (for SPACHETE)

# NOTE currently requires all args and has no defaults
# Checking for the correct number of arguments
if len(sys.argv) < 7:
    sys.stderr.write("Not enough arguments passed to main denovo python file (need 6)")
    sys.stderr.write("\tNumber of arguments passed: "+str(len(sys.argv)-1))
    sys.stderr.write("\tExample call: python denovo_pipeline.py $ALIGN_PARDIR $DATASET_NAME $MODE $NUM_FLANKING $NTRIM $DENOVOCIRC")
    sys.exit(1)

# Assigning arguments to their correct variables
parent_dir = sys.argv[1]
dataset_name = sys.argv[2]
mode = sys.argv[3]                      #NOTE only allows hg19 use and any mode will default to hg19
splice_flank_len = sys.argv[4]
ntrim = sys.argv[5]                     #NOTE not currently used
circle_prfs = sys.argv[6]               #NOTE not currently used
out_dir = "/scratch/PI/horence/rob/spork_outputs/"+dataset_name+"/"
if len(sys.argv) >= 8:
    out_dir = sys.argv[7]+"/spork_out/"

stem_name = None
if len(sys.argv) >= 9:
    stem_name = sys.argv[8]

#############
# Constants #
#############
bin_size = 50                           #size of bins in bps to split ref into
group_member_cutoff = 2                 #minimum number of reads that need to map to a bin pair
consensus_score_cutoff = 0.5            #relates to the number of mismatches in consensus
min_bases_per_col = 1                   #only consider column if it has at least n bases
read_gap_score = "--rdg 50,50"          #read gap set to be very high to not allow gaps
ref_gap_score = "--rfg 50,50"           #read gap set to be very high to not allow gaps
min_score = "--score-min L,0,-0.10"     #minimum allowed Bowtie2 score
allowed_mappings = "1"                  #allowed mappings of a given read. Currently not implemented

splice_finding_min_score = min_score    #can make more stringent for splice site finding
splice_finding_allowed_mappings = "1"   #just making this adjustable, probably want to stay at 1
splice_finding_allowed_mismatches = "1" #used to find high quality splice sites
num_threads = "8"                       #threads to use in each Bowtie2 call
use_prior = True                        #if set to true will use previously made files if available
timer_file_path = out_dir+"timed_events.txt"
start_entire_time = time.time()

########################################
# Searching Parent Dir for Input Files #
########################################
# Find the correct reference index based on the mode
reference, gtf = get_reference_and_gtf_from_mode(mode)

# Loop through the directory to use only R1 fastq files
# NOTE maybe move this to the utils script
unaligned_path = parent_dir + dataset_name + "/orig/unaligned/forDenovoIndex"
all_file_names = os.listdir(unaligned_path)
file_names = []
R1_patterns = []
R2_patterns = []
for file_name in all_file_names:
    if ".fq" in file_name or ".fastq" in file_name:
        # Pass the R2 files. Helps if a file looks like: "Human_simulated_reads1_2.fq"
        if "R2" in file_name or "r2" in file_name or "_2.fq" in file_name:
            continue
        if "_1.fq" in file_name and (not stem_name or stem_name in file_name):
            R1_patterns.append("_1.fq")
            R2_patterns.append("_2.fq")
            file_names.append(file_name)
        elif "_1.fastq" in file_name and (not stem_name or stem_name in file_name):
            R1_patterns.append("_1.fastq")
            R2_patterns.append("_2.fastq")
            file_names.append(file_name)
        elif "R1" in file_name and (not stem_name or stem_name in file_name):
            R1_patterns.append("R1")
            R2_patterns.append("R2")
            file_names.append(file_name)
        elif "r1" in file_name and (not stem_name or stem_name in file_name):
            R1_patterns.append("r1")
            R2_patterns.append("r2")
            file_names.append(file_name)
        elif "read1" in file_name and (not stem_name or stem_name in file_name):
            R1_patterns.append("read1")
            R2_patterns.append("read2")
            file_names.append(file_name)
        elif "reads1" in file_name and (not stem_name or stem_name in file_name):
            R1_patterns.append("reads1")
            R2_patterns.append("reads2")
            file_names.append(file_name)

# Making a dict that stores the input and constant values to make argument passing easier
# This is just a hash that stores the constants set above to allow them easy to pass
constants_dict = {"parent_dir":parent_dir,"dataset_name":dataset_name,"mode":mode,"splice_flank_len":splice_flank_len,
                  "ntrim":ntrim,"circle_prfs":circle_prfs,"bin_size":bin_size,"group_member_cutoff":group_member_cutoff,
                  "consensus_score_cutoff":consensus_score_cutoff,"min_score":min_score,"read_gap_score":read_gap_score,
                  "splice_finding_min_score":splice_finding_min_score,"read_gap_score":read_gap_score,"min_bases_per_col":min_bases_per_col,
                  "splice_finding_allowed_mismatches":splice_finding_allowed_mismatches,"unaligned_path":unaligned_path,
                  "splice_finding_allowed_mappings":splice_finding_allowed_mappings,"ref_gap_score":ref_gap_score,
                  "allowed_mappings":allowed_mappings,"num_threads":num_threads,"reference":reference,"gtf":gtf}

###################################
# Loop through each R1 input file #
###################################
for file_name in file_names:
    # Create out_dir specific for each file_name
    start_full_file_time = time.time()
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    #out_dir += "spork_"+file_name.split(".")[0]+"/"
    #if not os.path.isdir(out_dir):
    #    os.mkdir(out_dir)
    constants_dict["out_dir"] = out_dir
    
    # Get the full path of the input fq reads file
    # Also setup the timer file to see how long each portion of SPORK takes
    full_path = unaligned_path+"/"+file_name
    timer_file_path = out_dir+"timed_events.txt"
    constants_dict["timer_file_path"] = timer_file_path
    genome_sam_dir = parent_dir+dataset_name+"/orig/genome/"
    constants_dict["genome_sam_dir"] = genome_sam_dir

    # Identify whether there is an R2 present in the same directory as the R1
    start_find_R2 = time.time()
    R1_pattern = R1_patterns[file_names.index(file_name)]
    R2_pattern = R2_patterns[file_names.index(file_name)]
    base,extension = full_path.split(R1_pattern)
    R1_file = full_path
    R2_file = base+R2_pattern+extension
    full_paths_in = [R1_file]
    if not os.path.isfile(R2_file):
        R2_file = None
    if R2_file:
        print "Found R2 file:",R2_file
        full_paths_in.append(R2_file)
    else:
        print "Not found R2 file: Proceeding with just R1"
    sys.stdout.flush()
    write_time("Time to find R2 "+R2_file.split("/")[-1],start_find_R2,timer_file_path,False) #false overwrites the timer file

    # Create file locations for fq files w/out spaces in headers
    R1_file_name = file_name
    R1_file_path = out_dir+R1_file_name
    R2_file_name = R1_file_name.split(R1_pattern)[0]+R2_pattern+R1_file_name.split(R1_pattern)[1]
    R2_file_path = out_dir+R2_file_name
    full_paths_out = [R1_file_path,R2_file_path]
    R1_R2_names = [R1_file_name,R2_file_name]

    # Remove spaces from the header lines of the fq file (also remove the last split ind)
    # NOTE creates a no_spaces version of each file removing spaces in header lines
    # NOTE also removes the last "/" and everything to the right
    # NOTE this should be moved to denovo utils
    for file_ind in range(len(full_paths_in)):
        file_path_in = full_paths_in[file_ind]
        file_path_out = full_paths_out[file_ind]
        f_name = file_path_in.split("/")[-1]
        start_remove_spaces = time.time()
        with open(file_path_in,"r") as fq_read_file:
            with open(file_path_out,"w") as fq_write_file:
                fq_line = fq_read_file.readline()
                num_reads = 0
                while fq_line:
                    # Remove spaces and things to the right of the last space
                    num_reads += 1
                    split_fq_line = fq_line.split(" ")
                    if len(split_fq_line) > 1:
                        fq_line = "_".join(split_fq_line[:-1])+"\n"
                    else:
                        fq_line = split_fq_line[0]

                    # Remove forward slashes and things to to right of the last "/"
                    split_fq_line = fq_line.split("/")
                    if len(split_fq_line) > 1:
                        fq_line = "-".join(split_fq_line[:-1])+"\n"
                    else:
                        fq_line = split_fq_line[0]

                    fq_write_file.write(fq_line)
                    fq_write_file.write(fq_read_file.readline()) #Skip the seq line
                    fq_write_file.write(fq_read_file.readline()) #Skip the + line
                    fq_write_file.write(fq_read_file.readline()) #Skip the quality line
                    fq_line = fq_read_file.readline() #Get to the next header line
        write_time("Time to remove spaces "+f_name,start_remove_spaces,timer_file_path)
        write_time("--> Had "+str(num_reads)+" reads",start_remove_spaces,timer_file_path)

    #Do the SPORK main code for both R1 and R2 if they are both available
    for r_ind in range(len(full_paths_out)):
        # Process the file to split each read into a 5' and 3' fastq file
        write_time("Starting main portion R"+str(r_ind+1),time.time(),timer_file_path)
        start_split_reads = time.time()
        R_file_path = full_paths_out[r_ind]
        R_file_name = R1_R2_names[r_ind]
        five_prime_fq_name = out_dir + "R"+str(r_ind+1)+"_5prime_" + R_file_name.split(".")[0] + ".fq"
        three_prime_fq_name = out_dir + "R"+str(r_ind+1)+"_3prime_" + R_file_name.split(".")[0] + ".fq"
        if use_prior and os.path.isfile(five_prime_fq_name) and os.path.isfile(three_prime_fq_name):
            write_time("Using previous split unaligned read files "+R_file_name,start_split_reads,timer_file_path)
        else:
            with open(five_prime_fq_name, "w") as five_prime_file:
                with open(three_prime_fq_name, "w") as three_prime_file:
                    #NOTE THIS is where only R1 is being used instead of both R1 and R2
                    with open(R_file_path, "r") as f_in:
                        read_id = f_in.readline()
                        while read_id:
                            seq = f_in.readline()
                            plus_line = f_in.readline()
                            quality = f_in.readline()
                            fastq_read = FastQEntry(read_id, seq, plus_line, quality)
                            fragment_5, fragment_3 = fastq_read.get_edge_thirds()
                            five_prime_file.write(str(fragment_5))
                            three_prime_file.write(str(fragment_3))
                            read_id = f_in.readline()
            write_time("Time to make split unaligned read files "+R_file_name,start_split_reads,timer_file_path)

        # Map the 5' and 3' split files to the reference to generate the sam files
        # NOTE this can take a very long time for large R1 files
        start_split_mapping = time.time()
        five_prime_mapped_name = five_prime_fq_name.split(".")[0] + ".sam"
        three_prime_mapped_name = three_prime_fq_name.split(".")[0] + ".sam"
        if use_prior and os.path.isfile(five_prime_mapped_name) and os.path.isfile(three_prime_mapped_name):
            write_time("Using previous split unaligned reads "+R_file_name,start_split_mapping,timer_file_path)
        else:
            five_prime_mapped = open(five_prime_mapped_name, "w")
            three_prime_mapped = open(three_prime_mapped_name, "w")
            start_split_mapping = time.time()
            subprocess.call(
                ["bowtie2", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
                 "--un", out_dir + "R"+str(r_ind+1)+"_unaligned_5prime_" + R_file_name, "-x", reference, five_prime_fq_name], stdout=five_prime_mapped)
            write_time("Time to map split unaligned reads 5'"+R_file_name,start_split_mapping,timer_file_path)
            start_split_mapping = time.time()
            subprocess.call(
                ["bowtie2", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
                 "--un", out_dir + "R"+str(r_ind+1)+"_unaligned_3prime_" + R_file_name, "-x", reference, three_prime_fq_name], stdout=three_prime_mapped)
            write_time("Time to map split unaligned reads 3'"+R_file_name,start_split_mapping,timer_file_path)
            five_prime_mapped.close()
            three_prime_mapped.close()

        ##########################
        # Build Denovo Jct Fasta #
        ##########################
        # Make the fasta file to build the bowtie index from in the following steps:
        # [1] Store all 5' mappings in a dict by shared read id
        # [2] Walk through the 3' mappings and add them and their 5' counterpart to a bin_pair list
        # [3] Find the group ranges within the bin_pair list
        # [4] Build the denovo jcts
        # [5] Collapse the denovo jcts
        # [6] Get GTF annotation data for each jct
        # [7] Write out the denovo jcts
        fasta_for_bowtie_index_name = out_dir+"R"+str(r_ind+1)+"_novel_junctions.fasta"
        gtfs = generate_gtfs(gtf)

        if use_prior and os.path.isfile(fasta_for_bowtie_index_name):
            write_time("Using prior jcts: "+fasta_for_bowtie_index_name,time.time(),timer_file_path)
        else:
            # Store all five prime mappings by base_read_id (read_id w/out 5_prime or 3_prime)
            # There should not be two identical base_read_id's
            id_to_sam_dict = {}
            with open(five_prime_mapped_name,"r") as five_prime_mapped:
                sam_line = five_prime_mapped.readline()
                while sam_line and "@" == sam_line[0]: #Read past the header lines
                    sam_line = five_prime_mapped.readline()

                while sam_line:
                    sam_entry = SAMEntry(sam_line)
                    base_read_id = sam_entry.read_id.split("/")[0]
                    if base_read_id in id_to_sam_dict:
                        sys.stderr.write("ERROR: Found duplicate base_read_id in 5_prime mappings\n")
                        sys.exit(1)
                    # Filter out the strange chromosomes: (e.g. chrUn_gl000220)
                    if "_" not in sam_entry.chromosome:
                        id_to_sam_dict[base_read_id] = sam_entry
                    sam_line = five_prime_mapped.readline()
                    
            # Now walk through the three prime mappings creating bin pairs from all shared ids
            bin_pairs = []
            with open(three_prime_mapped_name,"r") as three_prime_mapped:
                sam_line = three_prime_mapped.readline()
                while sam_line and "@" == sam_line[0]: #Read past the header lines
                    sam_line = three_prime_mapped.readline()

                while sam_line:
                    three_prime_sam = SAMEntry(sam_line)
                    base_read_id = three_prime_sam.read_id.split("/")[0]
                    # Filter out the strange chromosomes: (e.g. chrUn_gl000220)
                    if "_" not in three_prime_sam.chromosome and base_read_id in id_to_sam_dict:
                        five_prime_sam = id_to_sam_dict[base_read_id]
                        five_prime_bin = five_prime_sam.start/bin_size
                        three_prime_bin = three_prime_sam.start/bin_size
                        bin_pair = BinPair(five_prime_sam,three_prime_sam,five_prime_bin,three_prime_bin)
                        bin_pairs.append(bin_pair)
                    sam_line = three_prime_mapped.readline()
            id_to_sam_dict = {} #clearing the dictionary to free up space

            # Sort the bin_pairs by bin_pair id to form list w/ groups adjacent
            bin_pairs.sort()

            # Save the bin_pairs to a file to see how they look
            bin_pair_out_file_name = out_dir+"R"+str(r_ind+1)+"_bin_pairs.txt"
            with open(bin_pair_out_file_name, "w") as f:
                for bin_pair in bin_pairs:
                    f.write(str(bin_pair))

            # Group the bin pairs by finding their index ranges within bin_pairs
            start_build_group_ranges = time.time()
            bin_pair_group_ranges = find_bin_pair_group_ranges(bin_pairs,constants_dict)
            write_time("-Time to build group ranges",start_build_group_ranges,timer_file_path)

            # Build the junctions from bin pairs and their ends
            denovo_junctions = []
            start_build_junctions = time.time()
            denovo_junctions = build_junction_sequences(bin_pairs,bin_pair_group_ranges,R_file_path,constants_dict)
            write_time("-Time to build junctions ",start_build_junctions,timer_file_path)
            bin_pairs = [] #clearing the bin pairs to free up space

            # Find the splice indicies of the junctions
            start_find_splice_inds = time.time()
            denovo_junctions,no_splice_jcts = find_splice_inds(denovo_junctions,constants_dict,r_ind)
            write_time("-Time to find splice indicies ",start_find_splice_inds,timer_file_path)



            # Collapse the junctions that have the same splice site
            # TODO I can't figure out how to get jct collapsing working so I'm just removing it for now
            # TODO there is some sort of offset error
            """
            start_collapse_junctions = time.time()
            denovo_junctions = collapse_junctions(denovo_junctions)
            write_time("-Time to collapse junctions ",start_collapse_junctions,timer_file_path)
            """

            # Get GTF information for the identified denovo_junctions
            start_get_jct_gtf_info = time.time()
            get_jct_gtf_info(denovo_junctions,gtfs)
            write_time("Time to get jct gtf info "+R_file_name,start_get_jct_gtf_info,timer_file_path)
            print len(denovo_junctions)

            # Identify fusions from the junctions
            start_identify_fusions = time.time()
            fusion_junctions = identify_fusions(denovo_junctions)
            sys.stderr.write("Len fusion jcts = "+str(len(fusion_junctions))+"\n")
            write_time("Time to identify fusions "+R_file_name,start_identify_fusions,timer_file_path)

            #########################################################
            # Write out the denovo_junction_sequences for each file #
            #########################################################
            start_save_denovo_junctions = time.time()
            machete_style_name = out_dir+"R"+str(r_ind+1)+"_novel_junctions_machete.fasta"
            jct_style_file_name = out_dir+"R"+str(r_ind+1)+"_novel_junctions.jct"
            fusions_file_name = out_dir+"R"+str(r_ind+1)+"_novel_fusions.fasta"

            # Open the three output files
            fasta_for_bowtie_index = open(fasta_for_bowtie_index_name, "w")
            machete_style_file = open(machete_style_name, "w")
            jct_style_file = open(jct_style_file_name, "w")
            fusions_file = open(fusions_file_name, "w")

            # Loop through the denovo junctions writing them where necessary
            for denovo_junction in denovo_junctions:
                jct_ind = denovo_junctions.index(denovo_junction)
                fasta_for_bowtie_index.write(denovo_junction.verbose_fasta_string())
                machete_style_file.write(denovo_junction.fasta_MACHETE())
                jct_style_file.write(str(denovo_junction)+"\n")
                if denovo_junction in fusion_junctions:
                    fusions_file.write(denovo_junction.fasta_MACHETE())

            # Close the three output files
            fasta_for_bowtie_index.close()
            machete_style_file.close()
            jct_style_file.close()
            fusions_file.close()
            write_time("Time to save denovo junctions "+R_file_name,start_save_denovo_junctions,timer_file_path)
            denovo_junctions = [] #NOTE clearing the denovo junctions to free up space
        
    #Combine the two machete style fastas for SPACHETE
    full_machete_style_name = out_dir+"novel_junctions_machete.fasta"
    with open(full_machete_style_name,"w") as full_machete_style:
        R1_machete_style_name = out_dir+"R1_novel_junctions_machete.fasta"
        if os.path.isfile(R1_machete_style_name):
            for R1_line in open(R1_machete_style_name,"r"):
                full_machete_style.write(R1_line)

        R2_machete_style_name = out_dir+"R2_novel_junctions_machete.fasta"
        if os.path.isfile(R2_machete_style_name):
            for R2_line in open(R2_machete_style_name,"r"):
                full_machete_style.write(R2_line)

           
write_time("Entire time",start_entire_time,timer_file_path)


