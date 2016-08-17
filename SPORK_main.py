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
import argparse
import time
import sys
import os

#Get the absolute path to this file for 7B
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

###########################
# Script specific Imports #
###########################
from SPORK_consensus_utils import *
from SPORK_utils import *
from SPORK_Junction import Junction
from SPORK_BinPair import BinPair
from SPORK_GTFEntry import GTFEntry
from SPORK_SAMEntry import SAMEntry
from SPORK_FastQEntry import FastQEntry


###########################
#  Parse Input Arguments  #
###########################
arg_description = "Required Args: --input-dir, --output-dir, --ref-dir"
arg_description += "\n"+"Optional Args: --stem-name, --flank-len, --mode"

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument("--input-dir",required=True,dest="input_dir",help="Directory containing KNIFE produced orig directory " \
                    "(i.e. '/scratch/PI/horence/rob/parent_dirs/fetal_lungs/fetal_lung_01')")
parser.add_argument("--output-dir",required=True,dest="output_dir",help="Output directory for resuts " \
                    "(i.e. '/scratch/PI/horence/rob/spachete_outputs/normal_fetal')")
parser.add_argument("--ref-dir",required=True,dest="ref_dir",help="Directory with bowtie2 indices of the desired genome, same as " \
                    "MACHETE CIRCREF (i.e. '/share/PI/horence/circularRNApipeline_Cluster/index')")
parser.add_argument("--stem-name",required=False,dest="stem_name",help="Stem name for SPORK use in SPACHETE " \
                    "to restrict SPORK from working on all the R1 files (i.e. 'SRR3192409')")
parser.add_argument("--flank-len",required=False,dest="flank_len",help="Flanking length on either side of the junction " \
                    "splice site. For example a flank len of 150 results in 300bp long junctions reported.",type=int)
parser.add_argument("--mode",required=False,dest="mode",help="Mode is the reference to use to align to, note that the proper " \
                    "bowtie2 indices and gtf file must exist in the indices folder. Defaults to hg19 (human genome)")

args = parser.parse_args()
input_dir = args.input_dir
output_dir = os.path.join(args.output_dir,"spork_out")
ref_dir = args.ref_dir
stem_name = args.stem_name
splice_flank_len = args.flank_len
mode = args.mode

#Check if there is an orig in the right spot
orig_dir = os.path.join(input_dir,"orig")
if not os.path.exists(orig_dir):
    sys.stderr.write("ERROR: No orig directory found in --input-dir\n")
    sys.exit(1)

#Set defaults if not given
if not splice_flank_len:
    splice_flank_len = 150
if not mode:
    mode = "hg19"

#############
# Constants #
#############
bin_size = 50                           #size of bins in bps to split ref into
group_member_cutoff = 2                 #minimum number of reads that need to map to a bin pair
consensus_score_cutoff = 0.5            #relates to the number of mismatches in consensus
span_cutoff = 1e6                       #minimum span distance to be classified a "fusion"
at_boundary_cutoff = 15                 #maximum distance to be away from a boundary for fusion classification
collapse_thresh = 5                     #distance used in the collapse junctions step
min_bases_per_col = 2                   #only consider column if it has at least n bases
read_gap_score = "--rdg 50,50"          #read gap set to be very high to not allow gaps
ref_gap_score = "--rfg 50,50"           #read gap set to be very high to not allow gaps
min_score = "--score-min L,0,-0.10"     #minimum allowed Bowtie2 score
allowed_mappings = "1"                  #allowed mappings of a given read. Currently not implemented

thirds_len = 36                         #length to cut the original reads into for the 5' and 3' pieces:
                                        #   |-------|---------------|-------|
                                        #       ^     unused middle     ^
                                        #       |                       |
                                        # len(5' piece) == thirds_len   |   
                                        #                          len(3' piece) == third_len

splice_finding_min_score = min_score    #can make more stringent for splice site finding
splice_finding_allowed_mappings = "1"   #just making this adjustable, probably want to stay at 1
splice_finding_allowed_mismatches = "1" #used to find high quality splice sites
num_threads = "8"                       #threads to use in each Bowtie2 call
use_prior = True                        #if set to true will use previously made files if available
timer_file_path = os.path.join(output_dir,"timed_events.txt")
start_entire_time = time.time()

########################################
# Searching Parent Dir for Input Files #
########################################
# Find the correct reference index based on the mode
reference,gtf_path = get_reference_and_gtf_from_mode(ref_dir,ABS_PATH,mode)

# Loop through the directory to use only R1 fastq files
# NOTE maybe move this to the utils script
unaligned_path = os.path.join(orig_dir,"unaligned","forDenovoIndex")
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
# Probably sloppy code
constants_dict = {"input_dir":input_dir,"mode":mode,"splice_flank_len":splice_flank_len,"bin_size":bin_size,"group_member_cutoff":group_member_cutoff,
                  "consensus_score_cutoff":consensus_score_cutoff,"min_score":min_score,"read_gap_score":read_gap_score,"thirds_len":thirds_len,
                  "splice_finding_min_score":splice_finding_min_score,"read_gap_score":read_gap_score,"min_bases_per_col":min_bases_per_col,
                  "splice_finding_allowed_mismatches":splice_finding_allowed_mismatches,"unaligned_path":unaligned_path,
                  "splice_finding_allowed_mappings":splice_finding_allowed_mappings,"ref_gap_score":ref_gap_score,"use_prior":use_prior,
                  "allowed_mappings":allowed_mappings,"num_threads":num_threads,"reference":reference,"gtf_path":gtf_path,"collapse_thresh":collapse_thresh,
                  "at_boundary_cutoff":at_boundary_cutoff,"span_cutoff":span_cutoff}

###################################
# Loop through each R1 input file #
###################################
for file_name in file_names:
    # Create output_dir specific for each file_name
    start_full_file_time = time.time()
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    #output_dir += "spork_"+file_name.split(".")[0]+"/"
    #if not os.path.isdir(output_dir):
    #    os.mkdir(output_dir)
    constants_dict["output_dir"] = output_dir
    
    # Get the full path of the input fq reads file
    # Also setup the timer file to see how long each portion of SPORK takes
    full_path = os.path.join(unaligned_path,file_name)
    timer_file_path = os.path.join(output_dir,"timed_events.txt")
    constants_dict["timer_file_path"] = timer_file_path
    genome_sam_dir = os.path.join(orig_dir,"genome/")
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
        sys.stdout.write("Not found R2 file: Proceeding with just R1\n")
    else:
        sys.stdout.write("Found R2 file :"+R2_file+"\n")
        full_paths_in.append(R2_file)
    write_time("Time to find R2 "+R2_file.split("/")[-1],start_find_R2,timer_file_path,False) #false overwrites the timer file

    # Create file locations for fq files w/out spaces in headers
    R1_file_name = file_name
    R1_file_path = os.path.join(output_dir,R1_file_name)
    R2_file_name = R1_file_name.split(R1_pattern)[0]+R2_pattern+R1_file_name.split(R1_pattern)[1]
    R2_file_path = os.path.join(output_dir,R2_file_name)
    full_paths_out = [R1_file_path,R2_file_path]
    R1_R2_names = [R1_file_name,R2_file_name]

    #Combine the R1 and R2 fastqs adding R1 and R1 to the headers
    #to disambiguate in case of identical matching
    if stem_name:
        combined_out_name = os.path.join(output_dir,stem_name+"_combined_reads.fq")
    else:
        combined_out_name = os.path.join(output_dir,"combined_reads.fq")
    with open(combined_out_name,"w") as combined_out:
        if len(full_paths_in) >= 1:
            with open(full_paths_in[0],"r") as R1_in_file:
                header_line = R1_in_file.readline()
                while header_line:
                    combined_out.write(header_line.strip()+" R1\n")
                    combined_out.write(R1_in_file.readline()) #Write out the seq line
                    combined_out.write(R1_in_file.readline()) #Write out the + line
                    combined_out.write(R1_in_file.readline()) #Write out the quality line
                    header_line = R1_in_file.readline()

        if len(full_paths_in) >= 2:
            with open(full_paths_in[1],"r") as R2_in_file:
                header_line = R2_in_file.readline()
                while header_line:
                    combined_out.write(header_line.strip()+" R2\n")
                    combined_out.write(R2_in_file.readline()) #Write out the seq line
                    combined_out.write(R2_in_file.readline()) #Write out the + line
                    combined_out.write(R2_in_file.readline()) #Write out the quality line
                    header_line = R2_in_file.readline()


    # Process the file to split each read into a 5' and 3' fastq file
    write_time("Starting main portion",time.time(),timer_file_path)
    start_split_reads = time.time()
    R_file_path = combined_out_name
    R_file_name = stem_name if stem_name else "combined"
    five_prime_fq_name = os.path.join(output_dir,"5prime_"+R_file_name.split(".")[0]+".fq")
    three_prime_fq_name = os.path.join(output_dir,"3prime_"+R_file_name.split(".")[0]+".fq")
    if use_prior and os.path.isfile(five_prime_fq_name) and os.path.isfile(three_prime_fq_name):
        write_time("Using previous split unaligned read files "+R_file_name,start_split_reads,timer_file_path)
    else:
        with open(five_prime_fq_name, "w") as five_prime_file:
            with open(three_prime_fq_name, "w") as three_prime_file:
                with open(R_file_path, "r") as f_in:
                    read_id = f_in.readline()
                    while read_id:
                        seq = f_in.readline()
                        plus_line = f_in.readline()
                        quality = f_in.readline()
                        fastq_read = FastQEntry(read_id, seq, plus_line, quality)
                        fragment_5, fragment_3 = fastq_read.get_first_last_n(thirds_len)
                        #Check that both fragments are defined before writing them out
                        if fragment_5 and fragment_3:
                            five_prime_file.write(str(fragment_5))
                            three_prime_file.write(str(fragment_3))
                        read_id = f_in.readline()
        write_time("Time to make split unaligned read files "+R_file_name,start_split_reads,timer_file_path)

    # Map the 5' and 3' split files to the reference to generate the sam files
    # NOTE this can take a very long time for large R1 files
    start_split_mapping = time.time()
    five_prime_mapped_name = five_prime_fq_name.split(".")[0]+".sam"
    three_prime_mapped_name = three_prime_fq_name.split(".")[0]+".sam"
    if use_prior and os.path.isfile(five_prime_mapped_name) and os.path.isfile(three_prime_mapped_name):
        write_time("Using previous split unaligned reads "+R_file_name,start_split_mapping,timer_file_path)
    else:
        five_prime_mapped = open(five_prime_mapped_name, "w")
        three_prime_mapped = open(three_prime_mapped_name, "w")
        start_split_mapping = time.time()
        subprocess.call(
            ["bowtie2", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
             "--un", os.path.join(output_dir,"unaligned_5prime_"+R_file_name+".fq"), "-x", reference, five_prime_fq_name], stdout=five_prime_mapped)
        write_time("Time to map split unaligned reads 5'"+R_file_name,start_split_mapping,timer_file_path)
        start_split_mapping = time.time()
        subprocess.call(
            ["bowtie2", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
             "--un", os.path.join(output_dir,"unaligned_3prime_"+R_file_name+".fq"), "-x", reference, three_prime_fq_name], stdout=three_prime_mapped)
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
    fasta_for_bowtie_index_name = os.path.join(output_dir,"novel_junctions.fasta")
    gtfs = generate_gtfs(gtf_path)

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
                base_read_id = sam_entry.read_id.replace("/5_prime","")
                if base_read_id in id_to_sam_dict:
                    sys.stderr.write("ERROR: Found duplicate base_read_id in 5_prime mappings\n")
                    sys.stderr.write(base_read_id+"\n")
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
                base_read_id = three_prime_sam.read_id.replace("/3_prime","")
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
        bin_pair_out_file_name = os.path.join(output_dir,"bin_pairs.txt")
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

        # Write out the ids of the reads used to build junctions
        used_read_ids_name = os.path.join(output_dir,"used_read_ids.txt")
        with open(used_read_ids_name,"w") as used_read_ids:
            for jct in denovo_junctions:
                for bin_pair in jct.bin_pair_group:
                    used_read_id = bin_pair.five_prime_SAM.read_id
                    if " R1" in used_read_id:
                        trimmed_read_id = used_read_id.replace(" R1","")
                    if " R2" in used_read_id:
                        trimmed_read_id = used_read_id.replace(" R2","")                   
                    used_read_ids.write("@"+trimmed_read_id+"\n")

        # Find the splice indicies of the junctions
        start_find_splice_inds = time.time()
        denovo_junctions,no_splice_jcts = find_splice_inds(denovo_junctions,constants_dict)
        write_time("-Time to find splice indices ",start_find_splice_inds,timer_file_path)

        # Write out the denovo_junctions before collapsing for debugging
        start_write_pre_collapsed = time.time()
        machete_style_name = os.path.join(output_dir,"pre_collapse_novel_junctions_machete.fasta")
        machete_style_file = open(machete_style_name, "w")
        for denovo_junction in denovo_junctions:
            jct_ind = denovo_junctions.index(denovo_junction)
            machete_style_file.write(denovo_junction.fasta_MACHETE())
        machete_style_file.close()
        write_time("-Time to write pre-collapse junctions ",start_write_pre_collapsed,timer_file_path)
 
        # Collapse the junctions that have the same splice site
        start_collapse_junctions = time.time()
        group_file_name = os.path.join(output_dir,"collapsing_group_log.txt")
        singular_jcts,collapsed_jcts = collapse_junctions(denovo_junctions,R_file_path,constants_dict,group_file_name)
        sys.stdout.write("Num singular: ["+str(len(singular_jcts))+"], num collapsed: ["+str(len(collapsed_jcts))+"]\n")
        denovo_junctions = singular_jcts+collapsed_jcts
        write_time("-Time to collapse junctions ",start_collapse_junctions,timer_file_path)


        # Get GTF information for the identified denovo_junctions
        # NOTE currently trying forward and rev-comp junction to see which one is closer to gtfs
        start_get_jct_gtf_info = time.time()
        forward_jcts = []
        reverse_jcts = []
        for jct in denovo_junctions:
            forward_jct,reverse_jct = jct.yield_forward_and_reverse()
            forward_jcts.append(forward_jct)
            reverse_jcts.append(reverse_jct)
        
        start_time = time.time()
        get_jct_gtf_info(forward_jcts,gtfs,constants_dict)
        write_time("Time to get jct gtf info 1 "+R_file_name,start_time,timer_file_path)
        
        start_time = time.time()
        get_jct_gtf_info(reverse_jcts,gtfs,constants_dict)
        write_time("Time to get jct gtf info 2 "+R_file_name,start_time,timer_file_path)
        
        gtf_denovo_junctions = []
        for jct_ind in range(len(denovo_junctions)):
            forward_jct = forward_jcts[jct_ind]
            reverse_jct = reverse_jcts[jct_ind]

            forward_dist = forward_jct.boundary_dist("donor")+forward_jct.boundary_dist("acceptor")
            reverse_dist = reverse_jct.boundary_dist("donor")+reverse_jct.boundary_dist("acceptor")
            
            #This if statement is just checking a known case where forward and reverse matter
            if (forward_jct.donor_sam.str_gene() == "NUP214" or 
                    forward_jct.acceptor_sam.str_gene() == "NUP214" or
                    reverse_jct.donor_sam.str_gene() == "NUP214" or
                    reverse_jct.acceptor_sam.str_gene() == "NUP214"):
                sys.stdout.write("Found a NUP214\n")
                sys.stdout.write(forward_jct.verbose_fasta_string())
                sys.stdout.write(str(forward_jct.donor_sam.gtf)+"\n")
                sys.stdout.write(str(forward_jct.acceptor_sam.gtf)+"\n")
                sys.stdout.write(reverse_jct.verbose_fasta_string())
                sys.stdout.write(str(reverse_jct.donor_sam.gtf)+"\n")
                sys.stdout.write(str(reverse_jct.acceptor_sam.gtf)+"\n")
                sys.stdout.write("Forward dist: "+str(forward_dist)+" reverse_dist: "+str(reverse_dist)+"\n")
                sys.stdout.write("Forward donor gtf span:"+str(forward_jct.donor_sam.gtf.span)+"\n")
                sys.stdout.write("Forward acceptor gtf span:"+str(forward_jct.acceptor_sam.gtf.span)+"\n")
                sys.stdout.write("Reverse donor gtf span:"+str(reverse_jct.donor_sam.gtf.span)+"\n")
                sys.stdout.write("Reverse acceptor gtf span:"+str(reverse_jct.acceptor_sam.gtf.span)+"\n")
                sys.stdout.write("\n")
                sys.stdout.flush()

            if forward_dist < reverse_dist:
                gtf_denovo_junctions.append(forward_jct)
            else:
                gtf_denovo_junctions.append(reverse_jct)

        write_time("Time to get full jct gtf info "+R_file_name,start_get_jct_gtf_info,timer_file_path)
        sys.stdout.write(str(len(denovo_junctions))+"\n")

        # Identify fusions from the junctions
        start_identify_fusions = time.time()
        fusion_junctions = identify_fusions(denovo_junctions,constants_dict)
        sys.stderr.write("Len fusion jcts = "+str(len(fusion_junctions))+"\n")
        write_time("Time to identify fusions "+R_file_name,start_identify_fusions,timer_file_path)

        #########################################################
        # Write out the denovo_junction_sequences for each file #
        #########################################################
        start_save_denovo_junctions = time.time()
        machete_style_name = os.path.join(output_dir,"novel_junctions_machete.fasta")
        jct_style_file_name = os.path.join(output_dir,"novel_junctions.jct")
        fusions_file_name = os.path.join(output_dir,"novel_fusions.fasta")

        # Open the three output files
        fasta_for_bowtie_index = open(fasta_for_bowtie_index_name, "w")
        machete_style_file = open(machete_style_name, "w")
        jct_style_file = open(jct_style_file_name, "w")
        fusions_file = open(fusions_file_name, "w")

        # Loop through the denovo junctions writing them where necessary
        for denovo_junction in denovo_junctions:
            jct_ind = denovo_junctions.index(denovo_junction)
            #NOTE change back to verbose_fasta_string()
            fasta_for_bowtie_index.write(denovo_junction.fasta_string())
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
        
           
write_time("Entire time",start_entire_time,timer_file_path)



