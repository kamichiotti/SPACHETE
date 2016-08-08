# General Imports
import subprocess
import itertools
import time
import sys
import os
import re

# Specific Imports
sys.path.append('./SPORK_classes/')
from consensus_utils import *
from Junction import Junction
from BinPair import BinPair
from GTFEntry import GTFEntry
from SAMEntry import SAMEntry
from FastQEntry import FastQEntry

#######################################
#   Get Reference and GTF from Mode   #
#######################################
# Finding the correct reference index based on the mode
# Human defaults are used and currently the only supported mode
# Could imagine lots of if else statements with other supported references though
# Putting this in the utils file should allow for easy reference addition
def get_reference_and_gtf_from_mode(mode):
    """
    Goal: take in the desired mode and return the index and gtf path
    Arguments:
        the mode (currently only allows hg19, so this function is mostly for show)
        will become an 'elif' tree as more and more indicies are added

    Returns:
        a tuple of (reference_path,gtf_path)
        note that regardless of the mode it is currently
        returning hg19 information
    """
    index_path = "/scratch/PI/horence/rob/index/"
    reference = index_path+"hg19"
    gtf = "/scratch/PI/horence/rob/index/hg19_genes.gtf"

    return reference,gtf


#################################
#   Find Bin Pair Group Ranges  #
#################################
# Run through the bin pair groupings to find the end index of each group
# Example bin pair list:
# [0] chr10:45_chr10:55
# [2] chr10:45_chr10:55
# [3] chr10:98_chr10:99
# [4] chr10:145_chr10:155
# [5] chr10:145_chr10:155
# [6] chr10:145_chr10:155
# 
# Output: [2,3,6]
def find_bin_pair_group_ranges(bin_pairs,constants_dict):
    """
    Goal: take in a sorted list of bin_pairs and collapse them to ranges
    Arguments:
        bin_pairs is a list[BinPair] of the bin pairs
        the constants_dict is a dictionary of global constants

    Returns:
        bin_pair_group_ranges which is a list[[start0,stop0],[start1,stop1],...]
        also filters out groups that have less then a cutoff of members
    """
    bin_pair_group_ends = []
    prev_bin_pair = ""
    for bin_pair_ind in range(len(bin_pairs)):
        curr_bin_pair = bin_pairs[bin_pair_ind].bin_pair
        if curr_bin_pair != prev_bin_pair:
            bin_pair_group_ends.append(bin_pair_ind)
            prev_bin_pair = curr_bin_pair

    # Remove the 0 at the front of the bin_pair list
    bin_pair_group_ends = bin_pair_group_ends[1:]

    # Add on the last bin_pair end which will necessarily be the end of the bin_pair list
    bin_pair_group_ends.append(len(bin_pairs))

    # Find the bin pair ranges
    bin_pair_group_ranges = []
    for bin_pair_ind in range(len(bin_pair_group_ends)-1):
        start_ind = bin_pair_group_ends[bin_pair_ind]
        stop_ind = bin_pair_group_ends[bin_pair_ind+1]
        if stop_ind-start_ind >= (constants_dict["group_member_cutoff"]-1): #The "-1" is for inclusive counting
            bin_pair_group_ranges.append([start_ind,stop_ind])
    #print bin_pair_group_ends
    #print bin_pair_group_ranges
    #sys.stdout.flush()
    return bin_pair_group_ranges


################################
#   Build Junction Sequences   #
################################
# Run through the bin_pair_group_ends to perform:
# (1) checking if there are enough in each group
# (2) padding the sequences in the first bin
# (3) creating a consensus in the first bin
# (4) scoring the consensus in the first bin
# (5) repeat (2)-(4) for the second bin
# (6) average the two consensus scores and see if they are below a cutoff
def build_junction_sequences(bin_pairs,bin_pair_group_ranges,full_path_name,constants_dict):
    """
    Goal: convert the bin_pairs into a junction object for each bin_pair
    Arguments:
        bin_pairs is a list[BinPairs]
        bin_pair_group_ranges is a list[[start0,stop0],[start1,stop1],...] so its a list of list
        full_path_name is the 'spaces-removed' file that is looked at to get the full sequences
            -> this is used to build the consensus for the junction
        the constants_dict is a dictionary of global constants

    Returns:
        denovo_junctions which is a list[Junction]
    """
    group_member_cutoff = constants_dict["group_member_cutoff"]
    consensus_score_cutoff = constants_dict["consensus_score_cutoff"]
    bin_size = constants_dict["bin_size"]
    reference = constants_dict["reference"]
    denovo_junctions = []

    # Look back at the original full path to get seq lines
    unaligned_file = open(full_path_name,"r")
    unaligned_reads = unaligned_file.readlines() 
    unaligned_file.close()
    unaligned_reads = [unaligned_reads[ind] for ind in range(len(unaligned_reads)) if ind%4 == 0 or ind%4 == 1]
    id_to_seq = {}

    # Build the dictionary of read_id to the full read
    for ind in range(0,len(unaligned_reads),2):
        key = unaligned_reads[ind].replace("\n","")
        value = unaligned_reads[ind+1].replace("\n","")
        id_to_seq[key] = value

    # walk through each bin_pair_group
    write_time("Working on the jcts :"+str(len(bin_pair_group_ranges)),time.time(),constants_dict["timer_file_path"])
    for bin_pair_group_range in bin_pair_group_ranges:
        #junction_num = "("+str(bin_pair_group_ranges.index(bin_pair_group_range)+1)+"/"+str(len(bin_pair_group_ranges))+")"
        #print junction_num
        #sys.stdout.flush()

        #start_build_junction = time.time()
        start_ind = bin_pair_group_range[0]
        stop_ind = bin_pair_group_range[1]
        group_members = bin_pairs[start_ind:stop_ind]

        # If there are not enough group members then skip this group
        if len(group_members) < group_member_cutoff:
            #sys.stderr.write("Skipped group in build junction seqs\n") #expecting to have filtered out before this
            #sys.stderr.write("len(group_members) == "+str(len(group_members))+"\n")
            continue

        # Otherwise start thinking about getting strandedness
        five_prime_strand = group_members[0].five_prime_strand
        three_prime_strand = group_members[0].three_prime_strand

        # Find the consensus sequence and score
        # Takes just the 5' ends to get the pos
        # the full original unaligned seq
        mapped_reads = [member.five_prime_SAM for member in group_members]
        #start_build_consensus = time.time()
        bin_consensus,bin_score = build_and_score_consensus(mapped_reads,five_prime_strand,id_to_seq,bin_size,constants_dict)
        #write_time("--Time to build a single consensus n="+str(len(mapped_reads))+": "+junction_num,start_build_consensus,constants_dict["timer_file_path"])

        # TODO what do I do if the five and three prime strands are not the same?
        # If the reverse compliment was taken above then take the rev compliment of the consensus too
        took_reverse_compliment = False
        if five_prime_strand == "-" and three_prime_strand == "-":
            group_members = [member.take_reverse_compliment() for member in group_members]
            bin_consensus = reverse_compliment(bin_consensus)
            took_reverse_compliment = True

        # If the bin score is high enough then add it
        if bin_score < consensus_score_cutoff:
            #__init__(consensus,score,bin_pair_group,took_reverse_compliment,constants_dict):
            denovo_junction = Junction(bin_consensus,bin_score,group_members,took_reverse_compliment,constants_dict)
            denovo_junctions.append(denovo_junction)


        #write_time("-Time to completely process a single junction: "+junction_num,start_build_junction,constants_dict["timer_file_path"])
    return denovo_junctions


########################
#   Find Splice Inds   #
########################
# Runs bowtie on all of the possible splice sites of all possible junctions
# Returns a dict keyed by jct_id and valued by a list of cut sites
def find_splice_inds(denovo_junctions,constants_dict):
    """
    Goal: find where in the consensus sequence to make and upstream and downstream cut
    Arguments:
        denovo_junctions is a list[Junction]
        the constants_dict is a dictionary of global constants

    Returns:
        returns a tuple of (jcts_with_splice,jcts_without_splice)
        to allow for continuing with only jcts that had a splice site found
    """
    # Gather info from the constants dictionary
    splice_finder_temp_name = constants_dict["out_dir"]+"splice_finder_temp_"
    min_score = constants_dict["splice_finding_min_score"]
    max_mismatches = int(constants_dict["splice_finding_allowed_mismatches"])
    read_gap_score = constants_dict["read_gap_score"]
    ref_gap_score = constants_dict["ref_gap_score"]
    allowed_mappings = constants_dict["splice_finding_allowed_mappings"]
    num_threads = constants_dict["num_threads"]
    reference = constants_dict["reference"]
    use_prior = constants_dict["use_prior"]

    # Handle temporary files
    five_prime_mapped_name = splice_finder_temp_name+"5_prime.sam"
    three_prime_mapped_name = splice_finder_temp_name+"3_prime.sam"
    five_prime_fa_file = splice_finder_temp_name+"5_prime.fa"
    three_prime_fa_file = splice_finder_temp_name+"3_prime.fa"
    five_temp_file = open(five_prime_fa_file,"w")
    three_temp_file = open(three_prime_fa_file,"w")

    # Do all the aligning work only if there are no mapped files already
    if use_prior and os.path.isfile(five_prime_mapped_name) and os.path.isfile(three_prime_mapped_name):
        write_time("--Using existing files in splice ind id",time.time(),constants_dict["timer_file_path"])
    else:
        # Write out all the possible splice sites for every jct out to a 5' and 3' file
        for jct_ind in range(len(denovo_junctions)):
            sys.stdout.flush()
            junction = denovo_junctions[jct_ind]
            cons_len = len(junction.consensus)
            splice_map_size = len(junction.consensus)/3

            five_prime_list = [junction.consensus[ind:ind+splice_map_size] for ind in range(0,cons_len-2*splice_map_size+1)]
            three_prime_list = [junction.consensus[ind:ind+splice_map_size] for ind in range(splice_map_size,cons_len-splice_map_size+1)]

            five_prime_fa_list = [">jct_"+str(jct_ind)+"_ind_"+str(ind)+"\n"+five_prime_list[ind] for ind in range(len(five_prime_list))]
            three_prime_fa_list = [">jct_"+str(jct_ind)+"_ind_"+str(ind)+"\n"+three_prime_list[ind] for ind in range(len(three_prime_list))]

            five_temp_file.write("\n".join(five_prime_fa_list)+"\n")
            three_temp_file.write("\n".join(three_prime_fa_list)+"\n")

        # Don't forget to close the files :)
        five_temp_file.close()
        three_temp_file.close()
       
        # Map the temp files above to the reference and save in temp sam files
        # Need to specify the -f flag because the inputs are fasta files
        with open(five_prime_mapped_name,"w") as five_prime_mapped:
            subprocess.call(
                ["bowtie2", "-f", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
                 "-x", reference, five_prime_fa_file], stdout=five_prime_mapped)
        with open(three_prime_mapped_name,"w") as three_prime_mapped:
            subprocess.call(
                ["bowtie2", "-f", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
                 "-x", reference, three_prime_fa_file], stdout=three_prime_mapped)

        # Sort the temp output files after removing the header lines
        p1 = subprocess.Popen(["grep","-v","@",five_prime_mapped_name],stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["sort","-o",five_prime_mapped_name],stdin=p1.stdout)
        p2.communicate()

        p3 = subprocess.Popen(["grep","-v","@",three_prime_mapped_name],stdout=subprocess.PIPE)
        p4 = subprocess.Popen(["sort","-o",three_prime_mapped_name],stdin=p3.stdout)
        p4.communicate()

    #Open the mapped, sorted, headerless files
    five_prime_sam_file = open(five_prime_mapped_name,"r")
    three_prime_sam_file = open(three_prime_mapped_name,"r")

    # Try and find the best splice sites walking through both files at once
    best_splices = {}
    sam_five_list = []
    sam_three_list = []
    prev_jct_ind = 0
    five_sam_line = five_prime_sam_file.readline()
    three_sam_line = three_prime_sam_file.readline()
    while five_sam_line and three_sam_line:
        sys.stdout.flush()
        five_sam_entry = SAMEntry(five_sam_line)
        three_sam_entry = SAMEntry(three_sam_line)
        five_jct_ind = int(five_sam_entry.read_id.split("_")[1])
        three_jct_ind = int(three_sam_entry.read_id.split("_")[1])

        # Both the 5' and 3' sam are the prev_jct_ind
        if five_jct_ind == prev_jct_ind and three_jct_ind == prev_jct_ind:
            if five_sam_entry.num_mismatches <= max_mismatches:
                sam_five_list.append(five_sam_entry)
            if three_sam_entry.num_mismatches <= max_mismatches:
                sam_three_list.append(three_sam_entry)
            five_sam_line = five_prime_sam_file.readline()
            three_sam_line = three_prime_sam_file.readline()

        # Only the 5' sam is at the prev_jct_ind
        elif five_jct_ind == prev_jct_ind:
            if five_sam_entry.num_mismatches <= max_mismatches:
                sam_five_list.append(five_sam_entry)
            five_sam_line = five_prime_sam_file.readline()

        # Only the 3' sam is at the prev_jct_ind
        elif three_jct_ind == prev_jct_ind:
            if three_sam_entry.num_mismatches <= max_mismatches:
                sam_three_list.append(three_sam_entry)
            three_sam_line = three_prime_sam_file.readline()

        # Niether the 5' nor 3' sam is at the prev_jct_ind
        else:
            prev_consensus = denovo_junctions[prev_jct_ind].consensus
            best_splices[prev_jct_ind] = get_best_splice(sam_five_list,sam_three_list,prev_consensus,max_mismatches)
            sam_five_list = []
            sam_three_list = []
            prev_jct_ind = min(five_jct_ind,three_jct_ind)

    # Have to push the last jct lists into the shared_dict
    prev_consensus = denovo_junctions[prev_jct_ind].consensus
    best_splices[prev_jct_ind] = get_best_splice(sam_five_list,sam_three_list,prev_consensus,max_mismatches)
    
    # Close the 5' and 3' sam files
    five_prime_sam_file.close()
    three_prime_sam_file.close()

    # Loop through the jcts assigning the splice site info
    jcts_with_splice = []
    jcts_without_splice = []
    for jct_ind in range(len(denovo_junctions)):
        sys.stdout.flush()
        jct = denovo_junctions[jct_ind]

        # Use calculted best splice if was previously found
        if jct_ind in best_splices:
            upstream_sam,downstream_sam = best_splices[jct_ind]
            jct.upstream_sam = upstream_sam
            jct.downstream_sam = downstream_sam
            jcts_with_splice.append(jct)
        else:
            jcts_without_splice.append(jct)

    #Return the jcts w/ and w/out splice separately
    return jcts_with_splice,jcts_without_splice

#######################
#   Get Best Splice   #
#######################
# This is a helper function for the splice site finder
# If there are multiple found splice, return the best from the shared sams dict from the two sam lists
def get_best_splice(sam_five_list,sam_three_list,consensus,max_mismatches):
    """
    Goal: given a list of possible splice, return the best one for a junction
    Arguments:
        sam_five_list is a list[SAMEntry] for possible 5' cuts
        sam_three_list is a list[SAMEntry] for possible 3' cuts
        consensus is a str and the consensus sequence of this junction
        max_mismatches is the maximum allowed mismatches in total for the 5' and 3' sides

    Returns:
        a tuple of the (best_5'_sam,best_3'_sam) to then be stored in the junction
    """
    shared_dict = {}
    id_dict = {}
    best_sam_five = SAMEntry()
    best_sam_five_len = 0
    best_sam_three = SAMEntry()
    best_sam_three_len = 0

    # Populate the id_dict to find perfectly matched 5' and 3' splices
    # Also keep track of the best (longest) 5' mapping
    for sam_five in sam_five_list:
        id_dict[sam_five.read_id] = sam_five
        if len(sam_five.seq) > best_sam_five_len:
            best_sam_five_len = len(sam_five.seq)
            best_sam_five = sam_five

    # Check the id_dict to see if the 3' has a perfect 5' match
    # Also keep track of the best (longest) 3' mapping
    for sam_three in sam_three_list:
        if len(sam_three.seq) > best_sam_three_len:
            best_sam_three_len = len(sam_three.seq)
            best_sam_three = sam_three
        if sam_three.read_id in id_dict:
            sam_five = id_dict[sam_three.read_id]
            shared_dict[sam_three.read_id] = [sam_five,sam_three]

    # Now pick out the best sam's to use
    # If there is at least one shared perfect splice ind find the one w/ least mismatches
    # then choose the perfect pair that has the least mismatches
    if len(shared_dict) > 0:
        best_key = ""
        min_mismatches = max_mismatches+1
        for key in shared_dict:
            sam1,sam2 = shared_dict[key]
            num_mismatches = sam1.num_mismatches+sam2.num_mismatches
            if best_key == "" or num_mismatches < min_mismatches:
                best_key = key
                min_mismatches = num_mismatches
        best_five_sam,best_three_sam = shared_dict[best_key]
        return best_five_sam,best_three_sam

    # Otherwise there is a mapping for the left and right pieces, although there is space in between
    # NOTE is it possible that the best 5' and 3' seqs have overlap in the middle?
    else:
        return best_sam_five,best_sam_three


#####################
#   Generate GTFS   #
#####################
# Helper function to generate a list of gtf objects from a gtf_file
def generate_gtfs(gtf_file_name):
    """
    Goal: simply open the gtf_file and put gtf objects in a list
    Arguments:
        gtf_file_name is the full path to the gtf file

    Returns:
        gtfs is a list[GTFEntry]
    """
    gtfs = []
    with open(gtf_file_name,"r") as gtf_file:
        for gtf_line in gtf_file.readlines():
            gtf = GTFEntry(gtf_line)
            gtfs.append(gtf)

    return gtfs


##############################################
#Quick check that the gtfs are in fact sorted#
##############################################
#A sanity check function
def check_gtfs_sorted(gtf_file_name):
    """
    Goal: this function is not necessary, but checks to ensure a gtf_file is sorted
    Arguments:
        gtf_file_name is the full path to the gtf file

    Returns:
        gtf_sorted which is a boolean telling is the file was sorted or not
    """
    with open(gtf_file_name,"r") as gtf_file:
        #NOTE assumes at least one line in gtf_file (ok assumption)
        gtf_line = gtf_file.readline()
        gtf_line = gtf_line.strip()
        prev_gtf = GTFEntry(gtf_line)

        gtf_sorted = True
        gtf_line = gtf_file.readline()
        while gtf_line:
            gtf_line = gtf_line.strip()
            gtf = GTFEntry(gtf_line)
            if gtf.stop < prev_gtf.stop:
                gtf_sorted = False
                sys.stderr.write("ERROR: gtfs not sorted:\n")
                sys.stderr.write("\tprev-gtf: "+str(prev_gtf)+":\n")
                sys.stderr.write("\tcurr-gtf: "+str(gtf)+":\n\n")
            else:
                sys.stdout.write("correct: gtfs sorted:\n")
                sys.stdout.write("\tprev-gtf: "+str(prev_gtf)+":\n")
                sys.stdout.write("\tcurr-gtf: "+str(gtf)+":\n\n")
            prev_gtf = gtf
            gtf_line = gtf_file.readline()

    return gtf_sorted

#################################
#   Search Closest GTF wrapper  #
#################################
#NOTE this function is repetitive
#I could improve it by adding another function later
def find_closest_gtf(jct,chrom_gtfs_start,chrom_gtfs_stop):
    """
    Goal: make it easier to call find closest gtf of upstream and downstream w/out code duplication
    Arguments:
        takes in a single junction
        takes in chrom_gtfs_start which is a dict["chrom":list[GTFEntry]] sorted by start
        takes in chrom_gtfs_start which is a dict["chrom":list[GTFEntry]] sorted by stop

    Returns:
        closest_results which is a dictionary keyed by
        -> closest_results["upstream"] -> GTFEntry
        -> closest_results["downstream"] -> GTFEntry
    """
    closest_results = {"upstream":None,"downstream":None}

    #Find the closest upstream sam gtf checking start and stop of the gtfs
    if jct.upstream_sam.exists:
        query = jct.upstream_sam.stop
        chrom = jct.upstream_sam.chromosome
        if chrom in chrom_gtfs_start and chrom in chrom_gtfs_stop:
            closest_results["upstream"] = bin_find_closest_gtf_helper(query,chrom,chrom_gtfs_start,chrom_gtfs_stop)

    #Find the closest downstream sam gtf checking start and stop of the gtfs
    if jct.downstream_sam.exists:
        query = jct.downstream_sam.start
        chrom = jct.downstream_sam.chromosome
        if chrom in chrom_gtfs_start and chrom in chrom_gtfs_stop:
            closest_results["downstream"] = bin_find_closest_gtf_helper(query,chrom,chrom_gtfs_start,chrom_gtfs_stop)

    return closest_results

#############################################
#  Brute Force Search Closest GTF Helper    #
#############################################
#Checks every gtf on the correct chromosome
def brute_find_closest_gtf_helper(query,chrom,chrom_gtfs_start,chrom_gtfs_stop):
    """
    Goal: find the closest gtf to a given coordinate
    Arguments:
        query is an integer of the genomic coordinate
        chrom is a string of the chromosome
        takes in chrom_gtfs_start which is a dict["chrom":list[GTFEntry]] sorted by start
        takes in chrom_gtfs_start which is a dict["chrom":list[GTFEntry]] sorted by stop

    Returns:
        the single closest GTFEntry object
    """
    closest_start_gtf = None
    closest_stop_gtf = None
    closest_start_dist = -1
    closest_stop_dist = -1
    for gtf in chrom_gtfs_start[chrom]:
        if abs(gtf.start-query) < closest_start_dist or not closest_start_gtf:
            closest_start_dist = abs(gtf.start-query)
            closest_start_gtf = gtf
        if abs(gtf.stop-query) < closest_stop_dist or not closest_stop_gtf:
            closest_stop_dist = abs(gtf.stop-query)
            closest_stop_gtf = gtf

    #Decide which one is closer to return
    if closest_start_dist < closest_stop_dist:
        return closest_start_gtf
    else:
        return closest_stop_gtf


########################################
#  Binary Search Closest GTF Helper    #
########################################
#A wrapper function which calls the recursive gtf binary search function
def bin_find_closest_gtf_helper(query,chrom,chrom_gtfs_start,chrom_gtfs_stop):
    """
    Goal: prepare input and handle output from the recursive binary search function
    Arguments:
        query is an integer genomic location
        chrom is a string chromosome selection
        takes in chrom_gtfs_start which is a dict["chrom":list[GTFEntry]] sorted by start
        takes in chrom_gtfs_start which is a dict["chrom":list[GTFEntry]] sorted by stop

    Returns:
        the single closest GTFEntry object
    """
    start_lib = [gtf.start for gtf in chrom_gtfs_start[chrom]]
    stop_lib = [gtf.stop for gtf in chrom_gtfs_stop[chrom]]

    closest_start,its_start = bin_search_gtf(query,start_lib)
    closest_stop,its_stop = bin_search_gtf(query,stop_lib)

    start_gtf = chrom_gtfs_start[chrom][closest_start] 
    stop_gtf = chrom_gtfs_stop[chrom][closest_stop]
   
    start_dist = abs(query-closest_start)
    stop_dist = abs(query-closest_stop)

    ret_val = start_gtf if start_dist < stop_dist else stop_gtf
    return ret_val
   

#################################
#    Binary Search Closest GTF  #
#################################
#Recursive binary search function
def bin_search_gtf(query,library,start_ind=0,end_ind=-1,its=1,disp=False):
    """
    Goal: do a binary search through the library for the closest ind
    Arguments:
        query is an int
        library is a list[int] of genomic positions to look through
        start_ind is defaulted to 0 and keeps track of where to look
        stop_ind is defaulted to -1 and keeps track of where to look

    Returns:
        the index of the closest matching library value to the query
    """
    #First iteration
    if its == 1 and end_ind == -1:
        end_ind = len(library)-1

    #Base case
    if end_ind-start_ind <= 1:
        sys.stdout.flush()
        start_dist = abs(query-library[start_ind])
        end_dist = abs(library[end_ind]-query)
        ret_ind = start_ind if start_dist <= end_dist else end_ind
        if disp:
            print str(its)+")","Found closest to:[",library[ret_ind],"]"
        return ret_ind,its
    #Recursive case
    else:
        mid_ind = (start_ind+end_ind)/2
        #Determine whether or not to print this line
        if disp:
            print str(its)+")",library[start_ind],"--",library[mid_ind],"--",library[end_ind]
            
        #Check to see how to recurse
        if query < library[mid_ind]:
            return bin_search_gtf(query,library,start_ind,mid_ind,its+1,disp=disp)
        else:
            return bin_search_gtf(query,library,mid_ind,end_ind,its+1,disp=disp) 


########################
#   Get JCT GTF info   #
########################
# Reads the entire gtf line by line and checks against every junction
# Reading the genes in as all exons to try and get start and stop site
def get_jct_gtf_info(junctions,gtfs):
    """
    Goal: for each junction find the closest gtf for upstream and downstream
    Arguments:
        junctions is a list[Junction]
        gtfs is a list[GTF]

    Returns:
        nothing (just updates the junction objects as pass-by-reference)
    """
    # Separate the gtfs by chromosome into a dictionary
    chrom_gtfs = {}
    for gtf in gtfs:
        if gtf.chromosome not in chrom_gtfs:
            chrom_gtfs[gtf.chromosome] = [gtf]
        else:
            chrom_gtfs[gtf.chromosome].append(gtf)

    # Pre sort the gtfs into a start and stop oriented list by chromosome
    chrom_gtfs_start = {}
    chrom_gtfs_stop = {}
    for chrom in chrom_gtfs:
        chrom_gtfs_start[chrom] = sorted(chrom_gtfs[chrom],key=lambda gtf: gtf.start)
        chrom_gtfs_stop[chrom] = sorted(chrom_gtfs[chrom],key=lambda gtf: gtf.stop)

    # Loop through the collapsed gtfs to see if the junction is in the range
    for junction in junctions:
        jct_ind = junctions.index(junction)
        closest_results = find_closest_gtf(junction,chrom_gtfs_start,chrom_gtfs_stop)
        if closest_results["upstream"]:
            gtf = closest_results["upstream"]
            junction.upstream_sam.gtf = gtf
            
        if closest_results["downstream"]:
            gtf = closest_results["downstream"]
            junction.downstream_sam.gtf = gtf
            

#####################################
#         Identify Fusions          #
#####################################
#Takes junctions that already have gtf info
#If a junction has the following properties call it a 'fusion':
#   If upstream and downstream sams are at_boundary:
#       If upstream and downstream are on different chromosomes
#           [yes fusion]
#       ELIF upstream and downstream are on different strands
#       ELIF distance between upstream and downstream > threshold
#           [yes fusion]
#       ELSE
#           [no fusion]
#   Else:
#       [no fusion]
#
#Returns a list of junctions that are deemed 'fusions'
def identify_fusions(junctions,span_cutoff=1e6):
    """
    Goal: take the junctions and find the ones that could be fusions
    Arguments:
        junctions is a list[Junction] objects
        span_cutoff is an optional int or float defining min distance for a fusion
            on the same chromosome
    Returns:
        fusion_jcts as a list[Junction] with the junctions defined as fusions
    """
    fusion_jcts = []
    for jct in junctions:
        if jct.at_boundary("upstream") and jct.at_boundary("downstream"):
            if jct.upstream_sam.chromosome != jct.downstream_sam.chromosome:
                fusion_jcts.append(jct)
            elif jct.upstream_sam.strand != jct.downstream_sam.strand:
                fusion_jcts.append(jct)
            elif abs(jct.span()) > span_cutoff:
                fusion_jcts.append(jct)
    
    return fusion_jcts

#####################################
# Collapse junctions by splice site #
#####################################
# Currently keeps linear and non-linear separated even if they share a splice site
# NOTE currently does not work
def old_collapse_junctions(junctions):
    """
    Goal: take the junctions and collapse ones which represent the same splice site
    Arguments:
        junctions is a list[Junction] objects

    Returns:
        collapsed_junctions is a list[Junction] objects
        the returned list will always have equal to or fewer objects
        than the input list
    """
    splice_to_jct_dict = {}
    for junction in junctions:
        splice_key = junction.upstream_chromosome+"-"+junction.downstream_chromosome+":"+str(junction.splice_site)+":"+str(junction.linear)
        if splice_key not in splice_to_jct_dict:
            splice_to_jct_dict[splice_key] = [junction]
        else:
            splice_to_jct_dict[splice_key].append(junction)
            
    collapsed_junctions = []
    for jcts_by_splice in splice_to_jct_dict.itervalues():
        first_junction = jcts_by_splice[0]
        for shared_jct in jcts_by_splice[1:]:
            first_junction.add_constitutive_junction(shared_jct)
        first_junction.combine_constitutive_junctions()
        collapsed_junctions.append(first_junction)
    return collapsed_junctions


####################
#   Assign Class   #
####################
# Assigns a pair of R1 and R2 to a class based on certain factors. These are both SAMEntry objects.
# Possible classes. R1 always maps to a denovo jct, and R2 somewhere else.
# NOTE I currently have no confidence that this function works
# TODO update the logic to allow 'Fusion' classification
# TODO use regex to parse the chromosome rather than all this messy splitting (it looks terrible)
# [1] Linear
# [2] Linear Anomally
# [3] Circular
# [4] Circular Anomally
# [5] None <-- kind of in the place of fusions for now
def assign_class(sam_R1,sam_R2):
    """
    Goal: take the upstream and downstream sam and categorize them
    Arguments:
        sam_R1 is of type SAMEntry
        sam_R2 is of type SAMEntry

    Returns:
        a string of the generated type
    """
    jct_chrom_1 = sam_R1.chromosome.split("|_|")[0].split("|")[1]
    jct_chrom_2 = sam_R1.chromosome.split("|_|")[1].split("|")[0]

    # If the jct splices 2 chromosomes together just skip it for now
    if jct_chrom_1 != jct_chrom_2:
        return "None"
    # If the R1 and R2 are on different chromosomes just skip it for now
    if jct_chrom_1 != sam_R2.chromosome:
        return "None"

    span = int(sam_R1.chromosome.split("|_|")[1].split("|")[5].split(":")[1])

    #Linear case
    if span > 0:
        if sam_R1.strand != sam_R2.strand:
            return "Linear"
        else:
            return "Linear_Anomaly"

    #Circular case
    else:
        start_1 = int(sam_R1.chromosome.split("|_|")[0].split("|")[3].split(":")[1])
        stop_1 = int(sam_R1.chromosome.split("|_|")[0].split("|")[4].split(":")[1])
        start_2 = int(sam_R1.chromosome.split("|_|")[1].split("|")[2].split(":")[1])
        stop_2 = int(sam_R1.chromosome.split("|_|")[1].split("|")[3].split(":")[1])
        if sam_R1.strand == sam_R2.strand:
            return "Circular_Anomaly"
        elif start_2 <= sam_R2.start <= stop_1:
            return "Circular"
        else:
            return "Circular_Anomaly"


############################
#   Write GLM Class File   #
############################
# Simple function to print out GLM class file in the right format
# NOTE I don't think this function currently works
def write_glm_class_file(class_file_name,sam_list):
    """
    Goal: print out a class file for the GLM
    Arguments:
        class_file_name is the name of the save file for the GLM
        sam_list is a list[[SAMEntry,SAMEntry,str,list[string],list[string]],...]

    Returns:
        nothing (just prints to the file instead)
    """
    header = ""
    header += "id\t"
    header += "class\t"
    header += "pos\t"
    header += "qual\t"
    header += "aScore\t"
    header += "numN\t"
    header += "readLen\t"
    header += "junction\t"
    header += "strand\t"
    header += "posR2\t"
    header += "qualR2\t"
    header += "aScoreR2\t"
    header += "numNR2\t"
    header += "readLenR2\t"
    header += "junctionR2\t"
    header += "strandR2\n"
    with open(class_file_name,"w") as class_file:
        class_file.write(header)
        for read_pair in sam_list:
            sam_R1,sam_R2,pair_class,sam_R1_genes,sam_R2_genes = read_pair

            #Add general info to the out_line
            out_line = ""
            out_line += str(sam_R1.read_id.split("/")[0])+"\t"
            out_line += str(pair_class)+"\t"

            #Add R1 info to the out_line
            out_line += str(sam_R1.start)+"\t"
            out_line += str(sam_R1.mapping_quality)+"\t"
            out_line += str(sam_R1.alignment_score)+"\t"
            out_line += str(sam_R1.num_Ns)+"\t"
            out_line += str(len(sam_R1.seq))+"\t"
            out_line += str(sam_R1.junction())+"|"+sam_R1_genes+"\t"
            out_line += str(sam_R1.strand)+"\t"

            #Add R2 info to the out_line
            out_line += str(sam_R2.start)+"\t"
            out_line += str(sam_R2.mapping_quality)+"\t"
            out_line += str(sam_R2.alignment_score)+"\t"
            out_line += str(sam_R2.num_Ns)+"\t"
            out_line += str(len(sam_R2.seq))+"\t"
            out_line += str(sam_R2.junction())+"|"+sam_R2_genes+"\t"
            out_line += str(sam_R2.strand)+"\n"

            #Write the built up out_line to the glm input class file
            class_file.write(out_line)


##################
#   Write Time   #
##################
# Helper function to write out the timing of something
# Takes in the message, a start time in seconds, and a timer_file_path
# Appends to the timer file by default, first call should overwrite
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

##########################
#   Collapse Junctions   #
##########################
def collapse_junctions(jcts,full_path_name,constants_dict):
    """
    Goal: take in the junctions and collapse ones at or near the same
          splice sites
    Arguments:
        junctions is a list[Junction] object
        full_path_name points to the the combined fastq file
        constants dict is a dictionary of global constants

    Returns:
        a tuple with the first element as the singles "uncollapsed list"
        and the second being the collapsed list of junctions list[Junction]
    """
    #Get the collapsing threshold (radius of donor and acceptor to collapse in)
    collapse_thresh = constants_dict["collapse_thresh"]

    #Separate the jcts by chromosome pairs to make later O(N^2) less painful
    #so will have one entry for chr1:chr2, chr1:chr3 etc
    #it will be combinations, not permutations (chr1:chr2 == chr2:chr1)
    splices_by_chroms = {}
    for jct in jcts:
        chrom_1 = str(jct.upstream_sam.chromosome)
        chrom_2 = str(jct.downstream_sam.chromosome)
        if chrom_1+chrom_2 in splices_by_chroms:
            splices_by_chroms[chrom_1+chrom_2].append(jct)
        elif chrom_2+chrom_1 in splices_by_chroms:
            splices_by_chroms[chrom_2+chrom_1].append(jct)
        else:
            splices_by_chroms[chrom_1+chrom_2] = [jct]

    groupings = {}
    for chroms in splices_by_chroms:
        groupings[chroms] = []
        for jct in splices_by_chroms[chroms]:
            don = jct.upstream_sam.stop
            acc = jct.downstream_sam.start
            found_prev_group = False

            #Only compare jct to other jcts if both
            #don and acc are not None
            if don and acc:
                for prev_group in groupings[chroms]:
                    for prev_jct in prev_group:
                        prev_don = prev_jct.upstream_sam.stop
                        prev_acc = prev_jct.downstream_sam.start
                        #If any one of the acceptor/donors are None
                        if not prev_don or not prev_acc:
                            continue
                        if abs(don-prev_don) <= collapse_thresh and abs(acc-prev_acc) <= collapse_thresh:
                            sys.stderr.write("Found match for:\n")
                            sys.stderr.write(jct.verbose_fasta_string())
                            prev_group.append(jct)
                            found_prev_group = True
                            break

                    #If found a prev_group, don't need to look through
                    #other prev groups
                    if found_prev_group:
                        break

            #If it didn't find any of the prev_groups, start a new group
            if not found_prev_group:
                groupings[chroms].append([jct])


    #Separate singles from groups
    singles = []
    groups = []
    sys.stdout.write("STARTING REP PRINT\n")
    for chroms in groupings:
        for group in groupings[chroms]:
            if len(group) <= 1:
                singles += group
            else:
                sys.stderr.write("Group info:\n")
                sys.stderr.write("".join([m.verbose_fasta_string() for m in group]))
                sys.stderr.write("\n")
                counts = [len(member.bin_pair_group) for member in group]
                max_ind = counts.index(max(counts))
                sys.stdout.write(group[max_ind].verbose_fasta_string()+"\n")
                groups.append(group[max_ind])

    return singles,groups
    

##########################
#   Reverse Compliment   #
##########################
# Just a little helper function to give the reverse compliment of a sequence
def reverse_compliment(seq):
    """
    Goal: take a sequence and return the reverse compliment
    Arguments:
        seq is a string of A's,T's,C's,G's, and N's

    Returns:
        the reverse compliment string
    """
    comp_dict = {"A":"T",
                 "a":"t",
                 "T":"A",
                 "t":"a",
                 "C":"G",
                 "c":"g",
                 "G":"C",
                 "g":"c",
                 "N":"N",
                 "n":"n"}
    rev_comp_seq = "".join([comp_dict[base] for base in seq])[::-1]
    return rev_comp_seq



