# General Imports
import subprocess
import itertools
import zlib
import time
import sys
import os
import re

# Specific Imports
from SPORK_consensus_utils import *
from SPORK_Junction import Junction
from SPORK_BinPair import BinPair
from SPORK_GTFEntry import GTFEntry
from SPORK_SAMEntry import SAMEntry
from SPORK_FastQEntry import FastQEntry

#######################################
#   Get Reference and GTF from Mode   #
#######################################
# Finding the correct reference index based on the mode
# Human defaults are used and currently the only supported mode
# Could imagine lots of if else statements with other supported references though
# Putting this in the utils file should allow for easy reference addition
def get_reference_and_gtf_from_mode(ref_dir,abs_path,mode="hg19"):
    """
    Goal: take in the desired mode and return the index and gtf path
    Arguments:
        the path to the circ_ref directory that MACHETE also uses
        the mode (currently only allows hg19, so this function is mostly for show)
        will become an 'elif' tree as more and more indicies are added

    Returns:
        a tuple of (reference_path,gtf_path)
        note that regardless of the mode it is currently
        returning hg19 information
    """
    #index_path = "/scratch/PI/horence/rob/index/"
    #reference = index_path+"hg19"
    reference = ref_dir
    gtf_path = ""
    if mode == "hg19":
        gtf_path = os.path.join(abs_path,"gtfs","hg19_gtfs")
        reference = os.path.join(reference,"hg19_genome")

    elif mode == "mm10":
        gtf_path = os.path.join(abs_path,"gtfs","mm10_gtfs")
        reference = os.path.join(reference,"mm10_genome")

    return reference,gtf_path


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
        key = unaligned_reads[ind].strip()
        value = unaligned_reads[ind+1].strip()
        id_to_seq[key] = value

    # walk through each bin_pair_group
    write_time("Working on the bin-pairs :"+str(len(bin_pair_group_ranges)),time.time(),constants_dict["timer_file_path"])
    jct_ind = 0
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
        bin_consensus,bin_score = build_and_score_consensus(mapped_reads,five_prime_strand,id_to_seq,bin_size,constants_dict)
        took_reverse_compliment = False

        """
        # TODO what do I do if the five and three prime strands are not the same? (this represents a translocation)
        # If the reverse compliment was taken above then take the rev compliment of the consensus too
        # NOTE currently just taking reverse compliment whenever the 5' strand is negative to help groupings
        #   i.e. if have 5' - and 3' + of a jct and 5' + and 3' - of the same jct, they should be collapsed, but won't be
        #   unless I implement this
        #if five_prime_strand == "-" and three_prime_strand == "-":
        if five_prime_strand == "-":
            group_members = [member.take_reverse_compliment() for member in group_members]
            bin_consensus = reverse_compliment(bin_consensus)
            took_reverse_compliment = True

        """
        # If the bin score is good enough then add it
        if bin_score < consensus_score_cutoff:
            denovo_junction = Junction(bin_consensus,bin_score,group_members,jct_ind,took_reverse_compliment,constants_dict)
            denovo_junctions.append(denovo_junction)
            jct_ind += 1

    return denovo_junctions


########################
#   Find Splice Inds   #
########################
# Runs bowtie on all of the possible splice sites of all possible junctions
# Returns a dict keyed by jct_id and valued by a list of cut sites
def find_splice_inds(denovo_junctions,constants_dict):
    """
    Goal: find where in the consensus sequence to make and donor and acceptor cut
    Arguments:
        denovo_junctions is a list[Junction]
        the constants_dict is a dictionary of global constants

    Returns:
        returns a tuple of (jcts_with_splice,jcts_without_splice)
        to allow for continuing with only jcts that had a splice site found
    """
    # Gather info from the constants dictionary
    splice_finder_temp_name = os.path.join(constants_dict["output_dir"],"splice_finder_temp_")
    thirds_len = constants_dict["thirds_len"]
    min_score = constants_dict["splice_finding_min_score"]
    max_mismatches = int(constants_dict["splice_finding_allowed_mismatches"])
    read_gap_score = constants_dict["read_gap_score"]
    ref_gap_score = constants_dict["ref_gap_score"]
    allowed_mappings = constants_dict["splice_finding_allowed_mappings"]
    num_threads = constants_dict["num_threads"]
    reference = constants_dict["reference"]
    use_prior = constants_dict["use_prior"]
    timer_file_path = constants_dict["timer_file_path"]

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
        num_size_excluded = 0
        avg_cons_len = 0
        for jct_ind in range(len(denovo_junctions)):
            sys.stdout.flush()
            junction = denovo_junctions[jct_ind]
            cons_len = len(junction.consensus)
            avg_cons_len = (float(jct_ind*avg_cons_len)/(jct_ind+1))+(float(cons_len)/(jct_ind+1))
            #splice_map_size = len(junction.consensus)/3
            #NOTE found that I was forcinng splice sites to be too central if I used the thirds len
            #NOTE and lost BCR-ABL this way, so instead I'll stick with using 1/3 of the consensus
            splice_map_size = thirds_len
            #sys.stderr.write("Cons len "+str(cons_len)+" and thirds_len "+str(splice_map_size)+"\n")
 
            five_prime_list = [junction.consensus[ind:ind+splice_map_size] for ind in range(0,cons_len-2*splice_map_size+1)]
            three_prime_list = [junction.consensus[ind:ind+splice_map_size] for ind in range(splice_map_size,cons_len-splice_map_size+1)]

            #Oh this is my problem, sometimes the consensus length is too small
            #to try and find splice inds from, so I should just throw that out and iterate,
            #otherwise I'll have blank lines added to my splice ind fa files and they will be
            #incorrectly formatted
            if len(five_prime_list) == 0 or len(three_prime_list) == 0:
                num_size_excluded += 1
                continue

            five_prime_fa_list = [">jct_"+str(jct_ind)+"_ind_"+str(ind)+"\n"
                                  +five_prime_list[ind] for ind in range(len(five_prime_list))]
            three_prime_fa_list = [">jct_"+str(jct_ind)+"_ind_"+str(ind)+"\n"
                                   +three_prime_list[ind] for ind in range(len(three_prime_list))]

            five_temp_file.write("\n".join(five_prime_fa_list)+"\n")
            three_temp_file.write("\n".join(three_prime_fa_list)+"\n")

        # Don't forget to close the files :)
        five_temp_file.close()
        three_temp_file.close()
        sys.stdout.write("SPORK: Average consensus length ["+str(avg_cons_len)+"]\n")
        sys.stdout.write("SPORK: Splitting fqs, size excluded ["+str(num_size_excluded)+"] of ["+str(len(denovo_junctions))+"]\n")
       
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
            #RB 11/18/16 Filter out 5' and 3' sams if chroms are not the same as the unsplit
            #RB 11/21/16 Also imposing a radius on the same chromosome filter (start at 30, very tight)
            radius = 30
            pj = denovo_junctions[prev_jct_ind]
            sam_five_list = [sam for sam in sam_five_list
                             if (sam.chromosome == pj.donor_sam.chromosome and 
                                 abs(pj.donor_sam.start-sam.start) < radius)]

            sam_three_list = [sam for sam in sam_three_list
                              if (sam.chromosome == pj.acceptor_sam.chromosome and 
                                  abs(pj.acceptor_sam.start-sam.start) < radius)]
            
            # Get the best 5' and 3' pair of Sams
            prev_consensus = pj.consensus
            best_five,best_three = get_best_splice(sam_five_list,sam_three_list,prev_consensus,max_mismatches)

            best_splices[prev_jct_ind] = [best_five,best_three]
            sam_five_list = []
            sam_three_list = []
            prev_jct_ind = min(five_jct_ind,three_jct_ind)

    # Have to push the last jct lists into the shared_dict
    pj = denovo_junctions[prev_jct_ind]
    prev_consensus = pj.consensus
    sam_five_list = [sam for sam in sam_five_list if sam.chromosome == pj.donor_sam.chromosome]
    sam_three_list = [sam for sam in sam_three_list if sam.chromosome == pj.acceptor_sam.chromosome]
 
    best_five,best_three = get_best_splice(sam_five_list,sam_three_list,prev_consensus,max_mismatches)
    best_splices[prev_jct_ind] = [best_five,best_three]
    
    # Close the 5' and 3' sam files
    five_prime_sam_file.close()
    three_prime_sam_file.close()

    # Loop through the jcts assigning the splice site info
    jcts_with_splice = []
    jcts_without_splice = []
    for jct_ind in range(len(denovo_junctions)):
        sys.stdout.flush()
        jct = denovo_junctions[jct_ind]

        # Use calculted best splice if was previously found and both the
        # donor and acceptor elements exist
        if jct_ind in best_splices and all([sam.exists for sam in best_splices[jct_ind]]):
            donor_sam,acceptor_sam = best_splices[jct_ind]

            #Need to include the donor seq not used in splitting
            #gets complicated by + and - strand
            donor_ind = int(donor_sam.read_id.split("_ind_")[1])
            donor_len = len(donor_sam.seq)
            up_remaining = jct.consensus[:donor_ind]
            if donor_sam.strand == "+":
                donor_sam.seq = up_remaining+donor_sam.seq
                donor_sam.start -= len(up_remaining)
            elif donor_sam.strand == "-":
                donor_sam.seq = up_remaining+reverse_compliment(donor_sam.seq)
                donor_sam.stop += len(up_remaining)
            jct.donor_sam = donor_sam

            #Need to include the acceptor seq not used in splitting
            #gets complicated by + and - strand
            acceptor_ind = int(acceptor_sam.read_id.split("_ind_")[1])
            acceptor_len = len(acceptor_sam.seq)
            down_remaining = jct.consensus[acceptor_ind+2*acceptor_len:]
            if acceptor_sam.strand == "+":
                acceptor_sam.seq = acceptor_sam.seq+down_remaining
                acceptor_sam.stop += len(down_remaining)-1
            elif acceptor_sam.strand == "-":
                acceptor_sam.seq = reverse_compliment(acceptor_sam.seq)+down_remaining
                acceptor_sam.start -= len(down_remaining)
            jct.acceptor_sam = acceptor_sam
 
            jcts_with_splice.append(jct)

        # If either the donor or acceptor doesn't map (or both), add it to the
        # jcts_without_splice list instead NOTE this list is currently not used
        else:
            jcts_without_splice.append(jct)

    #Choose either the forward or reverse form of the junction that yields
    #the smaller donor site since this will help collapsing in the next step
    #(could have chosen larger donor etc, just to flip them all same way)
    
    write_time("Num Jcts w/ splice = "+str(len(jcts_with_splice)),time.time(),timer_file_path)
    write_time("Num Jcts w/out splice = "+str(len(jcts_without_splice)),time.time(),timer_file_path)

    small_don_jcts_with_splice = []
    for jct in jcts_with_splice:
        if jct.donor_sam.donor() < jct.acceptor_sam.acceptor():
            small_don_jcts_with_splice.append(jct)
        else:
            small_don_jcts_with_splice.append(jct.yield_reverse())

    #Return the jcts w/ and w/out splice separately
    return small_don_jcts_with_splice,jcts_without_splice

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


##########################
#   Filter Map Quality   #
##########################
#Function to take in Junctions after Splice Ind finding and map pieces right around ind
#to filter out degenerate or multiple mapping pieces. Currently uses a hard threshold
def filter_map_quality(jcts, constants_dict):
    """
    Goal: Map left and right 25-mers of splice ind and filter out jcts by mapping score
    Arguments:
        jcts is a list of type Junction
        constants dict has all the constants used in the program

    Returns:
        (pass_jcts,fail_jcts, anom_jcts) is a tuple of 
            (1) list of type Junction of jcts that passed
            (2) list of type Junction of jcts that failed
            (3) list of type Junction of jcts that had some other error
    """

    #Add map qualities and only keep those above the cutoff
    mq_cutoff = constants_dict["mq_cutoff"]
    mq_len = constants_dict["mq_len"]

    #mqmallg=unique(mallg[(!is.na(match( paste(mallg$junction),paste(mq[mq3+mq5>lower.mq]$junction))))])
    
    #Write out temp fasta for Bowtie2 calls
    temp_fasta_name = os.path.join(constants_dict["output_dir"],"mapq_temp.fasta")
    with open(temp_fasta_name,"w") as temp_fasta:
        for ind,jct in enumerate(jcts):
            seq = jct.consensus
            don_seq = seq[:jct.splice_ind()][-mq_len:] #Get the last mq_len bases before splice
            acc_seq = seq[jct.splice_ind():][:mq_len]  #Get the first mq_len bases after splice
            temp_fasta.write(">jct_"+str(ind)+"_don"+"\n"+don_seq+"\n")
            temp_fasta.write(">jct_"+str(ind)+"_acc"+"\n"+acc_seq+"\n")

    #Get bowtie2 parameter constants
    min_score = constants_dict["splice_finding_min_score"]
    read_gap_score = constants_dict["read_gap_score"]
    ref_gap_score = constants_dict["ref_gap_score"]
    num_threads = constants_dict["num_threads"]
    reference = constants_dict["reference"]
    use_prior = constants_dict["use_prior"]
    timer_file_path = constants_dict["timer_file_path"]
    mq_mapped_name = os.path.join(constants_dict["output_dir"],"mapq_mapped.sam")

    #Call bowtie2 to get map qualities
    with open(mq_mapped_name,"w") as mq_mapped:
        subprocess.call(
           ["bowtie2", "-f", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
            "-x", reference, temp_fasta_name], stdout=mq_mapped)


    #Go through the mapped output file
    pass_jcts = []
    fail_jcts = []
    anom_jcts = []
    with open(mq_mapped_name,"r") as mq_mapped:
        sam_line = mq_mapped.readline()

        while sam_line:
            #Skip header lines
            if "@" in sam_line:
                sam_line = mq_mapped.readline()
                continue

            #Get sam entry 1
            sam_entry_1 = SAMEntry(sam_line)
            x,jct_ind_1,don = sam_entry_1.read_id.split("_")
            jct_ind_1 = int(jct_ind_1)

            #Get sam entry 2
            sam_line = mq_mapped.readline()
            sam_entry_2 = SAMEntry(sam_line)
            if not sam_line:
                anom_jcts.append(jcts[jct_ind_1])
                continue
            x,jct_ind_2,acc = sam_entry_2.read_id.split("_")
            jct_ind_2 = int(jct_ind_2)

            #Check for different anomalous cases where a don/acc doesn't appear in output
            #If any anomally is hit we start w/ sam_entry_2 next loop, don't advance
            if don != "don" or acc != "acc" or jct_ind_1 != jct_ind_2:
                anom_jcts.append(jcts[jct_ind_1])
                continue

            #Otherwise if they are don/acc from the same jct check if they passed mq_cutoff
            jct_mapq = sam_entry_1.mapping_quality + sam_entry_2.mapping_quality
            jcts[jct_ind_1].mapq = jct_mapq
            if jct_mapq > mq_cutoff:
                pass_jcts.append(jcts[jct_ind_1])
            else:
                fail_jcts.append(jcts[jct_ind_1])

            #Move on to the next line
            sam_line = mq_mapped.readline()

    #Return the passed, failed, and anomalous jcts
    return (pass_jcts, fail_jcts, anom_jcts)


#####################
#   Generate GTFS   #
#####################
# Helper function to generate a list of gtf objects from a gtf path full of gtf files
def generate_gtfs(gtf_path,allowed_feature_types=["exon"]):
    """
    Goal: open all the gtf_files and put all gtf objects in a list from the given path
    Arguments:
        gtf_file_name is the full path to the gtf files
        allowed_feature_types is a list of string specifiying which feature types to add (default 'exon' only)

    Returns:
        gtfs is a list[GTFEntry]
    """
    gtfs = []
    gtf_file_names = [gtf_name for gtf_name in os.listdir(gtf_path) if "gtf" in gtf_name]
    for gtf_file_name in gtf_file_names:
        abs_gtf_file_path = os.path.join(gtf_path,gtf_file_name)
        #sys.stdout.write("Reading in GTF file "+abs_gtf_file_path+"\n")
        with open(abs_gtf_file_path,"r") as gtf_file:
            for gtf_line in gtf_file.readlines():
                gtf = GTFEntry(gtf_line)
                if gtf.feature in allowed_feature_types:
                    gtfs.append(gtf)

    return gtfs

########################
#   Get JCT GTF info   #
########################
def get_jct_gtf_info(junctions,gtfs,constants_dict):
    """
    Goal: for each junction find the closest gtf for donor and acceptor
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

    # Pre sort the gtfs into a donor and acceptor oriented list by chromosome and strand
    chrom_gtfs_don = {}
    chrom_gtfs_acc = {}
    chrom_don_libs = {}
    chrom_acc_libs = {}
    for chrom in chrom_gtfs:
        chrom_gtfs_don[chrom+"+"] = []
        chrom_gtfs_don[chrom+"-"] = []
        chrom_don_libs[chrom+"+"] = []
        chrom_don_libs[chrom+"-"] = []
        for don_gtf in sorted(chrom_gtfs[chrom],key=lambda gtf: gtf.donor):
            chrom_gtfs_don[chrom+don_gtf.strand].append(don_gtf)
            chrom_don_libs[chrom+don_gtf.strand].append(don_gtf.donor)

        chrom_gtfs_acc[chrom+"+"] = []
        chrom_gtfs_acc[chrom+"-"] = []
        chrom_acc_libs[chrom+"+"] = []
        chrom_acc_libs[chrom+"-"] = []
        for acc_gtf in sorted(chrom_gtfs[chrom],key=lambda gtf: gtf.acceptor):
            chrom_gtfs_acc[chrom+acc_gtf.strand].append(acc_gtf)
            chrom_acc_libs[chrom+acc_gtf.strand].append(acc_gtf.acceptor)

    # Find the closest gtfs to donor and acceptor
    for junction in junctions:
        jct_ind = junctions.index(junction)
        #message = "Finding gtf info ("+str(jct_ind)+"/"+str(len(junctions))+")"
        #write_time(message,time.time(),constants_dict["timer_file_path"])
        closest_results = find_closest_gtf(junction,chrom_gtfs_don,chrom_gtfs_acc,chrom_don_libs,chrom_acc_libs)
        if closest_results["donor"]:
            gtf = closest_results["donor"]
            junction.donor_sam.gtf = gtf
            
        if closest_results["acceptor"]:
            gtf = closest_results["acceptor"]
            junction.acceptor_sam.gtf = gtf
            

#################################
#        Find Closest GTF       #
#################################
def find_closest_gtf(jct,chrom_gtfs_don,chrom_gtfs_acc,chrom_don_libs,chrom_acc_libs):
    """
    Goal: make it easier to call find closest gtf of donor and acceptor w/out code duplication
    Arguments:
        takes in a single junction
        takes in chrom_gtfs_don which is a dict["chrom":list[GTFEntry]] sorted by donor
        takes in chrom_gtfs_acc which is a dict["chrom":list[GTFEntry]] sorted by acceptor

    Returns:
        closest_results which is a dictionary keyed by
        -> closest_results["donor"] -> GTFEntry
        -> closest_results["acceptor"] -> GTFEntry
    """
    closest_results = {"donor":None,"acceptor":None}

    #Find the closest donor sam gtf
    if jct.donor_sam.exists:
        query = jct.donor_sam.donor()
        chrom = jct.donor_sam.chromosome
        chrom_key = chrom+jct.donor_sam.strand
        if chrom_key in chrom_gtfs_don and chrom_key in chrom_don_libs:
            gtfs_don = chrom_gtfs_don[chrom_key]
            don_lib = chrom_don_libs[chrom_key]
            #closest_don_ind,its = bin_search_gtf(query,don_lib)
            closest_don_ind = brute_search_gtf(query,don_lib) #<-- RB trying brute force
            closest_results["donor"] = gtfs_don[closest_don_ind]

    #Find the closest acceptor sam gtf
    if jct.acceptor_sam.exists:
        query = jct.acceptor_sam.acceptor()
        chrom = jct.acceptor_sam.chromosome
        chrom_key = chrom+jct.acceptor_sam.strand
        if chrom_key in chrom_gtfs_acc and chrom_key in chrom_acc_libs:
            gtfs_acc = chrom_gtfs_acc[chrom_key]
            acc_lib = chrom_acc_libs[chrom_key]
            #closest_acc_ind,its = bin_search_gtf(query,acc_lib)
            closest_acc_ind = brute_search_gtf(query,acc_lib) #<-- RB trying brute force
            closest_results["acceptor"] = gtfs_acc[closest_acc_ind]

    return closest_results


#################################
#    Binary Search Closest GTF  #
#################################
#Recursive binary search function
#NOTE this could probably be sped up using bisectleft, but I doubt by much
def bin_search_gtf(query,library,start_ind=0,end_ind=-1,its=1,disp=False):
    """
    Goal: do a binary search through the library for the closest ind
    Arguments:
        query is an int
        library is a list[int] of genomic positions to look through
        start_ind is defaulted to 0 and keeps track of where to look
        end_ind is defaulted to -1 and keeps track of where to look

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
            sys.stdout.write(str(its)+")"+" Found closest to:["+library[ret_ind]+"]\n")
        return ret_ind,its
    #Recursive case
    else:
        mid_ind = (start_ind+end_ind)/2
        #Determine whether or not to print this line
        #(almost never want to except when debugging)
        if disp:
            print str(its)+")",library[start_ind],"--",library[mid_ind],"--",library[end_ind]
            
        #Check to see how to recurse
        if query < library[mid_ind]:
            return bin_search_gtf(query,library,start_ind,mid_ind,its+1,disp=disp)
        else:
            return bin_search_gtf(query,library,mid_ind,end_ind,its+1,disp=disp) 

#####################################
#    Brute Force Closest GTF        #
#####################################
def brute_search_gtf(query,library):
    """
    Goal: do brute search through the library for the closest ind
    Arguments:
        query is an int
        library is a list[int] of genomic positions to look through

    Returns:
        the index of the closest matching library value to the query
    """

    closest_ind = 0
    closest_dist = abs(library[0]-query)
    for ind,val in enumerate(library[1:]):
        if abs(val-query) < closest_dist:
            closest_ind = ind+1 #<-- have to add 1 since skipping first ind
            closest_dist = abs(val-query)

    return closest_ind

#####################################
#         Identify Fusions          #
#####################################
#Takes junctions that already have gtf info
#If a junction has the following properties call it a 'fusion':
#   IF donor and acceptor sams are at_boundary:
#       IF donor and acceptor are on different chromosomes
#           [yes fusion]
#       ELIF donor and acceptor are on different strands
#       ELIF distance between donor and acceptor > threshold
#           [yes fusion]
#       ELSE
#           [no fusion]
#   ELSE:
#       [no fusion]
#
#Returns a list of junctions that are deemed 'fusions'
def identify_fusions(junctions,constants_dict):
    """
    Goal: take the junctions and find the ones that could be fusions
    Arguments:
        junctions is a list[Junction] objects
        span_cutoff is an optional int or float defining min distance for a fusion
            on the same chromosome
    Returns:
        fusion_jcts as a list[Junction] with the junctions defined as fusions
    """
    span_cutoff = constants_dict["span_cutoff"]
    fusion_jcts = []
    for jct in junctions:
        fusion_type = jct.get_fusion_type()
        if "fusion" in fusion_type and "no_fusion" not in fusion_type:
            fusion_jcts.append(jct)
    return fusion_jcts


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
    Goal: take the donor and acceptor sam and categorize them
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
def collapse_junctions(jcts,full_path_name,constants_dict,group_out_file_name=None):
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
        chrom_1 = str(jct.donor_sam.chromosome)
        chrom_2 = str(jct.acceptor_sam.chromosome)
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
            don = jct.donor_sam.donor()
            acc = jct.acceptor_sam.acceptor()
            found_prev_group = False

            #Only compare jct to other jcts if both
            #don and acc are not None
            if don and acc:
                for prev_group in groupings[chroms]:
                    for prev_jct in prev_group:
                        prev_don = prev_jct.donor_sam.donor()
                        prev_acc = prev_jct.acceptor_sam.acceptor()
                        #If any one of the acceptor/donors are None
                        if not prev_don or not prev_acc:
                            continue
                        if abs(don-prev_don) <= collapse_thresh and abs(acc-prev_acc) <= collapse_thresh:
                            #sys.stderr.write("Found match for:\n")
                            #sys.stderr.write(jct.verbose_fasta_string())
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
    group_out_file = None
    if group_out_file_name:
        group_out_file = open(group_out_file_name,"w")

    for chroms in groupings:
        for group in groupings[chroms]:
            if len(group) <= 1:
                singles += group
            else:
                if group_out_file_name:
                    group_out_file.write("Group info:\n")
                    group_out_file.write("".join([m.verbose_fasta_string() for m in group])+"\n")
                counts = [len(member.bin_pair_group) for member in group]
                max_ind = counts.index(max(counts))
                repr_jct = group[max_ind]
                groups.append(repr_jct)

    if group_out_file_name:
        group_out_file.close()

    return singles,groups
 
   
##########################
#     Badfj3 fusions     #
##########################
def badfj3_fusions(fusion_junctions,constants_dict):

    #Bowtie params
    min_score = constants_dict["splice_finding_min_score"]
    read_gap_score = constants_dict["read_gap_score"]
    ref_gap_score = constants_dict["ref_gap_score"]
    num_threads = constants_dict["num_threads"]
    reference = constants_dict["reference"]

    #Build up file stems
    badfj3_stem = os.path.join(constants_dict["output_dir"],"badfj3_")
    don_fasta = badfj3_stem+"don.fasta"
    acc_fasta = badfj3_stem+"acc.fasta"
    badfj3_mapped = badfj3_stem+"mapped.sam"

    #Open the R1 and R2
    don_fasta_f = open(don_fasta,"w")
    acc_fasta_f = open(acc_fasta,"w")

    #Build up the "paired end" files
    for fusion in fusion_junctions:
        header = ">"+str(fusion.jct_ind)+"\n"
        break_point = fusion.splice_ind()
        don = fusion.consensus[:break_point]
        acc = fusion.consensus[break_point:]
        don_fasta_f.write(header+don+"\n")
        acc_fasta_f.write(header+acc+"\n")

    #Close the R1 and R2
    don_fasta_f.close()
    acc_fasta_f.close()

    badfj3_gap = "500000" #If the don/acc can map within 1/2 Mb, then jct not fusion

    #Run bowtie2 on the R1 and R2
    with open(badfj3_mapped,"w") as badfj3_mapped_f:
        subprocess.call(
            ["bowtie2", "-f", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
             "-x", reference, "-X", badfj3_gap, "--ff", "-1", don_fasta, "-2", acc_fasta], stdout=badfj3_mapped_f)


    still_fusions = []
    now_jcts = []

    #Read back the sam file
    with open(badfj3_mapped,"r") as badfj3_mapped_f:
        for line in badfj3_mapped_f:
            print line #NOTE!!! need to see how PE output looks

    return still_fusions,now_jcts


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


################################################################
#   Print out the constants dict for better version tracking   #
################################################################
def write_constants_dict(constants_dict,params_out_name):
    """
    Goal: print out the constants dict to inform what parameters were used on the run
    Arguments:
        constants_dict is a dict keyed by string and valued by numeric or string
        params_out_name is the file path to store the constants dict info
    Returns:
        None, everything is written to the file
    """
    uniform_len = max([len(str(key)) for key in constants_dict])
    spaces = " "*uniform_len

    with open(params_out_name,"w") as params_out:
        #Header info with date and time
        params_out.write("Parameter file used for SPORK run\n")
        params_out.write("\tDate "+time.strftime("%d/%m/%y")+" (day/month/year)\n")
        params_out.write("\tTime "+time.strftime("%H:%M:%S")+"\n")
        params_out.write("\nParameters:\n")
        params_out.write("-"*(uniform_len+4)+"\n")

        #Loop through the parameters printing nicely
        len_sorted_params = sorted(constants_dict.keys(),key=lambda k: len(k))
        for param in len_sorted_params:
            padded_param = param[:uniform_len]+spaces[:uniform_len-len(param)]
            param_val = str(constants_dict[param])
            params_out.write(padded_param+"    :    "+param_val+"\n")

###################################
#   Track NUP214 as a test case   #
###################################
def follow_nup214(forward_jct,reverse_jct):
    """
    Goal: check if forward or reverse nup214 makes more sense in terms
          of donor and acceptor sites
    Arguments:
        forward_jct is of type Junction
        reverse_jct is of type Junction and is the rev-comp of forward_jct

    Returns:
        nothing, just writes out info to stdout
    """
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


###################################
#   Track NUP214 as a test case   #
###################################
def get_seq_complexity(seq):
    """
    Goal: take in a string and return the sequence complexity. Makes use of zlib to compress
          the string and see how much compression occurred. If lots of compression, then seq was low complexity

    Arguments:
        seq is a sequence (or really any string)
    
    Return:
        a float for the string complexity between 0 and 1 (0 is least complex, 1 is most)
    """
    uncompressed = sys.getsizeof(seq)
    compressed = sys.getsizeof(zlib.compress(seq))
    compression = float(compressed)/uncompressed #0 is worst, 1 is best (visually making 0.675 cutoff)
    return compression

