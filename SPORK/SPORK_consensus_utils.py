#Imports
import sys

#build and score consensus function
#(1) pads the left and right side of each sequence depending on where it mapped in the bin
#(2) finds a consensus using a majority vote
#(3) calculates a consensus score as the (number of mismatches)/(total_possible_mismatches)
def build_and_score_consensus(mapped_reads,strand,id_to_seq,bin_size,constants_dict):
    #Add padding to the left sides of the reads based on where they fell in the bin pair
    #Note that for the plus strand, more padding should be added the larger ther pos%bin_size value
    #while this should have less padding on a minus strand read. This is implemented w/ the ternary expression
    read_num_to_read_id = constants_dict["read_num_to_read_id"]
    padded_seqs = []
    left_padded_seqs = []
    for mapped_read in mapped_reads:
        id_key = mapped_read.read_id
        id_key = id_key.replace("/5_prime","")
        id_key = id_key.replace("/3_prime","")
        id_key = read_num_to_read_id[id_key]
        #id_key = id_key.replace("_"," ")
        if id_key in id_to_seq:
            full_seq = id_to_seq[id_key]
        else:
            sys.stderr.write("SPORK ERROR: Couldn't find sequence "+id_key+" in consensus building\n")
            sys.exit(1)
        mapped_read.read_id = id_key
        left_padding = int(mapped_read.start%bin_size) if strand == "+" else int(bin_size-mapped_read.start%bin_size)
        left_padded_seq = " "*left_padding+full_seq

        #Don't allow exact sequence duplicates
        if left_padded_seq not in left_padded_seqs:
            left_padded_seqs.append(left_padded_seq)

    #Add padding on the right sides so that every sequence is the same length
    #Handled the same for plus and minus strand
    max_len_seq = max([len(seq) for seq in left_padded_seqs])
    for left_padded_seq in left_padded_seqs:
        left_padded_len = len(left_padded_seq)
        padding_to_add = max_len_seq-left_padded_len
        padded_seq = left_padded_seq+" "*padding_to_add
        padded_seqs.append(padded_seq)

    #Go through each position and get the majority vote as the consensus base
    #The dictionaries are just to help convert letters into indices and back
    min_bases_per_col = constants_dict["min_bases_per_col"]
    consensus = ""
    empty_spaces = 0
    num_possible_discrepancies = 0
    num_discrepancies = 0
    base_dict = {"A":0,"C":1,"G":2,"T":3}
    rev_base_dict = {0:"A",1:"C",2:"G",3:"T"}
    for seq_ind in range(max_len_seq):
        counts = [0,0,0,0]
        for seq in padded_seqs:
            base = seq[seq_ind].upper()
            if base in base_dict:
                counts[base_dict[base]] += 1
        total_bases = sum(counts)
        if total_bases >= min_bases_per_col:
            num_possible_discrepancies += total_bases
            max_count = max(counts)
            num_discrepancies += total_bases-max_count
            max_index = counts.index(max_count)
            consensus += rev_base_dict[max_index]
        else: empty_spaces += 1

    #Give a terrible score if there are 0 possible discrepancies
    #This arises most often when a junction only has one unique sequence
    if int(num_possible_discrepancies) <= 0:
        consensus_score = 999999
        #sys.stderr.write("Null consensus:\n"+"\n".join(padded_seqs)+"\n")
    else:
        consensus_score = float(num_discrepancies)/int(num_possible_discrepancies)

    #Print out the consensus's in a nice way
    print_consensus = False
    if print_consensus:
        print "Strand: "+strand
        print "Discrepancies: "+str(num_discrepancies)
        print "="*(len(seq)+2)
        for seq in padded_seqs:
            print "|"+seq+"|"
        print "="*(len(seq)+2)
        print "|"+" "*(len(seq)-len(consensus)-min_bases_per_col)+consensus+" "*min_bases_per_col+"|"
        left_score_str = "|Score: "+str(consensus_score)+" "
        right_score_str = " "*(len(seq)+2-len(left_score_str)-1)+"|"
        print left_score_str+right_score_str
        print "="*(len(seq)+2)
        print ""

    return consensus,consensus_score

if __name__ == "__main__":
    print "Running through test cases (once I make some)"
