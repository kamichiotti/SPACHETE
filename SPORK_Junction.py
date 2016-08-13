#Identified junction class

#Imports
from SPORK_SAMEntry import SAMEntry
import sys

#Junction class
class Junction(object):
    __slots__ = ["consensus","score","bin_pair","bin_pair_group","took_reverse_compliment","constants_dict",
                 "upstream_sam","downstream_sam","jct_ind"]

    def __init__(self,consensus,score,bin_pair_group,took_reverse_compliment,constants_dict):
        """
        Goal: Initialization function of junction
        Arguments:
            consensus                -- str
            score                    -- float
            bin_pair_group           -- list[bin_pair]
            took_reverse_coompliment -- bool
            constants_dict           -- dict[str->multiple types]

        Returns:
            nothing
        """
        #Read in arguments
        self.consensus = consensus
        self.score = score
        self.bin_pair_group = bin_pair_group
        self.took_reverse_compliment = took_reverse_compliment
        self.constants_dict = constants_dict

        #Find chromosome, bin_pair and strand info from the first mapped read
        rep_bin_pair = self.bin_pair_group[0]
        self.bin_pair = rep_bin_pair.bin_pair
        self.upstream_sam = SAMEntry()
        self.downstream_sam = SAMEntry()


    #Use the sam's to find the splice index in reference to the concensus
    def splice_ind(self):
        """
        Goal: Get the splice site in consensus coordinates [0,len(consensus)-1]
        Arguments:
            none

        Returns:
            the 3' edge of the upstream sequence if both sams are defined
            otherwise returns the middle index of the consensus as a guess
        """
        if self.upstream_sam.exists and self.downstream_sam.exists:
            #shouldn't have to do this anymore since I fixed the reads in splice utils
            #donor_offset = int(self.upstream_sam.read_id.split("_")[-1])
            #return donor_offset+len(self.upstream_sam.seq)
            return len(self.upstream_sam.seq)
        else:
            return len(self.consensus)/2


    #Use the sam's again to find the size of the gap between the two pieces
    def splice_gap(self):
        """
        Goal: Find the distance between the upstream and downstream splice sites
              in consensus coordinates [0,len(consensus)-1]
        Arguments:
            none

        Returns:
            the distance between the 3' end of the upstream and 5' end of the downstream
            if one or both of the sam's are undefined return -1
        """
        if self.upstream_sam.exists and self.downstream_sam.exists:
            upstream_pos = self.consensus.index(self.upstream_sam.seq)+len(self.upstream_sam.seq)
            downstream_pos = self.consensus.index(self.downstream_sam.seq)

            upstream_pos = -1
            downstream_pos = -1

            #comp = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
            #rev_upstream_seq = "".join([comp[x] for x in self.upstream_sam.seq])[::-1]
            #rev_downstream_seq = "".join([comp[x] for x in self.downstream_sam.seq])[::-1]

            #don_seq = self.consensus[:self.splice_ind()]
            #acc_seq = self.consensus[self.splice_ind():]
            #sys.stderr.write(don_seq+"|"+acc_seq+"\n")
            #sys.stderr.write(self.upstream_sam.seq+"\n")
            #sys.stderr.write(rev_upstream_seq+"\n")
            #sys.stderr.write(self.downstream_sam.seq+"\n")
            #sys.stderr.write(rev_downstream_seq+"\n\n")
            #Check if mapping needs to have rev comp for upstream
            #if self.upstream_sam.seq in self.consensus:
            #    upstream_pos = self.consensus.index(self.upstream_sam.seq)+len(self.upstream_sam.seq)
            #else:
            #    upstream_pos = self.consensus.index(rev_upstream_seq)+len(self.upstream_sam.seq)

            #Check if mapping needs to have rev comp for downstream
            #if self.downstream_sam.seq in self.consensus:
            #    downstream_pos = self.consensus.index(self.downstream_sam.seq)
            #else:
            #    downstream_pos = self.consensus.index(rev_downstream_seq)
            
            #Return the difference (should be positive for linear)
            return upstream_pos-downstream_pos

        #If at least one of the two is undefined then just return -1
        else:
            return None


    #Use the sam's again to find the size of the gap between the two pieces
    def span(self):
        """
        Goal: find the genomic span between the sams
        Arguments:
            none

        Returns:
            if both exist subtract the stop of the downstream
            from the upstream start position (can be negative)
            if one or both don't exist just return -1
        """
        if self.upstream_sam.exists and self.downstream_sam.exists:
            return self.upstream_sam.start-self.downstream_sam.stop
        else:
            return -1


    #Give a name to the splice type for this junction
    def splice_type(self):
        """
        Goal: get the type of splice this junction represents
        Arguments:
            none

        Returns:
            "Full" if both sams exist and have zero gaps in the split
            "Gapped" if both sams exist but there is space in the middle
            "Five_Only" if only the upstream sam exists
            "Three_Only" if only the downstream sam exists
            "None" if niether sam exists
        """
        if self.upstream_sam.exists and self.downstream_sam.exists:
            if self.splice_gap() == 0:
                return "Full"
            else:
                return "Gapped"
        elif self.upstream_sam.exists:
            return "Five_Only"
        elif self.downstream_sam.exists:
            return "Three_Only"
        else:
            return "None"


    #Check to see if this jct represents a fusion
    def check_fusion(self,span_cutoff=1e6):
        """
        Goal: check if this junction represents a fusion
        Arguments:
            none

        Returns:
            bool of whether or not the upstream and downstream have different genes
            if one or more don't exists then return False
        """
        if self.at_boundary("upstream") and self.at_boundary("downstream"):
            if self.upstream_sam.chromosome != self.downstream_sam.chromosome:
                return True
            elif self.upstream_sam.strand != self.downstream_sam.strand:
                return True
            elif abs(self.span()) > span_cutoff:
                return True

        #If none of the previous ifs are hit, return false
        return False


    #Get distance to closest splice boundary
    def boundary_dist(self,stream):
        """
        Goal: get the distance of the specified stream from the closest exon
        Arguments:
            stream which is a string and can be either "upstream" or "downstream"

        Returns:
            the min distance to any exon of the specified stream
            if the specified stream is not specified, then returns -1
        """
        sam = self.upstream_sam if stream == "upstream" else self.downstream_sam
        if sam.gtf:
            boundary_start = sam.gtf.start
            boundary_stop = sam.gtf.stop
            sam_pos = sam.stop if stream == "upstream" else sam.start
            start_dist = abs(sam_pos-boundary_start)
            stop_dist = abs(sam_pos-boundary_stop)
            return min(start_dist,stop_dist)
        else:
            return -1

    #Return whether or not an upstream and downstream is at a boundary
    def at_boundary(self,stream,radius=8):
        """
        Goal: check to see if the specified sam is at an exon boundary
        Arguments:
            stream of type string. should be "upstream" or "downstream"
            to specify which sam to check

            radius is optional and signifies the maximum distance from
            an exon boundary to consider a sam. Default is 3

        Returns:
            a boolean of whether or not the specified sam is within
            'radius' distance of any exon boundary
        """
        dist = self.boundary_dist(stream)
        if 0 <= dist <= radius:
            return True
        else:
            return False

    #Returns whether or not this junction is linear
    def linear(self):
        """
        Goal: check to see if this junction in linear
        Arguments:
            none

        Returns:
            a boolean of whether the junction is linear or not
        """
        five_prime_bin,three_prime_bin,strand_info = self.bin_pair.split("_")
        five_prime_chr = five_prime_bin.split(":")[0]
        five_prime_bin = five_prime_bin.split(":")[1]
        three_prime_chr = three_prime_bin.split(":")[0]
        three_prime_bin = three_prime_bin.split(":")[1]
        linear = True if int(five_prime_bin) <= int(three_prime_bin) else False
        linear = not linear if self.took_reverse_compliment else linear
        return linear

    #Format the junction for MACHETE in fasta form
    #NOTE only call this function on 'fusion' identified junctions
    def fasta_MACHETE(self):
        """
        Goal: produce a fasta_string for MACHETE
        Arguments:
            none
        Returns:
            a fasta formatted string (with a newline between header and sequence)
        """
        #Make the necessary variables
        chrom1 = self.upstream_sam.chromosome
        chrom2 = self.downstream_sam.chromosome
        genes1 = self.upstream_sam.str_gene()
        genes2 = self.downstream_sam.str_gene()
        pos1 = self.upstream_sam.stop
        pos2 = self.downstream_sam.start
        strand1 = self.upstream_sam.strand
        strand2 = self.downstream_sam.strand
        fusion = "fusion" if self.check_fusion() else "no_fusion"

        #Start building the fasta string
        fasta_str = ""
        fasta_str += ">"
        fasta_str += str(chrom1)+":"+str(genes1)+":"+str(pos1)+":"+str(strand1)+"|"
        fasta_str += str(chrom2)+":"+str(genes2)+":"+str(pos2)+":"+str(strand2)+"|"
        fasta_str += fusion
        fasta_str += ",num="+str(len(self.bin_pair_group))
        fasta_str += ",score="+str(self.score)
        fasta_str += ",gap="+str(self.splice_gap())
        fasta_str += "\n"

        #Add the actual padded consensus to the output string
        splice_flank_len = int(self.constants_dict["splice_flank_len"])
        full_consensus = self.format_consensus(splice_flank_len)
        fasta_str += str(full_consensus)+"\n"

        return fasta_str


    #Format the junction to print in fasta form
    def fasta_string(self,jct_ind=False):
        """
        Goal: produce a fasta_string
        Arguments:
            optionally include a junction index.
            if it is included, it will be printed out

        Returns:
            a fasta formatted string (with a newline between header and sequence)
        """
        fasta_str = ""
        fasta_str += ">|"+str(self.upstream_sam.chromosome)+"|"
        fasta_str += str(self.upstream_sam.str_gene())+"|"
        fasta_str += str(self.upstream_sam.start)+"-"
        fasta_str += str(self.upstream_sam.stop)+"|"
        fasta_str += "strand1:"+str(self.upstream_sam.strand)+"|"
        fasta_str += "boundary_dist1:"+str(self.boundary_dist("upstream"))+"|"
        fasta_str += "at_boundary1:"+str(self.at_boundary("upstream"))+"|_"

        fasta_str += "|"+str(self.downstream_sam.chromosome)+"|"
        fasta_str += str(self.downstream_sam.str_gene())+"|"
        fasta_str += str(self.downstream_sam.start)+"-"
        fasta_str += str(self.downstream_sam.stop)+"|"
        fasta_str += "strand2:"+str(self.downstream_sam.strand)+"|"
        fasta_str += "boundary_dist2:"+str(self.boundary_dist("downstream"))+"|"
        fasta_str += "at_boundary2:"+str(self.at_boundary("downstream"))+"|_"

        fasta_str += "|splice:"+str(self.splice_ind())+"|"
        fasta_str += "score:"+str(self.score)+"|"
        #fasta_str += "fusion:"+str(self.check_fusion())+"|" #seems like unnecessary info at this point
        fasta_str += "num:"+str(len(self.bin_pair_group))+"|"
        fasta_str += "splice:"+str(self.splice_type())+"|"
        fasta_str += "jct_ind:"+str(jct_ind)+"|\n" if jct_ind else "\n"

        splice_flank_len = int(self.constants_dict["splice_flank_len"])
        full_consensus = self.format_consensus(splice_flank_len)
        fasta_str += str(full_consensus)+"\n"
        fasta_str += str(self.upstream_sam.seq)+"\n"
        fasta_str += " "*self.splice_ind()+str(self.downstream_sam.seq)+"\n"

        #Also printing out gtf information
        #fasta_str += "Upstream_gtf:"+str(self.upstream_sam.gtf)+"\n"
        #fasta_str += "Downstream_gtf:"+str(self.downstream_sam.gtf)+"\n"
        return fasta_str


    #Format the junction to print in fasta form
    def verbose_fasta_string(self):
        """
        Goal: produce a fasta formatted string of this junction with lots of header info
        Arguments:
            none

        Returns:
            a fasta string (with a newline between the header and sequence)
        """
        fasta_str = ""
        fasta_str += ">|chromosome1:"+str(self.upstream_sam.chromosome)+"|"
        fasta_str += "genes1:"+str(self.upstream_sam.str_gene())+"|"
        fasta_str += "start1:"+str(self.upstream_sam.start)+"|"
        fasta_str += "stop1:"+str(self.upstream_sam.stop)+"|"
        fasta_str += "strand1:"+str(self.upstream_sam.strand)+"|"
        fasta_str += "boundary_dist1:"+str(self.boundary_dist("upstream"))+"|"
        fasta_str += "at_boundary1:"+str(self.at_boundary("upstream"))+"|_"

        fasta_str += "|chromosome2:"+str(self.downstream_sam.chromosome)+"|"
        fasta_str += "genes2:"+str(self.downstream_sam.str_gene())+"|"
        fasta_str += "start2:"+str(self.downstream_sam.start)+"|"
        fasta_str += "stop2:"+str(self.downstream_sam.stop)+"|"
        fasta_str += "strand2:"+str(self.downstream_sam.strand)+"|"
        fasta_str += "boundary_dist2:"+str(self.boundary_dist("downstream"))+"|"
        fasta_str += "at_boundary2:"+str(self.at_boundary("downstream"))+"|_|"

        fasta_str += "splice:"+str(self.splice_ind())+"|"
        fasta_str += "span:"+str(self.span())+"|"
        fasta_str += "score:"+str(self.score)+"|"
        fasta_str += "fusion:"+str(self.check_fusion())+"|"
        fasta_str += "num:"+str(len(self.bin_pair_group))+"|"
        fasta_str += "splice-gap:"+str(self.splice_gap())+"|"
        fasta_str += "splice-type:"+str(self.splice_type())+"|"
        fasta_str += "took-rev-comp:"+str(self.took_reverse_compliment)+"|\n"

        # Add N padding to the consensus to get a uniform len
        splice_flank_len = int(self.constants_dict["splice_flank_len"])
        full_consensus = self.format_consensus(splice_flank_len)
        fasta_str += str(full_consensus)+"\n"
        return fasta_str

    
    #Add N padding to the consensus to get a uniform len
    #With the splice site in the middle
    def format_consensus(self,splice_flank_len):
        """
        Goal: return the consensus properly formatted centered and uniform len
        Arguments:
            splice_flank_len is an int deciding how long either side should be
            from the consensus
        Returns:
            a string of either the full consensus of None if there is no splice ind
        """
        full_consensus = None
        if self.splice_ind() != -1:
            splice_flank_len = int(self.constants_dict["splice_flank_len"])
            left_padding = "N"*(splice_flank_len-self.splice_ind())
            right_padding = "N"*(splice_flank_len-(len(self.consensus)-self.splice_ind()))
            five_consensus = self.consensus[:self.splice_ind()]
            three_consensus = self.consensus[self.splice_ind():]
            full_consensus = left_padding+five_consensus+three_consensus+right_padding
        return str(full_consensus)
   

    #Give back the R1 readIDs used to make this junction
    def get_read_ids(self):
        """
        Goal: return a list of the read ids (strings) that made this junction
        Arguments:
            none
        Returns:
            a list[string] of the read-ids for this junction
        """
        read_ids = []
        for bin_pair in self.bin_pair_group:
            upstream_id = bin_pair.five_prime_SAM.read_id.replace("/5_prime","")
            downstream_id = bin_pair.five_prime_SAM.read_id.replace("/3_prime","")
            if upstream_id == downstream_id:
                read_ids.append(upstream_id)
            else:
                sys.stderr.write("ERROR, nonmatching ids in jct: ["+upstream_id+"] vs ["+downstream_id+"]\n")
                sys.exit(1)

        return read_ids

    #More human readable format
    def __str__(self):
        """
        Goal: output the junction in an expanded human readable form
        Arguments:
            none

        Returns:
            the string to be printed out
        """
        out_str = ""
        out_str += "Junction with bin pair ["+self.bin_pair+"] with ["+str(len(self.bin_pair_group))+"] reads mapped\n"
        out_str += "Linear " if self.linear() else "Non-Linear "
        out_str += "Upstream on the "+str(self.upstream_sam.strand)+" strand and downstream on the "+str(self.downstream_sam.strand)+"\n"
        out_str += "5' map position ["+str(self.upstream_sam.start)+"-"+str(self.upstream_sam.stop)+"]\n"
        out_str += "3' map position ["+str(self.downstream_sam.start)+"-"+str(self.downstream_sam.stop)+"]\n"
        out_str += "Consensus with score ["+str(self.score)+"] and upstream splice site ["+str(self.upstream_sam.stop)+"]:\n"
        out_str += str(self.consensus)+"\n"
        out_str += str(self.upstream_sam.seq)+"\n"
        out_str += " "*len(str(self.upstream_sam.seq))+str(self.downstream_sam.seq)+"\n"
        out_str += "Upstream genes ["+str(self.upstream_sam.str_gene())+"]\n"
        out_str += "Downstream genes ["+str(self.downstream_sam.str_gene())+"]\n"
        return out_str

    #Rank junctions in order of bin_pairs when sorted
    def __lt__(self,other):
        """
        Goal: give a comparison operator for the Junction class
        Arguments:
            other junction to compare to

        Returns:
            a boolean of whether or not this bin_pair
            is smaller than the other bin_pair
        """
        return self.bin_pair < other.bin_pair


