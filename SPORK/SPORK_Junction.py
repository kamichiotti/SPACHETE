#Identified junction class

#Imports
from SPORK_SAMEntry import SAMEntry
import copy
import sys

#Junction class
class Junction(object):
    __slots__ = ["consensus","score","bin_pair","bin_pair_group","took_reverse_compliment","constants_dict",
                 "donor_sam","acceptor_sam","mapq","jct_ind"]

    def __init__(self,consensus,score,bin_pair_group,jct_ind,took_reverse_compliment,constants_dict):
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
        self.jct_ind = jct_ind
        self.took_reverse_compliment = took_reverse_compliment
        self.constants_dict = constants_dict
        self.mapq = 0

        #Find chromosome, bin_pair and strand info from the first mapped read
        rep_bin_pair = self.bin_pair_group[0]
        self.bin_pair = rep_bin_pair.bin_pair
        self.donor_sam = SAMEntry()
        self.acceptor_sam = SAMEntry()

        #Get some information for the donor sam from the five_prime sam of the bin pair
        self.donor_sam.chromosome = rep_bin_pair.five_prime_chr
        self.donor_sam.start = rep_bin_pair.five_prime_SAM.start
        self.donor_sam.stop = rep_bin_pair.five_prime_SAM.stop

        #Get some information for the acceptor sam from the three_prime sam of the bin pair
        self.acceptor_sam.chromosome = rep_bin_pair.three_prime_chr
        self.acceptor_sam.start = rep_bin_pair.three_prime_SAM.start
        self.acceptor_sam.stop = rep_bin_pair.three_prime_SAM.stop



    #Use the sam's to find the splice index in reference to the concensus
    def splice_ind(self):
        """
        Goal: Get the splice site in consensus coordinates [0,len(consensus)-1]
        Arguments:
            none

        Returns:
            the 3' edge of the donor sequence if both sams are defined
            otherwise returns the middle index of the consensus as a guess
        """
        if self.donor_sam.exists and self.acceptor_sam.exists:
            #NOTE currently doesn't handle gaps well (just returns the donor side index of gap)
            return len(self.donor_sam.seq)
        else:
            return len(self.consensus)/2


    #Use the sam's again to find the size of the gap between the two pieces
    def splice_gap(self):
        """
        Goal: Find the distance between the donor and acceptor splice sites
              in consensus coordinates [0,len(consensus)-1]
        Arguments:
            none

        Returns:
            the distance between the 3' end of the donor and 5' end of the acceptor
            if one or both of the sam's are undefined return -1
        """
        if self.donor_sam.exists and self.acceptor_sam.exists:
            donor_pos = self.consensus.index(self.donor_sam.seq)+len(self.donor_sam.seq)
            acceptor_pos = self.consensus.index(self.acceptor_sam.seq)
            return donor_pos-acceptor_pos
        else:
            return None


    #Use the sam's again to find the size of the gap between the two pieces
    def span(self):
        """
        Goal: find the genomic span between the sams
        Arguments:
            none

        Returns:
            if both exist subtract the donor and acceptor sites
            this distance will always be positive
            if one or both don't exist just return -1
        """
        if self.donor_sam.exists and self.acceptor_sam.exists:
            span = self.donor_sam.donor()-self.acceptor_sam.acceptor()
            return abs(span)
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
            "Five_Only" if only the donor sam exists
            "Three_Only" if only the acceptor sam exists
            "None" if niether sam exists
        """
        if self.donor_sam.exists and self.acceptor_sam.exists:
            if self.splice_gap() == 0:
                return "Full"
            else:
                return "Gapped"
        elif self.donor_sam.exists:
            return "Five_Only"
        elif self.acceptor_sam.exists:
            return "Three_Only"
        else:
            return "None"


    #Check to see if this jct represents a fusion
    def get_fusion_type(self,span_cutoff=1e6):
        """
        Goal: check if this junction represents a fusion
        Arguments:
            none

        Returns:
            bool of whether or not the donor and acceptor have different genes
            if one or more don't exists then return False
        """
        anonat = "" #Can be 'bot', 'donor', 'acceptor', or 'none'
        chroms = "" #Can be 'interchrom', 'distant-intrachrom', or 'local-intrachrom'
        strand = "" #Can be 'inversion', 'plus', or 'minus'
        revreg = "" #Can be 'rev', 'reg', or 'invert'
        
        #Get the anonat type
        if self.at_boundary("donor") and self.at_boundary("acceptor"):
            anonat = "both"
        elif self.at_boundary("donor"):
            anonat = "donor"
        elif self.at_boundary("acceptor"):
            anonat = "acceptor"
        else:
            anonat = "niether"

        #Get the chromosomes type
        if self.donor_sam.chromosome != self.acceptor_sam.chromosome:
            chroms = "interchrom"
        elif self.span() >= span_cutoff:
            chroms = "distant-intrachrom"
        else:
            chroms = "local-intrachrom"

        #Get the strand type
        if self.donor_sam.strand != self.acceptor_sam.strand:
            strand = "inversion"
        elif self.donor_sam.strand == "+":
            strand = "plus"
        elif self.donor_sam.strand == "-":
            strand = "minus"

        #Get the revreg type
        if strand == "inversion":
            revreg = "invert"
        elif self.donor_sam.donor < self.acceptor_sam.acceptor and self.donor_sam.strand == "+":
            revreg = "reg"
        elif self.donor_sam.donor > self.acceptor_sam.acceptor and self.donor_sam.strand == "-":
            revreg = "reg"
        else:
            revreg = "rev"
            
        #Should this be considered a fusion?
        fusion = "no_fusion"

        if anonat == "both":
            if self.splice_gap() != None and abs(self.splice_gap()) <= self.constants_dict["fusion_max_gap"]:
                if chroms != "local-intrachrom":
                    fusion = "fusion"

        #Concatenate them into one string
        fusion_type = fusion+"-"+anonat+"_"+chroms+"_"+strand+"_"+revreg
        return fusion_type


    #Get distance to closest splice boundary
    def boundary_dist(self,splice_site,bowtie_style=True):
        """
        Goal: get the distance of the specified splice site from the closest exon
        Arguments:
            splice_site which is a string and can be either "donor" or "acceptor"
            bowtie_style is an optional boolean argument
                if True (default), then if a donor/acceptor falls to the 'right' of the gtf-site,
                regardless of strand, it will be a positive distance

                if False, then the strand does matter, and being 5' of gtf is negative and 3' is positive


        Returns:
            the distance to the closest gtf of the specified splice_site
            the distance being positive or negative means different things based on the bowtie_style parameter
            explained above
        """
        #If donor distance is requested
        if splice_site == "donor" and self.donor_sam.gtf:
            donor_dist = 0
            if not bowtie_style:
                if self.donor_sam.strand == "+":
                    donor_dist = self.donor_sam.donor()-self.donor_sam.gtf.donor
                elif self.donor_sam.strand == "-":
                    donor_dist = self.donor_sam.gtf.donor-self.donor_sam.donor()
                else:
                    sys.stderr.write("SPORK ERROR: in Junction boundary dist, incorrect strand option \n")
                    sys.exit(1)
            elif bowtie_style:
                donor_dist = self.donor_sam.donor()-self.donor_sam.gtf.donor

            return donor_dist

        #If acceptor distance is requested
        elif splice_site == "acceptor" and self.acceptor_sam.gtf:
            acceptor_dist = 0
            if not bowtie_style:
                if self.acceptor_sam.strand == "+":
                    acceptor_dist = self.acceptor_sam.acceptor()-self.acceptor_sam.gtf.acceptor
                elif self.acceptor_sam.strand == "-":
                    acceptor_dist = self.acceptor_sam.gtf.acceptor-self.acceptor_sam.acceptor()
                else:
                    sys.stderr.write("SPORK ERROR: in Junction boundary dist, incorrect strand option \n")
                    sys.exit(1)
            elif bowtie_style:
                acceptor_dist = self.acceptor_sam.acceptor()-self.acceptor_sam.gtf.acceptor
                
            return acceptor_dist

        #If a different string was passed in or the specified gtf doesn't exist
        else:
            sys.stderr.write("SPORK ERROR: in Junction boundary dist, incorrect str or gtf doesn't exist\n")
            sys.exit(1)

    #Return whether or not an donor and acceptor is at a boundary
    def at_boundary(self,splice_site):
        """
        Goal: check to see if the specified sam is at an exon boundary
        Arguments:
            splice_site of type string. should be "donor" or "acceptor"
            to specify which sam to check

            radius is optional and signifies the maximum distance from
            an exon boundary to consider a sam. Default is 3

        Returns:
            a boolean of whether or not the specified sam is within
            'radius' distance of any exon boundary
        """
        dist = self.boundary_dist(splice_site)
        if abs(dist) <= self.constants_dict["at_boundary_cutoff"]:
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

    #Returns this junction and a reverse compliment of this junction
    #to facilitate finding the gtf's of each and seeing which form is better
    def yield_forward_and_reverse(self):
        """
        Goal: return a copy of self and a reverse compliment of self
        Arguments:
            none

        Returns:
            a tuple of Junction where the first entry is self and the
            second is a reverse compliment of self
        """
        #sys.stdout.write("Before copy in yield_forward_and_reverse\n")
        rev_self = Junction(self.consensus,self.score,self.bin_pair_group,self.jct_ind,self.took_reverse_compliment,self.constants_dict)
        #rev_self = copy.deepcopy(self)
        #sys.stdout.write("After copy in yield_forward_and_reverse\n")
        rev_self.took_reverse_compliment = not rev_self.took_reverse_compliment

        comp = {"A":"T","a":"t","T":"A","t":"a",
                "G":"C","g":"c","C":"G","c":"g",
                "N":"N","n":"n"}

        #Take the reverse compliments of the seqs and switch them between donor and acceptor
        rev_self.consensus = "".join([comp[base] for base in self.consensus])[::-1]
        rev_self.donor_sam.seq = "".join([comp[base] for base in self.donor_sam.seq])[::-1]
        rev_self.acceptor_sam.seq = "".join([comp[base] for base in self.acceptor_sam.seq])[::-1]
        rev_self.donor_sam.seq,rev_self.acceptor_sam.seq = rev_self.acceptor_sam.seq,rev_self.donor_sam.seq
        
        #Flip the strands of both SAMs
        #NOTE only switch the strands if both are + or -, don't do it otherwise
        #Interesting that it works this way, but I drew it out and I'm confident
        rev_self.donor_sam.strand = self.donor_sam.strand
        rev_self.acceptor_sam.strand = self.acceptor_sam.strand
        if rev_self.donor_sam.strand == rev_self.acceptor_sam.strand:
            rev_self.donor_sam.strand = "-" if rev_self.donor_sam.strand == "+" else "+"
            rev_self.acceptor_sam.strand = "-" if rev_self.acceptor_sam.strand == "+" else "+"

        #Trade starts and stops of donor and acceptor and chromosome
        rev_self.donor_sam.start,rev_self.acceptor_sam.start = self.acceptor_sam.start,self.donor_sam.start
        rev_self.donor_sam.stop,rev_self.acceptor_sam.stop = self.acceptor_sam.stop,self.donor_sam.stop
        rev_self.donor_sam.chromosome,rev_self.acceptor_sam.chromosome = self.acceptor_sam.chromosome,self.donor_sam.chromosome
        rev_self.donor_sam.exists = True
        rev_self.acceptor_sam.exists = True

        return self,rev_self

    #Returns this junction and a reverse compliment of this junction
    #to facilitate finding the gtf's of each and seeing which form is better
    def yield_reverse(self):
        """
        Goal: return a reversed self
        Arguments:
            none

        Returns:
            a Junction which is the reverse of self (note does change original)
        """
        self.took_reverse_compliment = not self.took_reverse_compliment

        comp = {"A":"T","a":"t","T":"A","t":"a",
                "G":"C","g":"c","C":"G","c":"g",
                "N":"N","n":"n"}

        #Take the reverse compliments of the seqs and switch them between donor and acceptor
        self.consensus = "".join([comp[base] for base in self.consensus])[::-1]
        self.donor_sam.seq = "".join([comp[base] for base in self.donor_sam.seq])[::-1]
        self.acceptor_sam.seq = "".join([comp[base] for base in self.acceptor_sam.seq])[::-1]
        self.donor_sam.seq,self.acceptor_sam.seq = self.acceptor_sam.seq,self.donor_sam.seq
        
        #Flip the strands of both SAMs
        #NOTE only switch the strands if both are + or -, don't do it otherwise
        #Interesting that it works this way, but I drew it out and I'm confident
        if self.donor_sam.strand == self.acceptor_sam.strand:
            self.donor_sam.strand = "-" if self.donor_sam.strand == "+" else "+"
            self.acceptor_sam.strand = "-" if self.acceptor_sam.strand == "+" else "+"

        #Trade starts and stops of donor and acceptor and chromosome
        self.donor_sam.start,self.acceptor_sam.start = self.acceptor_sam.start,self.donor_sam.start
        self.donor_sam.stop,self.acceptor_sam.stop = self.acceptor_sam.stop,self.donor_sam.stop
        self.donor_sam.chromosome,self.acceptor_sam.chromosome = self.acceptor_sam.chromosome,self.donor_sam.chromosome

        return self


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
        chrom1 = self.donor_sam.chromosome
        chrom2 = self.acceptor_sam.chromosome
        genes1 = self.donor_sam.str_gene()
        genes2 = self.acceptor_sam.str_gene()
        pos1 = self.donor_sam.donor()
        pos2 = self.acceptor_sam.acceptor()
        strand1 = self.donor_sam.strand
        strand2 = self.acceptor_sam.strand
        fusion_type = self.get_fusion_type()

        #Start building the fasta string
        fasta_str = ""
        fasta_str += ">"
        fasta_str += str(chrom1)+":"+str(genes1)+":"+str(pos1)+":"+str(strand1)+"|"
        fasta_str += str(chrom2)+":"+str(genes2)+":"+str(pos2)+":"+str(strand2)+"|"
        fasta_str += fusion_type
        fasta_str += ",num="+str(len(self.bin_pair_group))
        fasta_str += ",score="+str(self.score)
        fasta_str += ",gap="+str(self.splice_gap())
        fasta_str += ",don-dist:"+str(self.boundary_dist("donor"))
        fasta_str += ",acc-dist:"+str(self.boundary_dist("acceptor"))
        fasta_str += ",mapq="+str(self.mapq)
        fasta_str += ",jct_ind="+str(self.jct_ind)
        fasta_str += "\n"

        #Add the actual padded consensus to the output string
        splice_flank_len = int(self.constants_dict["splice_flank_len"])
        full_consensus = self.format_consensus(splice_flank_len)
        fasta_str += str(full_consensus)+"\n"

        return fasta_str


    #Format the junction to print in fasta-esque form
    def log_string(self):
        """
        Goal: produce a fasta_string
        Arguments:
            optionally include a junction index.
            if it is included, it will be printed out

        Returns:
            a description of the junction over multiple lines
        """
        fasta_str = ""
        fasta_str += ">|"+str(self.donor_sam.chromosome)+"|"
        fasta_str += str(self.donor_sam.str_gene())+" "
        fasta_str += str(self.donor_sam.gene_strand())+" strand|"
        fasta_str += str(self.donor_sam.start)+"-"
        fasta_str += str(self.donor_sam.stop)+"|"
        fasta_str += "strand1:"+str(self.donor_sam.strand)+"|"
        fasta_str += "boundary_dist1:"+str(self.boundary_dist("donor"))+"|"
        fasta_str += "at_boundary1:"+str(self.at_boundary("donor"))+"|\n"

        fasta_str += ">|"+str(self.acceptor_sam.chromosome)+"|"
        fasta_str += str(self.acceptor_sam.str_gene())+" "
        fasta_str += str(self.acceptor_sam.gene_strand())+" strand|"
        fasta_str += str(self.acceptor_sam.start)+"-"
        fasta_str += str(self.acceptor_sam.stop)+"|"
        fasta_str += "strand2:"+str(self.acceptor_sam.strand)+"|"
        fasta_str += "boundary_dist2:"+str(self.boundary_dist("acceptor"))+"|"
        fasta_str += "at_boundary2:"+str(self.at_boundary("acceptor"))+"|\n"

        fasta_str += ">|splice:"+str(self.splice_ind())+"|"
        fasta_str += "score:"+str(self.score)+"|"
        fasta_str += "fusion:"+str(self.get_fusion_type())+"|"
        fasta_str += "num:"+str(len(self.bin_pair_group))+"|"
        fasta_str += "splice:"+str(self.splice_type())+"|"
        fasta_str += "mapq="+str(self.mapq)+"|"
        fasta_str += "jct_ind:"+str(self.jct_ind)+"|\n"

        splice_flank_len = int(self.constants_dict["splice_flank_len"])
        full_consensus = self.format_consensus(splice_flank_len)
        fasta_str += str(full_consensus)+"\n"
        fasta_str += str(self.donor_sam.seq)+"\n"
        fasta_str += " "*self.splice_ind()+str(self.acceptor_sam.seq)+"\n"

        #Also printing out gtf information
        #fasta_str += "Donor_gtf:"+str(self.donor_sam.gtf)+"\n"
        #fasta_str += "Acceptor_gtf:"+str(self.acceptor_sam.gtf)+"\n"
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
        fasta_str += ">|chromosome1:"+str(self.donor_sam.chromosome)+"|"
        fasta_str += "genes1:"+str(self.donor_sam.str_gene())+"|"
        fasta_str += "start1:"+str(self.donor_sam.start)+"|"
        fasta_str += "stop1:"+str(self.donor_sam.stop)+"|"
        fasta_str += "strand1:"+str(self.donor_sam.strand)+"|"
        fasta_str += "boundary_dist1:"+str(self.boundary_dist("donor"))+"|"
        fasta_str += "at_boundary1:"+str(self.at_boundary("donor"))+"|_"

        fasta_str += "|chromosome2:"+str(self.acceptor_sam.chromosome)+"|"
        fasta_str += "genes2:"+str(self.acceptor_sam.str_gene())+"|"
        fasta_str += "start2:"+str(self.acceptor_sam.start)+"|"
        fasta_str += "stop2:"+str(self.acceptor_sam.stop)+"|"
        fasta_str += "strand2:"+str(self.acceptor_sam.strand)+"|"
        fasta_str += "boundary_dist2:"+str(self.boundary_dist("acceptor"))+"|"
        fasta_str += "at_boundary2:"+str(self.at_boundary("acceptor"))+"|_|"

        fasta_str += "jct_ind:"+str(self.jct_ind)+"|"
        fasta_str += "splice:"+str(self.splice_ind())+"|"
        fasta_str += "span:"+str(self.span())+"|"
        fasta_str += "score:"+str(self.score)+"|"
        fasta_str += "fusion:"+str(self.get_fusion_type())+"|"
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
            donor_id = bin_pair.five_prime_SAM.read_id.replace("/5_prime","")
            acceptor_id = bin_pair.five_prime_SAM.read_id.replace("/3_prime","")
            if donor_id == acceptor_id:
                read_ids.append(donor_id)
            else:
                sys.stderr.write("SPORK ERROR, nonmatching ids in jct: ["+donor_id+"] vs ["+acceptor_id+"]\n")
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
        out_str += "Donor on the "+str(self.donor_sam.strand)+" strand and acceptor on the "+str(self.acceptor_sam.strand)+"\n"
        out_str += "5' map position ["+str(self.donor_sam.start)+"-"+str(self.donor_sam.stop)+"]\n"
        out_str += "3' map position ["+str(self.acceptor_sam.start)+"-"+str(self.acceptor_sam.stop)+"]\n"
        out_str += "Consensus with score ["+str(self.score)+"] and donor splice site ["+str(self.donor_sam.stop)+"]:\n"
        out_str += str(self.consensus)+"\n"
        out_str += str(self.donor_sam.seq)+"\n"
        out_str += " "*len(str(self.donor_sam.seq))+str(self.acceptor_sam.seq)+"\n"
        out_str += "Donor genes ["+str(self.donor_sam.str_gene())+"]\n"
        out_str += "Acceptor genes ["+str(self.acceptor_sam.str_gene())+"]\n"
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


