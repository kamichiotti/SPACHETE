#Class to store GTF entries
#TODO give this better documentation

#Imports
import re
import sys

class GTFEntry(object):
    __slots__ = ["chromosome","source","feature",
                 "start","stop","score","strand",
                 "frame","gene_name",
                 "donor","acceptor","span"]

    def __init__(self,gtf_line='chr\tsrc\tfeat\t-1\t-1\t-1\t+\t0\tgene_name "default";'):
        """
        Goal: initialize a GTFEntry
        Arguments:
            a gtf_line from a standard gtf file, optional, if none given default made

        Returns:
            nothing
        """
        split_gtf_line = gtf_line.split("\t")
        self.chromosome = split_gtf_line[0]
        self.source = split_gtf_line[1]
        self.feature = split_gtf_line[2]
        self.start = int(split_gtf_line[3])
        self.stop = int(split_gtf_line[4])
        self.score = split_gtf_line[5]
        self.strand = split_gtf_line[6]
        self.frame = split_gtf_line[7]
        self.donor = self.stop if self.strand == "+" else self.start
        self.acceptor = self.start if self.strand == "+" else self.stop
        self.span = abs(self.donor-self.acceptor)
        group_info = split_gtf_line[8]
        gene_name_pattern = re.compile('gene_name "(.*?)";')
        gene_name = gene_name_pattern.findall(group_info)
        if len(gene_name) == 0:
            gene_id_pattern = re.compile('gene_id "(.*?)";')
            gene_name = gene_id_pattern.findall(group_info)
            if len(gene_name) == 0:
                sys.stdout.write("SPORK ERROR: in gtf init. No found gene_name or gene_id")
                sys.stderr.write("SPORK ERROR: in gtf init. No found gene_name or gene_id")
                sys.exit(1)
        self.gene_name = gene_name[0]

    def __str__(self):
        """
        Goal: yield a string representation of this GTFEntry
        Arguments:
            none

        Returns:
            a string of important info about this GTFEntry
        """
        ret_str = ""
        ret_str += "chromosome: "+self.chromosome+"\t"
        ret_str += "name: "+self.gene_name+"\t"
        ret_str += "donor: "+str(self.donor)+"\t"
        ret_str += "accep: "+str(self.acceptor)+"\t"
        ret_str += "strand: "+self.strand+"\t"
        ret_str += "start: "+str(self.start)+"\t"
        ret_str += "stop: "+str(self.stop)+"\t"
        return ret_str

    def __lt__(self,other):
        """
        Goal: allows comparison between GTFEntries
        Arguments:
            other is also a GTFEntry

        Returns:
            true if self's chromosome is smaller than other's,
            or if they are shared and self's start position is smaller
            false otherwise
        """
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome
        else:
            return self.start < other.start

