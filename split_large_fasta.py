import sys
import os

def split_big_fasta(big_fasta_name):
    path_prefix = big_fasta_name.split("_FarJunctions.fa")[0]
    name_prefix = path_prefix.split("/")[-1]
    name_suffix = "FarJunctions.fa"

    if not os.path.isdir(path_prefix):
        os.mkdir(path_prefix)

    chroms = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10",
              "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
              "chr21","chr22", "chrX", "chrY"]

    #Open all the out files
    out_files = {}
    for chrom in chroms:
        out_files[chrom] = open(path_prefix+"/"+name_prefix+"_"+chrom+name_suffix,"w")

    #Walk through the big fasta copying out the info
    #to the individual chromosomes files
    big_fasta = open(big_fasta_name,"r")
    header_line = big_fasta.readline()
    while header_line:
        seq_line = big_fasta.readline()
        chrom = header_line.split(":")[0][1:]
        if chrom in out_files:
            out_files[chrom].write(header_line)
            out_files[chrom].write(seq_line)

        #Advance to the next fasta
        header_line = big_fasta.readline()

    #Close all the out files    
    for chrom in chroms:
        out_files[chrom].close()



#Calling from command line
#The big_fasta_name will work best if it's an absolute path
if __name__ == "__main__":
    big_fasta_name = sys.argv[1]
    split_big_fasta(big_fasta_name)
