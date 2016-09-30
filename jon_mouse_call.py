import subprocess
import os

#Bowtie2 Parameters
read_gap_score = "--rdg 50,50"          #read gap set to be very high to not allow gaps
ref_gap_score = "--rfg 50,50"           #read gap set to be very high to not allow gaps
min_score = "--score-min L,0,-0.10"     #minimum allowed Bowtie2 score
num_threads = "8"                       #for multithreading
reference = 
#Input fasta file
input_fasta_name = ""

#Output directory (will create a new directory)
output_dir = ""

#Call Bowtie2 on the upstream primers
subprocess.call(
            ["bowtie2", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
             "--un", os.path.join(output_dir,"unaligned_5prime_"+R_file_name+".fq"), "-x", reference, five_prime_fq_name], stdout=five_prime_mapped)
