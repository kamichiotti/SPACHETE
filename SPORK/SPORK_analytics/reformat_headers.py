#Goal is to take headers in the form of:
#>|chr10|ERCC6|50658329-50658362|strand1:+|boundary_dist1:6129|at_boundary1:False|_|chr10|ERCC6|50658396-50658429|strand2:+|boundary_dist2:6095|at_boundary2:False|_|splice:33|score:0.0|num:2|splice:Gapped|
#And convert it to:
#>chr10:ERCC6:50649599:+|chr10:ERCC6:50649645:+|no_fusion

#Imports
import sys
import os

#Setting up paths
spork_output_dir = "/scratch/PI/horence/rob/spork_outputs/"
remaining_path = "ovarian_2_or_more/"
dir_path = spork_output_dir+remaining_path
spork_jcts_name = "novel_junctions.fasta"

#machete_format function
def machete_format(in_file_name,out_file_name):
    in_file = open(in_file_name,"r")
    out_file = open(out_file_name,"w")
    wrote_header = False
    for in_line in in_file:
        if ">" in in_line:
            upstream,downstream,info = in_line.split("|_")
            up_split = upstream.split("|")
            down_split = downstream.split("|")
            chrom1 = up_split[1]
            chrom2 = down_split[1]
            gene1 = up_split[2]
            gene2 = down_split[2]
            pos1 = up_split[3].split("-")[0]
            pos2 = down_split[3].split("-")[0]
            strand1 = up_split[4].split(":")[1]
            strand2 = down_split[4].split(":")[1]

            out_str = ">"
            out_str += chrom1+":"+gene1+":"+pos1+":"+strand1+"|"
            out_str += chrom2+":"+gene2+":"+pos2+":"+strand2+"|"

            if "None" not in out_str:
                wrote_header = True
                fusion = "no_fusion"
                if chrom1 != chrom2:
                    fusion = "fusion"
                elif abs(int(pos1)-int(pos2)) > 1e6:
                    fusion = "fusion"
                elif strand1 != strand2:
                    fusion = "fusion"
                    
                out_str += fusion+"\n"
                out_file.write(out_str)
        else:
            if wrote_header:
                out_file.write(in_line)
                wrote_header = False

    in_file.close()
    out_file.close()

#Finding all sub files that have the spork_jcts_name
possible_dirs = []
sub_dirs = [dir_path]
while len(sub_dirs) > 0:
    curr_dir = sub_dirs[0]
    sub_dirs = sub_dirs[1:]
    sub_dirs += [curr_dir+sub_dir+"/" for sub_dir in os.listdir(curr_dir)
                 if os.path.isdir(curr_dir+sub_dir+"/")]
    possible_dirs.append(curr_dir)

#Making the call to the machete format function
spork_jct_names = [p_dir+f for p_dir in possible_dirs for f in os.listdir(p_dir) if f == spork_jcts_name]
for spork_jct_name in spork_jct_names:
    jct_dir = "/".join(spork_jct_name.split("/")[:-1])
    machete_formatted = jct_dir+"/machete_"+spork_jcts_name
    print machete_formatted
    machete_format(spork_jct_name,machete_formatted)


