from SPORK_utils import *

example = "chr17:NXN:862496:+|chr19:C19orf70:5678619:-|no_fusion,num=5,score=0.005,gap=-12,jct_ind=13253	0	0	-	0	5	-	0	0	-	0	0	-	0	0	-	5:0 0	0	1	1	419 -	-	-	-	-	-	-"
report_lines = [example]

reformatted_report_name = "reformatted_report.txt"

with open(reformatted_report_name,"w") as reformatted_report:
    for report_line in report_lines:
        if report_line and report_line[0] != "@":
            split_report_line = report_line.split("\t")
            jct_info = split_report_line[0]
            remaining_info = "\t".join(split_report_line[1:])+"\n"

            split_jct_info = jct_info.split("|")
            chr1,gene1,pos1,strand1 = split_jct_info[0].split(":")
            chr2,gene2,pos2,strand2 = split_jct_info[1].split(":")

            split_extra_info = split_jct_info[2].split(",")
            fusion_type = split_extra_info[0]
            remaining_extra_info = ",".join(split_extra_info[1:])

            pos1 = int(pos1)
            pos2 = int(pos2)

            
