# combine_depth.py
# python combine_depth.py cell1 C dir_bowtie_result
# to combine multiple depth.txt file into one total_depth.txt

import sys
import numpy as np
import os
#import re

# dir_input = sys.argv[1]
# cell_input="cell1"
# flag="C"

cell_input=sys.argv[1]
flag=sys.argv[2]

# dir_input="/home/ygli/gam_paper_data/gam_seq_mapped3/bowtie_result/"+cell_input+"/"
dir_input=sys.argv[3]+cell_input+"/"
# pattern=r'*.chr$'
key_list=["chr1","chr2","chr3","chr4","chr5",
        "chr6","chr7","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15",
        "chr16","chr17","chr18","chr19",
        "chrM","chrX"]

chr_total_reads=[]
chr_total_cover=[]
for root, dirs, files in os.walk(dir_input):
    for file_name in files:
        chr_reads={"chr1":0,"chr2":0,"chr3":0,"chr4":0,"chr5":0,
        "chr6":0,"chr7":0,"chr8":0,"chr9":0,"chr10":0,
        "chr11":0,"chr12":0,"chr13":0,"chr14":0,"chr15":0,
        "chr16":0,"chr17":0,"chr18":0,"chr19":0,
        "chrM":0,"chrX":0}
        chr_cover={"chr1":0,"chr2":0,"chr3":0,"chr4":0,"chr5":0,
        "chr6":0,"chr7":0,"chr8":0,"chr9":0,"chr10":0,
        "chr11":0,"chr12":0,"chr13":0,"chr14":0,"chr15":0,
        "chr16":0,"chr17":0,"chr18":0,"chr19":0,
        "chrM":0,"chrX":0}
        # name="B10F-NP260_S.depth.txt"
        if "_"+flag+".depth.txt" in file_name:
            print(file_name)
            chr_each_reads=[]
            chr_each_reads.append(file_name)
            chr_each_cover=[]
            chr_each_cover.append(file_name)
            with open(os.path.join(root, file_name),'r') as f:
                all_f=f.readlines()
                for each_line in all_f:
                    chr_name,head,tail,reads,base,lenth,cover=each_line.replace("\n","").split("\t")
                    if chr_reads[chr_name]==0:
                        chr_reads[chr_name]=reads
                        chr_cover[chr_name]=cover
            for key in key_list:
                chr_each_reads.append(chr_reads[key])
            for key in key_list:
                chr_each_cover.append(chr_cover[key])
            chr_total_reads.append(chr_each_reads)
            chr_total_cover.append(chr_each_cover)

with open(os.path.join(root, cell_input+"_"+flag+"_chr_reads.list"),'w') as f:
    f.write("chromatin\t")
    for key in key_list:
        f.write(key)
        f.write("\t")
    f.write("\n")
    for each_line in chr_total_reads:
        for i in each_line:
            f.write(str(i)+"\t")
        f.write("\n")

with open(os.path.join(root, cell_input+"_"+flag+"_chr_cover.list"),'w') as f:
    f.write("chromatin\t")
    for key in key_list:
        f.write(key)
        f.write("\t")
    f.write("\n")
    for each_line in chr_total_cover:
        for i in each_line:
            f.write(str(i)+"\t")
        f.write("\n")
