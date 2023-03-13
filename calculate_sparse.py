# calculate_sparse.py

# python calculate_sparse.py "Cell11_C_combine.txt" "/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result/combine/"
# to combine the sparsity of the combined txt
# output cover_begin\tcover_ent\tfrequence

import sys
import numpy as np
import os
#import re

# file_input = "Cell11_C_combine.txt"
# dir_input = "/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result/combine/"

file_input=sys.argv[1]
dir_input=sys.argv[2]

# pattern=r'*.chr$'

chr_index_all=[]
cover_num_all=[]

with open(dir_input+file_input,"r") as f:
    head=f.readline().split()
    head_later=head[0]+"\t"+head[1]+"\t"+head[2]
    all_line=f.readlines()
    for each_line in all_line:
        chr_i,chr_begin,chr_end=each_line.split()[0:3]
        chr_index=chr_i+"\t"+chr_begin+"\t"+chr_end
        cover_num=np.sum([x!="0" for x in each_line.replace("\n","").split()[3:]])
        chr_index_all.append(chr_index)
        cover_num_all.append(cover_num)

with open(dir_input+file_input.replace(".txt","_count.txt"),"w") as f:
    # f.write(head_later+"\t"+count)
    for i in range(len(chr_index_all)):
        f.write(chr_index_all[i]+"\t"+str(cover_num_all[i]))
        f.write("\n")










