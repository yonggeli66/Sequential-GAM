# build_3D_structure_snpstrain_with_parameter.py
# this is from build_3D_structure_snpstrain.py

from __future__ import division
from __future__ import print_function
from mpl_toolkits import mplot3d
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import math
import random
import copy
from scipy.optimize import minimize
import sys
import seaborn as sns; sns.set()
import pandas as pd
import os
import argparse

CUDA_VISIBLE_DEVICES=0

parser = argparse.ArgumentParser(description='parameter for sequential-GAM')
parser.add_argument('--chr_name', type=str, default = "chr1")
parser.add_argument('--cell_name', type=str, default= "6")
parser.add_argument('--method', type=str, default="overlap")
parser.add_argument('--my_step', type=int, default=180000)
parser.add_argument('--my_lr', type=float, default=1e-4)
parser.add_argument('--strain_plot', type=str, default="C")
parser.add_argument('--limit_LAD', type=float, default=1e-7)
parser.add_argument('--limit_res', type=float, default=1e+9)
parser.add_argument('--constant_LAD', type=float, default=1e-8)
parser.add_argument('--constant_res', type=float, default=1e-3)
parser.add_argument('--limit_power_LAD', type=int, default=3)
parser.add_argument('--limit_power_res', type=int, default=3)
parser.add_argument('--save_file_prefix',type=str,default="CS_depart_parameter_all_chr_10kb_")

args = parser.parse_args()

chr_name=args.chr_name
cell_name=args.cell_name
strain_plot=args.strain_plot
method=args.method
my_step=args.my_step
my_lr=args.my_lr
limit_LAD=args.limit_LAD
constant_LAD=args.constant_LAD
limit_res=args.limit_res
constant_res=args.constant_res
limit_power_LAD=1/float(args.limit_power_LAD)
limit_power_res=1/float(args.limit_power_res)
save_file_prefix=args.save_file_prefix
LAD_pattern="m_p_LAD"

# chr_name="chr17"
# method="overlap" # "hierarchical" "overlap56"
# cell_name="6" #"7"
# my_step=18000
# my_lr=1e-4
# strain_plot="C"
# LAD_pattern="m_p_LAD"

prefix = '/home/tianyu/'if os.path.exists('/home/tianyu') else '/home/'
png_dir=prefix+"/ygli/gam_paper_data/single_cell_3D/set_1cell/set_1cell_"+str(cell_name)+"/"+save_file_prefix+strain_plot+"/"

# png_dir=prefix+"/ygli/gam_paper_data/single_cell_3D/set_1cell/set_1cell_6/all_chr_10kb_"+strain_plot+"/"
R=4.3975
resolution=1e+4 # consist with input data

#################################################################################
# def load LAD data function
#################################################################################

# def load_LAD_data(chr_name):
#     chr_name_lad=[]
#     lad_location_head=[]
#     lad_location_tail=[]
#     chr_needed=[chr_name]
#     with open(prefix+"/ygli/gam_paper_data/lad_data/ES_non-allelic_LAD_coordinates.bed",'r') as f:
#         all_lines=f.readlines()
#         for each_line in all_lines:
#             chr_name_this_line=each_line.replace("\n","").split("\t")[0]
#             if chr_name_this_line in chr_needed:
#                 chr_name_lad.append(chr_name_this_line)
#                 lad_location_head.append(int(each_line.replace("\n","").split("\t")[1]))
#                 lad_location_tail.append(int(each_line.replace("\n","").split("\t")[2]))
#     print("LAD num in %s is %s"%(chr_needed,len(lad_location_head)))
#     return chr_name_lad,lad_location_head,lad_location_tail

def load_LAD_data_m_p_LAD(chr_name,which_strain):
    chr_name_lad=[]
    lad_location_head=[]
    lad_location_tail=[]
    chr_needed=[chr_name]
    if which_strain=="C":
        # dir_LAD=prefix+"/ygli/gam_paper_data/lad_data/ES_paternal_LAD_coordinates.bed"
        dir_LAD=prefix+"/ygli/gam_paper_data/lad_data/CAST_paternal.bed"
        # dir_LAD=prefix+"/ygli/gam_paper_data/lad_data/129S1_paternal.tsv"
    elif which_strain=="S":
        # dir_LAD=prefix+"/ygli/gam_paper_data/lad_data/ES_maternal_LAD_coordinates.bed"
        dir_LAD=prefix+"/ygli/gam_paper_data/lad_data/129S1_maternal.bed"
        # dir_LAD=prefix+"/ygli/gam_paper_data/lad_data/CAST_maternal.bed"
    with open(dir_LAD,'r') as f:
        all_lines=f.readlines()
        for each_line in all_lines:
            chr_name_this_line=each_line.replace("\n","").split("\t")[0]
            if chr_name_this_line in chr_needed:
                chr_name_lad.append(chr_name_this_line)
                lad_location_head.append(int(each_line.replace("\n","").split("\t")[1]))
                lad_location_tail.append(int(each_line.replace("\n","").split("\t")[2]))
    lad_location_head_sorted=np.sort(lad_location_head)
    sort_order=sorted(range(len(lad_location_head)), key=lambda k: lad_location_head[k])
    lad_location_tail_sorted=[lad_location_tail[i] for i in sort_order]

    # # del the repeat LAD
    # index_repetition=[]
    # right_pointer=lad_location_tail_sorted[0]
    # for i in range(1,len(lad_location_tail_sorted)):
    #     if lad_location_tail_sorted[i]<=right_pointer:
    #         index_repetition.append(i)
    #     else:
    #         right_pointer=lad_location_tail_sorted[i]

    # lad_location_head_after_del=[i for j, i in enumerate(lad_location_head_sorted) if j not in index_repetition]
    # lad_location_tail_after_del=[i for j, i in enumerate(lad_location_tail_sorted) if j not in index_repetition]

# del repeat LAD and merge overlap LAD
    lad_location_head_after_del=[]
    lad_location_tail_after_del=[]
    left_pointer=lad_location_head_sorted[0]
    right_pointer=lad_location_tail_sorted[0]
    for i in range(1,len(lad_location_tail_sorted)):
        if lad_location_head_sorted[i]<=right_pointer:
            # if lad_location_tail_sorted[i]<=right_pointer: #dont add anything
            if lad_location_tail_sorted[i]>right_pointer:
                right_pointer=lad_location_tail_sorted[i]
        else:
            lad_location_head_after_del.append(left_pointer)
            lad_location_tail_after_del.append(right_pointer)
            left_pointer=lad_location_head_sorted[i]
            right_pointer=lad_location_tail_sorted[i]

    print("LAD num in %s is %s"%(chr_needed,len(lad_location_head_after_del)))
    return chr_name_lad,lad_location_head_after_del,lad_location_tail_after_del

#################################################################################
#1. def load covered_point_y data function
#################################################################################

# # all_chr_name=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
# #         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"]
# # all_chr_name=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
# #         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"]
# all_chr_name="chr1"
# # chr_name='chr18'
# strain="C"

def load_NP_data(chr_name,strain):
    chr_name_loci=[]
    loci_location_head=[]
    loci_location_tail=[]
    covered_point_y=[]
    slice_distance=[]

    dic_loci_location_head_covered_point_y={}
    dic_loci_location_head_loci_index={}

#     with open(prefix+"/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/set_1cell/set_1cell_5/all_chr_10kb_"+strain+"/"+chr_name+"_10kb_"+strain+"_distance_data_1.txt",'r') as f:
    with open(prefix+"/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/set_1cell/set_1cell_"+str(cell_name)+"/all_chr_10kb_"+strain+"/"+chr_name+"_10kb_"+strain+"_distance_data_1.txt",'r') as f:
        all_lines=f.readlines()
        index=0
        for each_line in all_lines:
            covered_point_y.append([R-float(x) for x in each_line.replace("\n","").split("\t")[3:]])
            chr_name_loci.append(each_line.replace("\n","").split("\t")[0])
            loci_location_head.append(int(each_line.replace("\n","").split("\t")[1]))
            loci_location_tail.append(int(each_line.replace("\n","").split("\t")[2]))            
            dic_loci_location_head_covered_point_y[int(each_line.replace("\n","").split("\t")[1])]=(R-float(each_line.replace("\n","").split("\t")[3]))
            dic_loci_location_head_loci_index[int(each_line.replace("\n","").split("\t")[1])]=index
            index+=1

    for slice_distance_i in range(len(covered_point_y[1])):
        slice_distance_each=[(covered_point_y[i][slice_distance_i]) for i in range(len(covered_point_y))]
#         slice_distance_each=[np.abs(covered_point_y[i][slice_distance_i]) for i in range(len(covered_point_y))]
        slice_distance.append(slice_distance_each)

    sample_num=len(covered_point_y)
    sequential_num=len(covered_point_y[0])

    # print(covered_point_y[1])
    print("loci num:",sample_num)

    # slice_distance_0=[covered_point_y[i][0] for i in range(len(covered_point_y))]

    print("sequential_num:",sequential_num)
    print("len(slice_distance)",len(slice_distance))
    return chr_name_loci,loci_location_head,loci_location_tail,covered_point_y,slice_distance,dic_loci_location_head_covered_point_y,dic_loci_location_head_loci_index,sample_num,sequential_num

#################################################################################
#2. load data covered_point_y
#################################################################################
# all_chr_name=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
#         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"]
# all_chr_name=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
#         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"]

strain="C"
# chr_name_loci,loci_location_head,loci_location_tail,covered_point_y,slice_distance,dic_loci_location_head_covered_point_y
chr_name_loci_C,loci_location_head_C,loci_location_tail_C,covered_point_y_C,slice_distance_C,dic_loci_location_head_covered_point_y_C,dic_loci_location_head_loci_index_C,sample_num_C,sequential_num_C=load_NP_data(chr_name,strain)
strain="S"
chr_name_loci_S,loci_location_head_S,loci_location_tail_S,covered_point_y_S,slice_distance_S,dic_loci_location_head_covered_point_y_S,dic_loci_location_head_loci_index_S,sample_num_S,sequential_num_S=load_NP_data(chr_name,strain)
if LAD_pattern==[]:
    chr_name_lad,lad_location_head,lad_location_tail=load_LAD_data(chr_name)
elif LAD_pattern=="m_p_LAD":
    chr_name_lad_C,lad_location_head_C,lad_location_tail_C=load_LAD_data_m_p_LAD(chr_name,"C")
    chr_name_lad_S,lad_location_head_S,lad_location_tail_S=load_LAD_data_m_p_LAD(chr_name,"S")
#################################################################################
#3. build the 10kb loci set, point_location_head, point_h, point_phi, point_r
#
# inclued len_this_chr: captured lens, len_this_chr_array: array lens
# point_location_head_set,point_location_tail_set,point_h_set,point_phi_set,point_r_set
# dic_point_location_head_index: D[point_location_head]=point_index
#################################################################################

def build_basic_point_info(resolution,loci_location_head,dic_loci_location_head_covered_point_y):
    len_this_chr=np.max(loci_location_head)-np.min(loci_location_head)
    len_this_chr_array=int(np.ceil(len_this_chr/resolution))+1

    # print(len_this_chr_array)
    point_location_head_set=np.zeros([len_this_chr_array])
    point_h_set=np.zeros([len_this_chr_array])
    point_phi_set=np.zeros([len_this_chr_array])
    point_r_set=np.zeros([len_this_chr_array])

    point_location_head_set=np.linspace(np.min(loci_location_head),np.max(loci_location_head),len_this_chr_array)
    point_location_tail_set=point_location_head_set+resolution

    for i in range(len_this_chr_array):
        if int(point_location_head_set[i]) in dic_loci_location_head_covered_point_y:
            point_h_set[i]=dic_loci_location_head_covered_point_y[int(point_location_head_set[i])]

    dic_point_location_head_index={}
    i_index=0
    for i_loci_location in point_location_head_set:
        dic_point_location_head_index[i_loci_location]=i_index
        i_index+=1

    # print(point_location_head_set[0:10])
    # print(dic_loci_location_head_covered_point_y.keys()[0])
    # print(dic_loci_location_head_covered_point_y[7000])

    # print(ii)
    print("len_this_chr_captured:",len_this_chr)
    print("len_this_chr_array:",len_this_chr_array)
    return len_this_chr,len_this_chr_array,point_location_head_set,point_location_tail_set,point_h_set,point_phi_set,point_r_set,dic_point_location_head_index

len_this_chr_C,len_this_chr_array_C,point_location_head_set_C,point_location_tail_set_C,point_h_set_C,point_phi_set_C,point_r_set_C,dic_point_location_head_index_C=build_basic_point_info(resolution,loci_location_head_C,dic_loci_location_head_covered_point_y_C)
len_this_chr_S,len_this_chr_array_S,point_location_head_set_S,point_location_tail_set_S,point_h_set_S,point_phi_set_S,point_r_set_S,dic_point_location_head_index_S=build_basic_point_info(resolution,loci_location_head_S,dic_loci_location_head_covered_point_y_S)


#################################################################################
#4. calculate linear distance for each loci
#################################################################################
def calculate_nearest_distance(loci_location_head,loci_location_tail,lad_location_head,lad_location_tail):
    nearest_lad_index=[]
    nearest_head_lad_distance=[] # if 0, loci in this lad
    nearest_tail_lad_distance=[]
    dic_loci_location_head_nearest_lad_distance={}

    begin_in_loci=0 #the index of last loci in the first LAD
    last_head_lad_index=0
#     for each_loci_num in range(10):
    for each_loci_num in range(len(loci_location_head)):
        flag_get_tail_lad=0
        head_lad_index=0
        for lad_index in range(last_head_lad_index,len(lad_location_head)):
            #to find the first lad that loci_tail<lad_head, regard as tail_lad
            if loci_location_tail[each_loci_num]-lad_location_head[lad_index]<=0:
                flag_get_tail_lad=1
                break
        if lad_index==0: # if the tail_lad is the first lad
            begin_in_loci=each_loci_num
            nearest_head_lad_distance_i=sys.maxint
            nearest_tail_lad_distance_i=lad_location_head[lad_index]-loci_location_tail[each_loci_num]
        elif flag_get_tail_lad==0: #if the last LAD is still not the tail_lad
    #         print("loci in lad's right index:",each_loci_num)
            head_lad_index=lad_index
            last_head_lad_index=head_lad_index
            # if the loci head locate ahead of lad tail
            if loci_location_head[each_loci_num]-lad_location_tail[head_lad_index]<=0:
                nearest_head_lad_distance_i=nearest_tail_lad_distance_i=0
            # elif the loci head locate behind lad tail
            else:
                nearest_tail_lad_distance_i=sys.maxint
                nearest_head_lad_distance_i=loci_location_head[each_loci_num]-lad_location_tail[head_lad_index]
        else: # find out the tail_lad and tail_lad>=1
            head_lad_index=lad_index-1
            last_head_lad_index=head_lad_index
            nearest_tail_lad_distance_i=lad_location_head[lad_index]-loci_location_tail[each_loci_num]
            if lad_location_tail[head_lad_index]-loci_location_head[each_loci_num]<=0:
                nearest_head_lad_distance_i=loci_location_head[each_loci_num]-lad_location_tail[head_lad_index]
            else:
                nearest_head_lad_distance_i=0
            # if loci_location_tail[each_loci_num]-lad_location_tail[head_lad_index]<=0:
            #     if loci_location_head[each_loci_num]-lad_location_head[head_lad_index]>=0:
            #         nearest_head_lad_distance_i=nearest_tail_lad_distance_i=0
            #     else:
            #         # last_head_lad_index=head_lad_index
            #         nearest_head_lad_distance_i=0
            #         nearest_tail_lad_distance_i=lad_location_head[lad_index]-loci_location_tail[each_loci_num]
            # else:
            #     last_head_lad_index=head_lad_index
            #     nearest_head_lad_distance_i=loci_location_head[each_loci_num]-lad_location_tail[head_lad_index]
            #     nearest_tail_lad_distance_i=lad_location_head[lad_index]-loci_location_tail[each_loci_num]
            #     # if nearest_head_lad_distance_i<0:
                #     print("!!nearest_head_lad_distance_i<0, head_lad_index is %s, each_loci_num is %s"%(head_lad_index,each_loci_num))
                #     print("loci_location_head[each_loci_num]:%s, lad_location_tail[head_lad_index]:%s, nearest_tail_lad_distance_i:%s"%(loci_location_head[each_loci_num],lad_location_tail[head_lad_index],nearest_tail_lad_distance_i))
        nearest_lad_index.append(head_lad_index)
        nearest_head_lad_distance.append(nearest_head_lad_distance_i)
        nearest_tail_lad_distance.append(nearest_tail_lad_distance_i)
        
        dic_loci_location_head_nearest_lad_distance[loci_location_head[each_loci_num]]=np.min([nearest_head_lad_distance_i,nearest_tail_lad_distance_i],0)
#         print(nearest_head_lad_distance_i,nearest_tail_lad_distance_i)
#         print(np.min([nearest_head_lad_distance_i,nearest_tail_lad_distance_i],0))
        
    return nearest_lad_index,nearest_head_lad_distance,nearest_tail_lad_distance,dic_loci_location_head_nearest_lad_distance

# for each loci, search the LAD from last_head_index(=last time find nearest lad) to end
if LAD_pattern==[]:
    nearest_lad_index_C,nearest_head_lad_distance_C,nearest_tail_lad_distance_C,dic_loci_location_head_nearest_lad_distance_C=calculate_nearest_distance(loci_location_head_C,loci_location_tail_C,lad_location_head,lad_location_tail)
    nearest_lad_index_S,nearest_head_lad_distance_S,nearest_tail_lad_distance_S,dic_loci_location_head_nearest_lad_distance_S=calculate_nearest_distance(loci_location_head_S,loci_location_tail_S,lad_location_head,lad_location_tail)
elif LAD_pattern=="m_p_LAD":
    nearest_lad_index_C,nearest_head_lad_distance_C,nearest_tail_lad_distance_C,dic_loci_location_head_nearest_lad_distance_C=calculate_nearest_distance(loci_location_head_C,loci_location_tail_C,lad_location_head_C,lad_location_tail_C)
    nearest_lad_index_S,nearest_head_lad_distance_S,nearest_tail_lad_distance_S,dic_loci_location_head_nearest_lad_distance_S=calculate_nearest_distance(loci_location_head_S,loci_location_tail_S,lad_location_head_S,lad_location_tail_S)


#################################################################################
#5. calculate the point_r_set for all points in 10kb
#################################################################################
if LAD_pattern==[]:
    nearest_lad_index_point_C,nearest_head_lad_distance_point_C,nearest_tail_lad_distance_point_C,dic_loci_location_head_nearest_lad_distance_point_C=calculate_nearest_distance(point_location_head_set_C,point_location_tail_set_C,lad_location_head,lad_location_tail)
    nearest_lad_index_point_S,nearest_head_lad_distance_point_S,nearest_tail_lad_distance_point_S,dic_loci_location_head_nearest_lad_distance_point_S=calculate_nearest_distance(point_location_head_set_S,point_location_tail_set_S,lad_location_head,lad_location_tail)
elif LAD_pattern=="m_p_LAD":
    nearest_lad_index_point_C,nearest_head_lad_distance_point_C,nearest_tail_lad_distance_point_C,dic_loci_location_head_nearest_lad_distance_point_C=calculate_nearest_distance(point_location_head_set_C,point_location_tail_set_C,lad_location_head_C,lad_location_tail_C)
    nearest_lad_index_point_S,nearest_head_lad_distance_point_S,nearest_tail_lad_distance_point_S,dic_loci_location_head_nearest_lad_distance_point_S=calculate_nearest_distance(point_location_head_set_S,point_location_tail_set_S,lad_location_head_S,lad_location_tail_S)

point_r_set_kb_linear_C=np.min([nearest_head_lad_distance_point_C,nearest_tail_lad_distance_point_C],0)
point_r_set_kb_linear_S=np.min([nearest_head_lad_distance_point_S,nearest_tail_lad_distance_point_S],0)

#################################################################################
#6. calculate the point_r_set_np for all loci in 10kb
# include dic_slice_loci_index,dic_loci_in_lad,
# dic_slice_loci_index:[h]=loci_index; 
# dic_loci_in_lad:[h]=loci_index in LAD
# point_r_set_kb_np
#################################################################################

# build the dictionary of loci in each NP, dic_slice_loci_index[h]=[each_np]
def cal_point_r_set_kb_np(sample_num,len_this_chr_array,covered_point_y,dic_loci_location_head_nearest_lad_distance_point,loci_location_head,loci_location_tail,dic_point_location_head_index):

    dic_slice_loci_index={}
    for i in range(sample_num):
        if covered_point_y[i][0] in dic_slice_loci_index:
            dic_slice_loci_index[covered_point_y[i][0]].append(i)
        else:
            dic_slice_loci_index[covered_point_y[i][0]]=[i]

    # build the dictionary of loci in LAD in each NP, dic_slice_loci_index[h]=[each_np]
    dic_loci_in_lad={}
    for i in range(sample_num):
        if dic_loci_location_head_nearest_lad_distance_point[loci_location_head[i]]==0:
            if covered_point_y[i][0] in dic_loci_in_lad:
                dic_loci_in_lad[covered_point_y[i][0]].append(i)
            else:
                dic_loci_in_lad[covered_point_y[i][0]]=[i]

    point_r_set_kb_np=[sys.maxint for i in range(len_this_chr_array)]
    # calculate point_r_set_np for all loci
    for each_h in dic_slice_loci_index.keys():
        point_location_head_set_np=[ loci_location_head[each_loci] for each_loci in dic_slice_loci_index[each_h] ]
        point_location_tail_set_np=[ loci_location_tail[each_loci] for each_loci in dic_slice_loci_index[each_h] ]
        
        if each_h in dic_loci_in_lad:
            lad_location_head_np=[ loci_location_head[each_loci] for each_loci in dic_loci_in_lad[each_h] ]
            lad_location_tail_np=[ loci_location_tail[each_loci] for each_loci in dic_loci_in_lad[each_h] ]
            nearest_lad_index_np,nearest_head_lad_distance_np,nearest_tail_lad_distance_np,dic_loci_location_nearest_lad_distance_np=calculate_nearest_distance(point_location_head_set_np,point_location_tail_set_np,lad_location_head_np,lad_location_tail_np)
        else:
            lad_location_head_np=[]
            lad_location_tail_np=[]
            nearest_lad_index_np,nearest_head_lad_distance_np,nearest_tail_lad_distance_np,dic_loci_location_nearest_lad_distance_np=[],[],[],{}
        
        for each_loci_head in dic_loci_location_nearest_lad_distance_np.keys():
            point_r_set_kb_np[dic_point_location_head_index[each_loci_head]]=dic_loci_location_nearest_lad_distance_np[each_loci_head]
    
    return dic_slice_loci_index,dic_loci_in_lad,point_r_set_kb_np


dic_slice_loci_index_C,dic_loci_in_lad_C,point_r_set_kb_np_C=cal_point_r_set_kb_np(sample_num_C,len_this_chr_array_C,covered_point_y_C,dic_loci_location_head_nearest_lad_distance_point_C,loci_location_head_C,loci_location_tail_C,dic_point_location_head_index_C)
dic_slice_loci_index_S,dic_loci_in_lad_S,point_r_set_kb_np_S=cal_point_r_set_kb_np(sample_num_S,len_this_chr_array_S,covered_point_y_S,dic_loci_location_head_nearest_lad_distance_point_S,loci_location_head_S,loci_location_tail_S,dic_point_location_head_index_S)

# print("======C======")
# print("np.max(point_r_set_kb_linear_C):",np.max(point_r_set_kb_linear_C))
# print("np.max(point_r_set_kb_np_C):",np.max(point_r_set_kb_np_C))

# print("point_r_set_kb_linear_C:",point_r_set_kb_linear_C[0])
# print("point_r_set_kb_linear_C:",point_r_set_kb_linear_C[0:20])
# print("point_r_set_kb_np:",point_r_set_kb_np[0:20])

# print("point_r_set_kb_np_C:",point_r_set_kb_np_C[0])
# print("point_r_set_kb_np_C:",point_r_set_kb_np_C[0:20])

# print("======S======")
# print("np.max(point_r_set_kb_linear_S):",np.max(point_r_set_kb_linear_S))
# print("np.max(point_r_set_kb_np_S):",np.max(point_r_set_kb_np_S))

# print("point_r_set_kb_linear_S:",point_r_set_kb_linear_S[0])
# print("point_r_set_kb_linear_S:",point_r_set_kb_linear_S[0:20])
# # print("point_r_set_kb_np:",point_r_set_kb_np[0:20])

# print("point_r_set_kb_np_S:",point_r_set_kb_np_S[0])
# print("point_r_set_kb_np_S:",point_r_set_kb_np_S[0:20])


#################################################################################
#calculate point_r_set_kb_linear -> point_r_set_linear and point_r_set_kb_np -> point_r_set_np 
# and compare point_r_set_linear and point_r_set_np to get max
#################################################################################
# def the functions to calculate distances

def get_spatial_distance(x):
    intercept=0.45
    para=0.00459
    inflexion=intercept+para*(5*resolution)**(1/3)
    if x<5*resolution:
        y=x/(5*resolution)*inflexion
    else:
        y=intercept+para*(x)**(1/3)
    return y

def get_maxmin_r_linear(R,d,kb_lad):
    l_lad=get_spatial_distance(kb_lad)
    d=np.abs(d)
#     min_dis_in=R**2-2*d*l_lad-l_lad**2
    if d>=R-l_lad:
        min_dis=0
    else:
        min_dis=((R-l_lad)**2-d**2)**(0.5)
    return min_dis

def get_maxmin_r_np(R,d,kb_lad):
    l_lad=get_spatial_distance(kb_lad)
    if R-l_lad<=0:
#         print("R-l_lad<=0")
        min_dis=0
    else:
        min_dis=R-l_lad
    return min_dis

# point_r_set_kb_linear, point_r_set_kb_np
# although only for captured loci, this distance changed, calculate all for convenient
def calculate_point_r_lower_bound(R,sample_num,point_h_set,point_r_set_kb_linear,point_r_set_kb_np):
    point_r_lower_bound=[]
    # for i_point in range(20):
    for i_point in range(sample_num):
        dis_lower_bound=0
        if point_r_set_kb_linear[i_point]!=point_r_set_kb_np[i_point]:
            d=point_h_set[i_point]
            dis_linear=get_maxmin_r_linear(R,d,point_r_set_kb_linear[i_point])
            dis_np=get_maxmin_r_np(R,d,point_r_set_kb_np[i_point])
    #         dis_lower_bound=0
            del dis_lower_bound
    #         print("dis_linear>dis_np?",dis_linear>dis_np)
            dis_lower_bound=max([dis_linear,dis_np])
        else:
            dis_lower_bound=0
    #     print("dis_linear:%s,dis_np:%s,dis_lower_bound:%s"%(dis_linear,dis_np,dis_lower_bound))
        point_r_lower_bound.append(dis_lower_bound)
    return point_r_lower_bound

point_r_lower_bound_C=calculate_point_r_lower_bound(R,sample_num_C,point_h_set_C,point_r_set_kb_linear_C,point_r_set_kb_np_C)
point_r_lower_bound_S=calculate_point_r_lower_bound(R,sample_num_S,point_h_set_S,point_r_set_kb_linear_S,point_r_set_kb_np_S)

#################################################################################
# load SNP official data to test the SNP found
#################################################################################
# def load_SNP_data(chr_name,resolution):
#     dir_SNP=prefix+"/ygli/gam_paper_data/f123_SNP/GSE48592_castx129_variants.vcf.txt"
#     chr_needed=[chr_name]
#     SNP_location_complete=[]
#     with open(dir_SNP,'r') as f:
#         all_lines=f.readlines()
#         for each_line in all_lines:
#             chr_name_this_l,pos_this_l,_,_,_,_,filter_thisl=each_line.replace("\n","").split("\t")[0:7]
#             if chr_name_this_l in chr_needed and filter_thisl=="PASS":
#                 SNP_location_complete.append(int(int(int(pos_this_l)/resolution)*resolution))
#     return SNP_location_complete

# SNP_location_complete=load_SNP_data(chr_name,resolution)

## 'feat_name', 'source_int', 'mapped_int', 'source_id', 'mapped_id', 'source_length', 'mapped_length', 'source_start', 'source_stop', 'source_strand', 'source_sub_start', 'source_sub_stop', 'mapped_start', 'mapped_stop', 'mapped_strand', 'coverage', 'recip', 'asm_unit', '
def load_SNP_C(chr_name,resolution):
    dir_SNP=prefix+"/ygli/gam_paper_data/f123_SNP/CAST_snps_indels.tsv"
    chr_needed=[chr_name]
    SNP_location_complete=[]
    i=0
    with open(dir_SNP,'r') as f:
        # _=f.readline().split("\t")
        all_lines=f.readlines()
        for each_line in all_lines:
            i+=1
            if each_line[0]!="#":
                this_line=each_line.replace('\n','').split()
                if this_line[2]!='NULL':
                    chr_name_this_l=this_line[3]
                    mapped_start,mapped_stop=this_line[12:14]
                    SNP_begin=int(int(int(mapped_start)/resolution)*resolution)
                    SNP_end=int(int(int(mapped_stop)/resolution)*resolution)
                    if chr_name_this_l in chr_needed:
                        for SNP_location in range(SNP_begin,SNP_end+1,int(resolution)):
                            SNP_location_complete.append(SNP_location)
    return SNP_location_complete
    
## 'feat_name', 'source_int', 'mapped_int', 'source_id', 'mapped_id', 'source_length', 'mapped_length', 'source_start', 'source_stop', 'source_strand', 'source_sub_start', 'source_sub_stop', 'mapped_start', 'mapped_stop', 'mapped_strand', 'coverage', 'recip', 'asm_unit', '
def load_SNP_S(chr_name,resolution):
    dir_SNP=prefix+"/ygli/gam_paper_data/f123_SNP/129S1_snps_indels.tsv"
    chr_needed=[chr_name]
    SNP_location_complete=[]
    i=0
    with open(dir_SNP,'r') as f:
        # _=f.readline().split("\t")
        all_lines=f.readlines()
        for each_line in all_lines:
            i+=1
            if each_line[0]!="#":
                this_line=each_line.replace('\n','').split()
                if this_line[2]!='NULL':
                    chr_name_this_l=this_line[3]
                    mapped_start,mapped_stop=this_line[12:14]
                    SNP_begin=int(int(int(mapped_start)/resolution)*resolution)
                    SNP_end=int(int(int(mapped_stop)/resolution)*resolution)
                    if chr_name_this_l in chr_needed:
                        for SNP_location in range(SNP_begin,SNP_end+1,int(resolution)):
                            SNP_location_complete.append(SNP_location)
    return SNP_location_complete

SNP_location_complete_C=load_SNP_C(chr_name,resolution)
SNP_location_complete_S=load_SNP_S(chr_name,resolution)

#################################################################################
# use SNP official data to build the loci captured
# for dic_loci_location_head_loci_index_C,dic_loci_location_head_loci_index_S
# contain C snp -> add to dic_SNP_authented_location_loci_index_C 
# contain S snp -> add to dic_SNP_authented_location_loci_index_S
# contain non snp in C -> add to dic_nonSNP_authented_location_loci_index_C
# contain non snp in S -> add to dic_nonSNP_authented_location_loci_index_S
#################################################################################

def build_SNP_loci(dic_loci_location_head_loci_index,SNP_location_complete):
    dic_SNP_authented_location_loci_index={}
    dic_nonSNP_authented_location_loci_index={}
    i_out=0
    for loci_location_i in dic_loci_location_head_loci_index:
        if loci_location_i in SNP_location_complete:
            dic_SNP_authented_location_loci_index[loci_location_i]=dic_loci_location_head_loci_index[loci_location_i]
        else:
            i_out+=1
            dic_nonSNP_authented_location_loci_index[loci_location_i]=dic_loci_location_head_loci_index[loci_location_i]
    print("all loci number is %s,nonSNP number is %s,percent is %s"%(len(dic_loci_location_head_loci_index.keys()),i_out,i_out/len(dic_loci_location_head_loci_index.keys())))
    return dic_SNP_authented_location_loci_index,dic_nonSNP_authented_location_loci_index

dic_SNP_authented_location_loci_index_C,dic_nonSNP_authented_location_loci_index_C=build_SNP_loci(dic_loci_location_head_loci_index_C,SNP_location_complete_C)
dic_SNP_authented_location_loci_index_S,dic_nonSNP_authented_location_loci_index_S=build_SNP_loci(dic_loci_location_head_loci_index_S,SNP_location_complete_S)


#################################################################################
# plot 2 group SNP linear location, h, r_lower_bound
#################################################################################

def get_SNP_info(dic_SNP_authented_location_loci_index,dic_loci_location_head_covered_point_y,point_r_lower_bound):
    SNP_linear_location=[]
    SNP_h=[]
    SNP_r_lower_bound=[]
    SNP_index=[]
    for SNP_location_i in dic_SNP_authented_location_loci_index:
        index=dic_SNP_authented_location_loci_index[SNP_location_i]
        SNP_linear_location.append(SNP_location_i)
        SNP_h.append(dic_loci_location_head_covered_point_y[SNP_location_i])
        SNP_r_lower_bound.append(point_r_lower_bound[index])
        SNP_index.append(index)
    return SNP_linear_location,SNP_h,SNP_r_lower_bound,SNP_index

SNP_linear_location_C,SNP_h_C,SNP_r_lower_bound_C,SNP_index_C=get_SNP_info(dic_SNP_authented_location_loci_index_C,dic_loci_location_head_covered_point_y_C,point_r_lower_bound_C)
SNP_linear_location_S,SNP_h_S,SNP_r_lower_bound_S,SNP_index_S=get_SNP_info(dic_SNP_authented_location_loci_index_S,dic_loci_location_head_covered_point_y_S,point_r_lower_bound_S)
SNP_label_C=["SNP_CAST"]*len(SNP_linear_location_C)
SNP_label_S=["SNP_129S1"]*len(SNP_linear_location_S)

# SNP_info=pd.DataFrame({"SNP_linear_location":np.array(SNP_linear_location_C+SNP_linear_location_S),
#                      "SNP_h":np.array(SNP_h_C+SNP_h_S),
#                      "SNP_r_lower_bound":np.array(SNP_r_lower_bound_C+SNP_r_lower_bound_S),
#                      "SNP label":np.array(SNP_label_C+SNP_label_S),
#                     })

SNP_info_C=pd.DataFrame({"SNP_linear_location":np.array(SNP_linear_location_C),
                     "SNP_h":np.array(SNP_h_C),
                     "SNP_r_lower_bound":np.array(SNP_r_lower_bound_C),
                     "SNP label":np.array(SNP_label_C),
                    })

SNP_info_S=pd.DataFrame({"SNP_linear_location":np.array(SNP_linear_location_S),
                     "SNP_h":np.array(SNP_h_S),
                     "SNP_r_lower_bound":np.array(SNP_r_lower_bound_S),
                     "SNP label":np.array(SNP_label_S),
                    })
# plot the 
# fig = plt.figure(figsize=(16, 16))
print("plot SNP_h with different SNP_linear_location ...")
fig = plt.figure(figsize=(16, 12))
ax=fig.add_subplot(211)
ax1=sns.scatterplot(SNP_info_C["SNP_linear_location"],SNP_info_C["SNP_h"],hue=SNP_info_C["SNP label"])
ax1.set_xlabel('SNP_linear_location')
ax1.set_ylabel('SNP_h')
plt.title("SNP_h with different SNP_linear_location for CAST")
# plt.xlim([0,loci_location_range])
# plt.ylim([0,np.max(mydata["slice_distance"])])
# plt.xlim([0.25e+8,0.5e+8])

ax=fig.add_subplot(212)
ax2=sns.scatterplot(SNP_info_S["SNP_linear_location"],SNP_info_S["SNP_h"],hue=SNP_info_S["SNP label"])
ax2.set_xlabel('SNP_linear_location')
ax2.set_ylabel('SNP_h')
plt.title("SNP_h with different SNP_linear_location for 129S1")

plt.savefig(png_dir+chr_name+"_10kb_SNP_h_SNP_linear_location.eps",dpi=600,format='eps')
plt.close()

# plot the 3D fig for SNP 
print("plot SNP_h with different SNP_linear_location_3D ...")
fig = plt.figure()
ax = plt.axes(projection='3d')
# ax.scatter3D(xdata, ydata, zdata, cmap='Greens');
ax.scatter3D(SNP_linear_location_C,SNP_h_C,SNP_r_lower_bound_C, c='r', label='SNP CAST')
ax.scatter3D(SNP_linear_location_S,SNP_h_S,SNP_r_lower_bound_S, c='g', label='SNP 129S1')
ax.legend(loc='best')
ax.set_xlabel('linear_location', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('h', fontdict={'size': 15, 'color': 'red'})
ax.set_zlabel('r_lower_bound', fontdict={'size': 15, 'color': 'red'})
# plt.show()
plt.savefig(png_dir+chr_name+"_10kb_SNP_h_SNP_linear_location_3D.eps",dpi=600,format='eps')
plt.close()

#################################################################################
# use fisher aka LDA to classify the nonSNP loci
#################################################################################

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.utils import shuffle
X = np.array([[SNP_linear_location_C[i],SNP_h_C[i]] for i in range(len(SNP_h_C))]+[[SNP_linear_location_S[i],SNP_h_S[i]] for i in range(len(SNP_h_S))])
y = np.array([0]*len(SNP_h_C)+[1]*len(SNP_h_S))
cc = list(zip(X, y)) 
X[:], y[:] = zip(*cc) 
X,y=shuffle(X,y,random_state=0)

clf = LinearDiscriminantAnalysis()
clf.fit(X, y)

# get the info for nonSNP
nonSNP_linear_location_C,nonSNP_h_C,nonSNP_r_lower_bound_C,nonSNP_index_C=get_SNP_info(dic_nonSNP_authented_location_loci_index_C,dic_loci_location_head_covered_point_y_C,point_r_lower_bound_C)
nonSNP_linear_location_S,nonSNP_h_S,nonSNP_r_lower_bound_S,nonSNP_index_S=get_SNP_info(dic_nonSNP_authented_location_loci_index_S,dic_loci_location_head_covered_point_y_S,point_r_lower_bound_S)
x_new_C=[[nonSNP_linear_location_C[i],nonSNP_h_C[i]] for i in range(len(nonSNP_h_C))]
x_new_S=[[nonSNP_linear_location_S[i],nonSNP_h_S[i]] for i in range(len(nonSNP_h_S))]

if len(x_new_C)>0:
    nonSNP_group_C=clf.predict(x_new_C)
else:
    nonSNP_group_C=[]

if len(x_new_S)>0:
    nonSNP_group_S=clf.predict(x_new_S)
else:
    nonSNP_group_S=[]

print(nonSNP_group_C[0:20])
print(nonSNP_group_S[0:20])
print("C in C percent %s"%(1-np.sum(nonSNP_group_C)/len(nonSNP_group_C)))
print("S in S percent %s"%(np.sum(nonSNP_group_S)/len(nonSNP_group_S)))

print("C")
print("loci captured total:",len(dic_loci_location_head_loci_index_C.keys()))
print("loci SNP total:",len(SNP_linear_location_C))
print("loci not SNP total:",len(nonSNP_linear_location_C))
# print("loci not SNP total:",len(dic_nonSNP_authented_location_loci_index_C.keys()))
print("loci in noSNP is S:",np.sum(nonSNP_group_C))

print("S")
print("loci captured total:",len(dic_loci_location_head_loci_index_S.keys()))
print("loci SNP total:",len(SNP_linear_location_S))
print("loci not SNP total:",len(nonSNP_linear_location_S))
# print("loci not SNP total:",len(dic_nonSNP_authented_location_loci_index_S.keys()))
print("loci in noSNP is S:",np.sum(nonSNP_group_S))

#################################################################################
# build the dictionary for nonSNP_loci_location to group, 0:C, 1:S
#################################################################################

def build_classify_dic_nonSNP(nonSNP_group,nonSNP_linear_location,nonSNP_h,nonSNP_r_lower_bound,nonSNP_index):
    dic_nonSNP_location_group={}
    dic_nonSNP_location_h={}
    dic_nonSNP_location_r={}
    dic_nonSNP_location_index={}
    for i in range(len(nonSNP_linear_location)):
        dic_nonSNP_location_group[nonSNP_linear_location[i]]=nonSNP_group[i]
        dic_nonSNP_location_h[nonSNP_linear_location[i]]=nonSNP_h[i]
        dic_nonSNP_location_r[nonSNP_linear_location[i]]=nonSNP_r_lower_bound[i]
        dic_nonSNP_location_index[nonSNP_linear_location[i]]=nonSNP_index[i]
    return dic_nonSNP_location_group,dic_nonSNP_location_h,dic_nonSNP_location_r,dic_nonSNP_location_index

dic_nonSNP_location_group_from_C,dic_nonSNP_location_h_from_C,dic_nonSNP_location_r_from_C,dic_nonSNP_location_index_C=build_classify_dic_nonSNP(nonSNP_group_C,nonSNP_linear_location_C,nonSNP_h_C,nonSNP_r_lower_bound_C,nonSNP_index_C)
dic_nonSNP_location_group_from_S,dic_nonSNP_location_h_from_S,dic_nonSNP_location_r_from_S,dic_nonSNP_location_index_S=build_classify_dic_nonSNP(nonSNP_group_S,nonSNP_linear_location_S,nonSNP_h_S,nonSNP_r_lower_bound_S,nonSNP_index_S)

# build the dictionary for total nonSNP_loci_location to h and r
# actually, the data for group_C group_S should be the same
# dic_nonSNP_loci_location_group = dic_nonSNP_location_group_from_C.copy()
# dic_nonSNP_loci_location_group.update(dic_nonSNP_location_group_from_S)
# dic_nonSNP_loci_location_h = dic_nonSNP_location_h_from_C.copy()
# dic_nonSNP_loci_location_h.update(dic_nonSNP_location_h_from_S)
# dic_nonSNP_loci_location_r=dic_nonSNP_location_r_from_C.copy()
# dic_nonSNP_loci_location_r.update(dic_nonSNP_location_r_from_S)

#################################################################################
# build the info for two group
#################################################################################

def get_SNP_info_dic(dic_SNP_authented_location_loci_index,dic_loci_location_head_covered_point_y,point_r_lower_bound):
    dic_SNP_linear_location_loci_index={}
    dic_SNP_linear_location_h={}
    dic_SNP_linear_location_r_lower_bound={}
    for SNP_location_i in dic_SNP_authented_location_loci_index:
        index=dic_SNP_authented_location_loci_index[SNP_location_i]
        dic_SNP_linear_location_loci_index[SNP_location_i]=index
        dic_SNP_linear_location_h[SNP_location_i]=dic_loci_location_head_covered_point_y[SNP_location_i]
        dic_SNP_linear_location_r_lower_bound[SNP_location_i]=point_r_lower_bound[index]
    return dic_SNP_linear_location_loci_index,dic_SNP_linear_location_h,dic_SNP_linear_location_r_lower_bound

dic_SNP_linear_location_loci_index_C,dic_SNP_linear_location_h_C,dic_SNP_linear_location_r_lower_bound_C=get_SNP_info_dic(dic_SNP_authented_location_loci_index_C,dic_loci_location_head_covered_point_y_C,point_r_lower_bound_C)
dic_SNP_linear_location_loci_index_S,dic_SNP_linear_location_h_S,dic_SNP_linear_location_r_lower_bound_S=get_SNP_info_dic(dic_SNP_authented_location_loci_index_S,dic_loci_location_head_covered_point_y_S,point_r_lower_bound_S)

# build the two group loci, include snp nonsnp. group flag: 0:C, 1:S
def build_loci_info_for_group(group_flag,dic_SNP_linear_location_loci_index,dic_SNP_linear_location_h,dic_SNP_linear_location_r_lower_bound,dic_nonSNP_location_index,dic_nonSNP_loci_location_group,dic_nonSNP_loci_location_h,dic_nonSNP_loci_location_r):
    dic_group_loci_location_wSNP={}
    dic_group_loci_location_loci_index=dic_SNP_linear_location_loci_index.copy()
    dic_group_loci_location_h=dic_SNP_linear_location_h.copy()
    dic_group_loci_location_r=dic_SNP_linear_location_r_lower_bound.copy()
    # for SNP
    for loci_location_i in dic_SNP_linear_location_loci_index:
        dic_group_loci_location_wSNP[loci_location_i]="Y"
    # for nonSNP
    for loci_location_i in dic_nonSNP_loci_location_group:
        if dic_nonSNP_loci_location_group[loci_location_i]==group_flag:
            dic_group_loci_location_wSNP[loci_location_i]="N"
            dic_group_loci_location_loci_index[loci_location_i]=dic_nonSNP_location_index[loci_location_i]
            dic_group_loci_location_h[loci_location_i]=dic_nonSNP_loci_location_h[loci_location_i]
            dic_group_loci_location_r[loci_location_i]=dic_nonSNP_loci_location_r[loci_location_i]
    return dic_group_loci_location_wSNP,dic_group_loci_location_loci_index,dic_group_loci_location_h,dic_group_loci_location_r

dic_group_loci_location_wSNP_C,dic_group_loci_location_loci_index_C,dic_group_loci_location_h_C,dic_group_loci_location_r_C=build_loci_info_for_group(0,dic_SNP_linear_location_loci_index_C,dic_SNP_linear_location_h_C,dic_SNP_linear_location_r_lower_bound_C,dic_nonSNP_location_index_C,dic_nonSNP_location_group_from_C,dic_nonSNP_location_h_from_C,dic_nonSNP_location_r_from_C)
dic_group_loci_location_wSNP_S,dic_group_loci_location_loci_index_S,dic_group_loci_location_h_S,dic_group_loci_location_r_S=build_loci_info_for_group(1,dic_SNP_linear_location_loci_index_S,dic_SNP_linear_location_h_S,dic_SNP_linear_location_r_lower_bound_S,dic_nonSNP_location_index_S,dic_nonSNP_location_group_from_S,dic_nonSNP_location_h_from_S,dic_nonSNP_location_r_from_S)

print("in group C, the loci num is:%s"%(len(dic_group_loci_location_wSNP_C.keys())))
print("in group S, the loci num is:%s"%(len(dic_group_loci_location_wSNP_S.keys())))

if strain_plot=="C":
    var_list_location=[]
    for loci_location_i in dic_group_loci_location_h_C:
        var_list_location.append(loci_location_i)
    sorted_var_location=np.sort(var_list_location)
    sorted_h=[]
    sorted_nearest_lad_distance_all=[]
    # sorted_r=[]
    for loci_location_i in sorted_var_location:
        sorted_h.append(dic_group_loci_location_h_C[loci_location_i])
        # sorted_r.append(dic_group_loci_location_r_C[loci_location_i])
        sorted_nearest_lad_distance_all.append(dic_loci_location_head_nearest_lad_distance_point_C[loci_location_i])
    d=np.where(np.array(sorted_nearest_lad_distance_all)==0)

    print("total loci num:%s,cover_percent:%s, LAD loci num:%s, percent:%s"%(len(sorted_nearest_lad_distance_all),len(sorted_nearest_lad_distance_all)/len_this_chr_array_C,len(d[0]),len(d[0])/len(sorted_nearest_lad_distance_all)))
else:
    var_list_location=[]
    for loci_location_i in dic_group_loci_location_h_S:
        var_list_location.append(loci_location_i)
    sorted_var_location=np.sort(var_list_location)
    sorted_h=[]
    sorted_nearest_lad_distance_all=[]
    # sorted_r=[]
    for loci_location_i in sorted_var_location:
        sorted_h.append(dic_group_loci_location_h_S[loci_location_i])
        # sorted_r.append(dic_group_loci_location_r_C[loci_location_i])
        sorted_nearest_lad_distance_all.append(dic_loci_location_head_nearest_lad_distance_point_S[loci_location_i])
    d=np.where(np.array(sorted_nearest_lad_distance_all)==0)
    print("total loci num:%s,cover_percent:%s, LAD loci num:%s, percent:%s"%(len(sorted_nearest_lad_distance_all),len(sorted_nearest_lad_distance_all)/len_this_chr_array_S,len(d[0]),len(d[0])/len(sorted_nearest_lad_distance_all)))

def get_spatial_distance_list(x):
    intercept=0.45
    para=0.00459
    inflexion=intercept+para*(5*resolution)**(1/3)
    y=np.zeros_like(x)
    for x_i in range(len(x)):
        if x[x_i]<5*resolution:
            y[x_i]=x[x_i]/(5*resolution)*inflexion
        else:
            y[x_i]=intercept+para*(x[x_i])**(1/3)
    return y

def get_spatial_distance_array(x):
    intercept=0.45
    para=0.00459
    inflexion=intercept+para*(5*resolution)**(1/3)
    y=np.zeros_like(x)
    for x_i in range(len(x)):
        for x_j in range(len(x[0])):
            if x[x_i][x_j]<5*resolution:
                y[x_i][x_j]=x[x_i][x_j]/(5*resolution)*inflexion
            else:
                y[x_i][x_j]=intercept+para*(x[x_i][x_j])**(1/3)
    return y

# #################################################################################
# begin function: with r as parameter, with limitation of LAD
# USING loss funciton: calculate the phi for each loci
# # M(r,phi,h)  (r*cos\phi, r*sin\phi, h)   M1M2=np.linalg.norm(M1-M2)
# #################################################################################
import tensorflow as tf
from scipy.spatial.distance import pdist
import math

# use_point_num=200
# use_begin=200

def run_phi_with_LAD(steps,lr,use_point_num,use_begin,dic_group_loci_location_h_cal,dic_loci_location_head_nearest_lad_distance_point_cal):

    var_list_location=[]
    for loci_location_i in dic_group_loci_location_h_cal:
        var_list_location.append(loci_location_i)

    sorted_var_location=np.sort(var_list_location)
    sorted_h=[]
    sorted_nearest_lad_distance_all=[]
    # sorted_r=[]
    for loci_location_i in sorted_var_location:
        sorted_h.append(dic_group_loci_location_h_cal[loci_location_i])
        # sorted_r.append(dic_group_loci_location_r_C[loci_location_i])
        sorted_nearest_lad_distance_all.append(dic_loci_location_head_nearest_lad_distance_point_cal[loci_location_i])

    location_matrix=np.array(sorted_var_location[use_begin:use_begin+use_point_num],dtype='float32')
    h_matrix=np.array(sorted_h[use_begin:use_begin+use_point_num],dtype='float32')
    sorted_nearest_lad_distance=np.array(sorted_nearest_lad_distance_all[use_begin:use_begin+use_point_num])
    # r_matrix=np.array(sorted_r[use_begin:use_begin+use_point_num],dtype='float32')

    location_matrix=location_matrix.reshape([1,use_point_num])
    h_matrix=h_matrix.reshape([1,use_point_num])
    # r_matrix=r_matrix.reshape([1,use_point_num])

    # cos_phi_sq=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=1.0, dtype=tf.float32, seed=1234))
    # sin_phi_sq=tf.ones([1,use_point_num])-cos_phi_sq

    phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=36))
#     phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=2020))
    phi_place=tf.placeholder(tf.float32,[1,use_point_num])
    phi_assign_opt=tf.assign(phi,phi_place)

    # r_matrix=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=R-0.5, dtype=tf.float32, seed=63))
    r_matrix=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=R-0.5, dtype=tf.float32, seed=63))
    r_place=tf.placeholder(tf.float32,[1,use_point_num])
    r_assign_opt=tf.assign(r_matrix,r_place)

    # phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=1995))
    # phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=1234))
    cos_phi=tf.cos(phi)
    sin_phi=tf.sin(phi)

    #calculate r*cos\phi
    rc=tf.multiply(r_matrix,cos_phi) #[1,use_point_num]
    M_rc=tf.convert_to_tensor(tf.tile(rc,(use_point_num,1))) #[use_point_num,use_point_num]
    M_rcT=tf.transpose(M_rc)
    M_rc_res_M=M_rc-M_rcT

    #calculate (r*sin\phi-r*sin\phi)**2
    rs=tf.multiply(r_matrix,sin_phi) #[1,use_point_num]
    M_rs=tf.convert_to_tensor(tf.tile(rs,(use_point_num,1))) #[use_point_num,use_point_num]
    M_rsT=tf.transpose(M_rs)
    M_rs_res_M=M_rs-M_rsT

    #calculate (h-h)**2
    h_r=np.tile(h_matrix,(use_point_num,1)) #[use_point_num,use_point_num]
    h_rT=np.transpose(h_r)
    h_r_res_M=h_r-h_rT

    #calculate (distance_12)**2
    L_r=np.tile(location_matrix,(use_point_num,1)) #[use_point_num,use_point_num]
    L_rT=np.transpose(L_r)
    # dis=get_spatial_distance(np.abs(L_r-L_rT))
    dis=get_spatial_distance_array(np.abs(L_r-L_rT))
    dis_square=np.square(np.array(dis))

    #calculate (M1-M2)**2-M12**2
    # M_rc_res=tf.reduce_sum(tf.math.square(M_rc_res_M))
    # M_rs_res=tf.reduce_sum(tf.math.square(M_rs_res_M))
    # h_r_res=np.sum(np.square(h_r_res_M))
    M_rc_res_square=tf.math.square(M_rc_res_M)
    M_rs_res_square=tf.math.square(M_rs_res_M)
    h_r_res_square=np.square(h_r_res_M)
    Res_matrix=tf.abs(M_rc_res_square+M_rs_res_square+h_r_res_square-dis_square)
    Res_matrix_original=tf.abs(M_rc_res_square+M_rs_res_square+h_r_res_square-dis_square)
    Res_matrix=tf.multiply(tf.divide(1,(np.abs(L_r-L_rT+constant_res)/limit_res)**(limit_power_res)),Res_matrix)

    #! r_z should=(h**2+r**2)**2
    #calculate abs( (R-r)**2 - ( min LAD distance )**2  )
    r_sphere_matrix=tf.math.sqrt(tf.math.square(r_matrix)+tf.math.square(h_matrix))
    R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_sphere_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    # R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    R_minus_matrix_square=tf.math.square(R_minus_matrix)

    # LAD_distance_array=get_spatial_distance(np.abs(sorted_nearest_lad_distance))
    LAD_distance_array=get_spatial_distance_list(np.abs(sorted_nearest_lad_distance))
    LAD_distance_matrix_square=np.square(np.tile(LAD_distance_array,(use_point_num,1))) #[use_point_num,use_point_num]
    LAD_distance_matrix_square=np.array(LAD_distance_matrix_square,dtype='float32')
    limit_r_by_LAD_square=tf.abs(R_minus_matrix_square-LAD_distance_matrix_square)
    limit_r_by_LAD_square=tf.convert_to_tensor(limit_r_by_LAD_square,dtype=tf.float32)
    limit_r_by_LAD_square_original=tf.convert_to_tensor(limit_r_by_LAD_square,dtype=tf.float32)
    limit_r_by_LAD_square=tf.multiply(tf.divide(1,((np.array(sorted_nearest_lad_distance,dtype='float32')+constant_LAD)/limit_LAD)**(limit_power_LAD)),limit_r_by_LAD_square)

    #calculate cost funciton (all is *2 )
    cost_func=tf.reduce_sum(Res_matrix+limit_r_by_LAD_square)

    # steps=500000
    # lr=1e-5
    # global_step = tf.Variable(0,trainable=False)  
    # lr=tf.train.piecewise_constant(global_step, boundaries=lr_boundaries, values=lr_values)
    opt = tf.train.AdamOptimizer(learning_rate=lr)
    opt_op = opt.minimize(cost_func)

    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        phi_now=sess.run(phi)
        r_now=sess.run(r_matrix)
        # _,_=sess.run([parameter_A_assign_opt,parameter_C_assign_opt],feed_dict={parameter_A_place:a_for_parameter,parameter_C_place:c_for_parameter})
        for i in range(steps):
    #         phi_now%(2*math.pi)
    #         _=sess.run(phi_assign_opt,feed_dice={phi_place:phi_now})
            for i_r in range(len(r_now[0])):
                if r_now[0][i_r]**2 + h_matrix[0][i_r]**2 > R**2:
                    r_now[0][i_r]=(R**2-h_matrix[0][i_r]**2)**(0.5)
            _,_,_,loss_now,phi_now,r_now=sess.run([phi_assign_opt,r_assign_opt,opt_op,cost_func,phi,r_matrix],feed_dict={phi_place:phi_now%(2*math.pi),r_place:r_now})
            if i%10000==0:
                print("Step %d,loss_: %f"%(i,loss_now))
        loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now=sess.run([cost_func,Res_matrix,rc,rs,r_matrix])
                # print("Step %d,loss_: %f, phi_now:%s"%(i,loss_now,phi_now))
        M_rc_res_square_now,M_rs_res_square_now,limit_r_by_LAD_square_now=sess.run([M_rc_res_square,M_rs_res_square,limit_r_by_LAD_square])
        Res_matrix_original_now,limit_r_by_LAD_square_original_now=sess.run([Res_matrix_original,limit_r_by_LAD_square_original])
    
    return loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square,dis_square,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now
    # return loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix



# #################################################################################
# #################################################################################
# build the result for whole chr
# #################################################################################
# #################################################################################
# begin_5_r,begin_5_phi is the 5 overlap point from end of last use_begin.
def run_phi_with_LAD_with_overlap(steps,lr,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_cal,dic_loci_location_head_nearest_lad_distance_point_cal):

    var_list_location=[]
    for loci_location_i in dic_group_loci_location_h_cal:
        var_list_location.append(loci_location_i)

    sorted_var_location=np.sort(var_list_location)
    sorted_h=[]
    sorted_nearest_lad_distance_all=[]
    # sorted_r=[]
    for loci_location_i in sorted_var_location:
        sorted_h.append(dic_group_loci_location_h_cal[loci_location_i])
        # sorted_r.append(dic_group_loci_location_r_C[loci_location_i])
        sorted_nearest_lad_distance_all.append(dic_loci_location_head_nearest_lad_distance_point_cal[loci_location_i])

    location_matrix=np.array(sorted_var_location[use_begin:use_begin+use_point_num],dtype='float32')
    h_matrix=np.array(sorted_h[use_begin:use_begin+use_point_num],dtype='float32')
    sorted_nearest_lad_distance=np.array(sorted_nearest_lad_distance_all[use_begin:use_begin+use_point_num])
    # r_matrix=np.array(sorted_r[use_begin:use_begin+use_point_num],dtype='float32')

    location_matrix=location_matrix.reshape([1,use_point_num])
    h_matrix=h_matrix.reshape([1,use_point_num])
    # r_matrix=r_matrix.reshape([1,use_point_num])

    # cos_phi_sq=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=1.0, dtype=tf.float32, seed=1234))
    # sin_phi_sq=tf.ones([1,use_point_num])-cos_phi_sq

    phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=2020))
    phi_place=tf.placeholder(tf.float32,[1,use_point_num])
    phi_assign_opt=tf.assign(phi,phi_place)

    r_matrix=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=R-0.5, dtype=tf.float32, seed=1995))
    r_place=tf.placeholder(tf.float32,[1,use_point_num])
    r_assign_opt=tf.assign(r_matrix,r_place)

    # phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=1995))
    # phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=1234))
    cos_phi=tf.cos(phi)
    sin_phi=tf.sin(phi)

    #calculate r*cos\phi
    rc=tf.multiply(r_matrix,cos_phi) #[1,use_point_num]
    M_rc=tf.convert_to_tensor(tf.tile(rc,(use_point_num,1))) #[use_point_num,use_point_num]
    M_rcT=tf.transpose(M_rc)
    M_rc_res_M=M_rc-M_rcT

    #calculate (r*sin\phi-r*sin\phi)**2
    rs=tf.multiply(r_matrix,sin_phi) #[1,use_point_num]
    M_rs=tf.convert_to_tensor(tf.tile(rs,(use_point_num,1))) #[use_point_num,use_point_num]
    M_rsT=tf.transpose(M_rs)
    M_rs_res_M=M_rs-M_rsT

    #calculate (h-h)**2
    h_r=np.tile(h_matrix,(use_point_num,1)) #[use_point_num,use_point_num]
    h_rT=np.transpose(h_r)
    h_r_res_M=h_r-h_rT

    #calculate (distance_12)**2
    L_r=np.tile(location_matrix,(use_point_num,1)) #[use_point_num,use_point_num]
    L_rT=np.transpose(L_r)
    # dis=get_spatial_distance(np.abs(L_r-L_rT))
    dis=get_spatial_distance_array(np.abs(L_r-L_rT))
    dis_square=np.square(np.array(dis))

    #calculate (M1-M2)**2-M12**2
    # M_rc_res=tf.reduce_sum(tf.math.square(M_rc_res_M))
    # M_rs_res=tf.reduce_sum(tf.math.square(M_rs_res_M))
    # h_r_res=np.sum(np.square(h_r_res_M))
    M_rc_res_square=tf.math.square(M_rc_res_M)
    M_rs_res_square=tf.math.square(M_rs_res_M)
    h_r_res_square=np.square(h_r_res_M)
    Res_matrix=tf.abs(M_rc_res_square+M_rs_res_square+h_r_res_square-dis_square)
    Res_matrix_original=tf.abs(M_rc_res_square+M_rs_res_square+h_r_res_square-dis_square)
    Res_matrix=tf.multiply(tf.divide(1,(np.abs(L_r-L_rT+constant_res)/limit_res)**(limit_power_res)),Res_matrix)

    #! r_z should=(h**2+r**2)**2
    #calculate abs( (R-r)**2 - ( min LAD distance )**2  )
    r_sphere_matrix=tf.math.sqrt(tf.math.square(r_matrix)+tf.math.square(h_matrix))
    R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_sphere_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    # R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    R_minus_matrix_square=tf.math.square(R_minus_matrix)
    # LAD_distance_array=get_spatial_distance(np.abs(sorted_nearest_lad_distance))
    LAD_distance_array=get_spatial_distance_list(np.abs(sorted_nearest_lad_distance))
    LAD_distance_matrix_square=np.square(np.tile(LAD_distance_array,(use_point_num,1))) #[use_point_num,use_point_num]
    LAD_distance_matrix_square=np.array(LAD_distance_matrix_square,dtype='float32')
    limit_r_by_LAD_square=tf.abs(R_minus_matrix_square-LAD_distance_matrix_square)
    limit_r_by_LAD_square=tf.convert_to_tensor(limit_r_by_LAD_square,dtype=tf.float32)
    limit_r_by_LAD_square_original=tf.convert_to_tensor(limit_r_by_LAD_square,dtype=tf.float32)
    limit_r_by_LAD_square=tf.multiply(tf.divide(1,((np.array(sorted_nearest_lad_distance,dtype='float32')+constant_LAD)/limit_LAD)**(limit_power_LAD)),limit_r_by_LAD_square)

    #calculate cost funciton (all is *2 )
    cost_func=tf.reduce_sum(Res_matrix+limit_r_by_LAD_square)

    # steps=500000
    # lr=1e-5
    # global_step = tf.Variable(0,trainable=False)  
    # lr=tf.train.piecewise_constant(global_step, boundaries=lr_boundaries, values=lr_values)
    opt = tf.train.AdamOptimizer(learning_rate=lr)
    opt_op = opt.minimize(cost_func)

    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        phi_now=sess.run(phi)
        phi_now[0][0:len(begin_5_phi)]=begin_5_phi
        r_now=sess.run(r_matrix)
        r_now[0][0:len(begin_5_r)]=begin_5_r
        # _,_=sess.run([parameter_A_assign_opt,parameter_C_assign_opt],feed_dict={parameter_A_place:a_for_parameter,parameter_C_place:c_for_parameter})
        for i in range(steps):
    #         phi_now%(2*math.pi)
    #         _=sess.run(phi_assign_opt,feed_dice={phi_place:phi_now})
            for i_r in range(len(r_now[0])):
                # if r_now[0][i_r]>R:
                #     r_now[0][i_r]=R
                if r_now[0][i_r]**2 + h_matrix[0][i_r]**2 > R**2:
                    r_now[0][i_r]=(R**2-h_matrix[0][i_r]**2)**(0.5)
            _,_,_,loss_now,phi_now,r_now=sess.run([phi_assign_opt,r_assign_opt,opt_op,cost_func,phi,r_matrix],feed_dict={phi_place:phi_now%(2*math.pi),r_place:r_now})
            phi_now[0][0:len(begin_5_phi)]=begin_5_phi
            r_now[0][0:len(begin_5_r)]=begin_5_r
            if i%10000==0:
                print("Step %d,loss_: %f"%(i,loss_now))
        loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now=sess.run([cost_func,Res_matrix,rc,rs,r_matrix])
                # print("Step %d,loss_: %f, phi_now:%s"%(i,loss_now,phi_now))
        M_rc_res_square_now,M_rs_res_square_now,limit_r_by_LAD_square_now=sess.run([M_rc_res_square,M_rs_res_square,limit_r_by_LAD_square])
        Res_matrix_original_now,limit_r_by_LAD_square_original_now=sess.run([Res_matrix_original,limit_r_by_LAD_square_original])

    return loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square,dis_square,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now
    # return loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix


# # fixed_point_loci_location,fixed_point_r,fixed_point_phi is the fixed point from hierarchically calculation
# ATTENTION: the fixed_point_loci_location is all in use_loci_locaiton
def run_phi_with_LAD_with_given_point_fixed_point(steps,lr,use_point_num,use_loci_locaiton,fixed_point_loci_location,fixed_point_r,fixed_point_phi,dic_group_loci_location_h_cal,dic_loci_location_head_nearest_lad_distance_point_cal):

    var_list_location=[]
    for loci_location_i in use_loci_locaiton:
        var_list_location.append(loci_location_i)

    sorted_var_location=np.sort(var_list_location)
    sorted_h=[]
    sorted_nearest_lad_distance_all=[]
    for loci_location_i in sorted_var_location:
        sorted_h.append(dic_group_loci_location_h_cal[loci_location_i])
        # sorted_r.append(dic_group_loci_location_r_C[loci_location_i])
        sorted_nearest_lad_distance_all.append(dic_loci_location_head_nearest_lad_distance_point_cal[loci_location_i])
    
    sorted_fixed_loci_index=[] # the index of fixed loci in sorted_var_location,output r etc
    for loci_location_i in fixed_point_loci_location:
        sorted_fixed_loci_index.append(np.where(sorted_var_location==loci_location_i)[0][0])

    location_matrix=np.array(sorted_var_location,dtype='float32')
    h_matrix=np.array(sorted_h,dtype='float32')
    sorted_nearest_lad_distance=np.array(sorted_nearest_lad_distance_all)

    location_matrix=location_matrix.reshape([1,use_point_num])
    h_matrix=h_matrix.reshape([1,use_point_num])

    phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=2020))
    phi_place=tf.placeholder(tf.float32,[1,use_point_num])
    phi_assign_opt=tf.assign(phi,phi_place)

    r_matrix=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=R, dtype=tf.float32, seed=1995))
    r_place=tf.placeholder(tf.float32,[1,use_point_num])
    r_assign_opt=tf.assign(r_matrix,r_place)

    cos_phi=tf.cos(phi)
    sin_phi=tf.sin(phi)

    #calculate r*cos\phi
    rc=tf.multiply(r_matrix,cos_phi) #[1,use_point_num]
    M_rc=tf.convert_to_tensor(tf.tile(rc,(use_point_num,1))) #[use_point_num,use_point_num]
    M_rcT=tf.transpose(M_rc)
    M_rc_res_M=M_rc-M_rcT

    #calculate (r*sin\phi-r*sin\phi)**2
    rs=tf.multiply(r_matrix,sin_phi) #[1,use_point_num]
    M_rs=tf.convert_to_tensor(tf.tile(rs,(use_point_num,1))) #[use_point_num,use_point_num]
    M_rsT=tf.transpose(M_rs)
    M_rs_res_M=M_rs-M_rsT

    #calculate (h-h)**2
    h_r=np.tile(h_matrix,(use_point_num,1)) #[use_point_num,use_point_num]
    h_rT=np.transpose(h_r)
    h_r_res_M=h_r-h_rT

    #calculate (distance_12)**2
    L_r=np.tile(location_matrix,(use_point_num,1)) #[use_point_num,use_point_num]
    L_rT=np.transpose(L_r)
    # dis=get_spatial_distance(np.abs(L_r-L_rT))
    dis=get_spatial_distance_array(np.abs(L_r-L_rT))
    dis_square=np.square(np.array(dis))

    #calculate (M1-M2)**2-M12**2
    M_rc_res_square=tf.math.square(M_rc_res_M)
    M_rs_res_square=tf.math.square(M_rs_res_M)
    h_r_res_square=np.square(h_r_res_M)
    Res_matrix=tf.abs(M_rc_res_square+M_rs_res_square+h_r_res_square-dis_square)
    Res_matrix_original=tf.abs(M_rc_res_square+M_rs_res_square+h_r_res_square-dis_square)
    Res_matrix=tf.multiply(tf.divide(1,(np.abs(L_r-L_rT+constant_res)/limit_res)**(limit_power_res)),Res_matrix)
    
    #! r_z should=(h**2+r**2)**2
    #calculate abs( (R-r)**2 - ( min LAD distance )**2  )
    r_sphere_matrix=tf.math.sqrt(tf.math.square(r_matrix)+tf.math.square(h_matrix))
    R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_sphere_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    # R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    R_minus_matrix_square=tf.math.square(R_minus_matrix)
    # LAD_distance_array=get_spatial_distance(np.abs(sorted_nearest_lad_distance))
    LAD_distance_array=get_spatial_distance_list(np.abs(sorted_nearest_lad_distance))
    LAD_distance_matrix_square=np.square(np.tile(LAD_distance_array,(use_point_num,1))) #[use_point_num,use_point_num]
    LAD_distance_matrix_square=np.array(LAD_distance_matrix_square,dtype='float32')
    limit_r_by_LAD_square=tf.abs(R_minus_matrix_square-LAD_distance_matrix_square)
    limit_r_by_LAD_square=tf.convert_to_tensor(limit_r_by_LAD_square,dtype=tf.float32)
    limit_r_by_LAD_square_original=tf.convert_to_tensor(limit_r_by_LAD_square,dtype=tf.float32)
    limit_r_by_LAD_square=tf.multiply(tf.divide(1,((np.array(sorted_nearest_lad_distance,dtype='float32')+constant_LAD)/limit_LAD)**(limit_power_LAD)),limit_r_by_LAD_square)

    #calculate cost funciton (all is *2 )
    cost_func=tf.reduce_sum(Res_matrix+limit_r_by_LAD_square)

    # global_step = tf.Variable(0,trainable=False)  
    # lr=tf.train.piecewise_constant(global_step, boundaries=lr_boundaries, values=lr_values)
    opt = tf.train.AdamOptimizer(learning_rate=lr)
    opt_op = opt.minimize(cost_func)

    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        #sorted_fixed_loci_index, fixed_point_loci_location,fixed_point_r,fixed_point_phi,
        phi_now=sess.run(phi)
        r_now=sess.run(r_matrix)
        for i_fixed_point in range(len(sorted_fixed_loci_index)):
            phi_now[0][sorted_fixed_loci_index[i_fixed_point]]=fixed_point_phi[i_fixed_point]
            r_now[0][sorted_fixed_loci_index[i_fixed_point]]=fixed_point_r[i_fixed_point]
        # _,_=sess.run([parameter_A_assign_opt,parameter_C_assign_opt],feed_dict={parameter_A_place:a_for_parameter,parameter_C_place:c_for_parameter})
        for i in range(steps):
            for i_r in range(len(r_now[0])):
                if r_now[0][i_r]**2 + h_matrix[0][i_r]**2 > R**2:
                    r_now[0][i_r]=(R**2-h_matrix[0][i_r]**2)**(0.5)
            _,_,_,loss_now,phi_now,r_now=sess.run([phi_assign_opt,r_assign_opt,opt_op,cost_func,phi,r_matrix],feed_dict={phi_place:phi_now%(2*math.pi),r_place:r_now})
            for i_fixed_point in range(len(sorted_fixed_loci_index)):
                phi_now[0][sorted_fixed_loci_index[i_fixed_point]]=fixed_point_phi[i_fixed_point]
                r_now[0][sorted_fixed_loci_index[i_fixed_point]]=fixed_point_r[i_fixed_point]
            if i%10000==0:
                print("Step %d,loss_: %f"%(i,loss_now))
        loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now=sess.run([cost_func,Res_matrix,rc,rs,r_matrix])
                # print("Step %d,loss_: %f, phi_now:%s"%(i,loss_now,phi_now))
        M_rc_res_square_now,M_rs_res_square_now,limit_r_by_LAD_square_now=sess.run([M_rc_res_square,M_rs_res_square,limit_r_by_LAD_square])
        Res_matrix_original_now,limit_r_by_LAD_square_original_now=sess.run([Res_matrix_original,limit_r_by_LAD_square_original])

    return loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square,dis_square,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now
    # return loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix

# #################################################################################
# run the 2 method
# #################################################################################

if method=="overlap":
    # run all the loci for my data (for method with overlap)
    loss_now_all=[]
    Res_matrix_now_all=[]
    rc_now_all=[]
    rs_now_all=[]
    r_matrix_now_all=[]
    dis_square_all=[]
    L_r_all=[]
    L_rT_all=[]
    h_matrix_all=[]

    calculate_point_num=200
    use_point_num=205
    overlap_num=use_point_num-calculate_point_num

    if strain_plot=="C":
        circle_num=int(len(dic_group_loci_location_h_C.keys())/calculate_point_num)+1
    else:
        circle_num=int(len(dic_group_loci_location_h_S.keys())/calculate_point_num)+1

    begin_list = np.linspace(0, circle_num*calculate_point_num, circle_num+1)[:-1]

    for use_begin in begin_list:
        use_begin=int(use_begin)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi(use_point_num,use_begin,dic_group_loci_location_h_C)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(600000,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_C)
    #     del rc_now,rs_now,h_matrix
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(1800000,1e-5,use_point_num,use_begin,dic_group_loci_location_wSNP_C,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(810000,1e-5,use_point_num,use_begin,dic_group_loci_location_wSNP_C,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
        if use_begin != 0:
            if strain_plot=="C":
                if len(dic_group_loci_location_h_C.keys())-use_begin < use_point_num:
                    use_point_num=len(dic_group_loci_location_h_C.keys())-use_begin
            else:
                if len(dic_group_loci_location_h_S.keys())-use_begin < use_point_num:
                    use_point_num=len(dic_group_loci_location_h_S.keys())-use_begin
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(1,1e-5,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)    
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(810000,1e-5,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)    
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(18000,1e-3,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)    
            if strain_plot=="C":
                loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(my_step,my_lr,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)    
            else:
                loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(my_step,my_lr,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_S,dic_loci_location_head_nearest_lad_distance_point_S)    
            begin_5_r=r_matrix_now[0][-overlap_num:]
            begin_5_phi=phi_now[0][-overlap_num:]

            loss_now_all.append(loss_now)
            Res_matrix_now_all.append(Res_matrix_now[overlap_num:])
            rc_now_all.append(rc_now[0][overlap_num:])
            rs_now_all.append(rs_now[0][overlap_num:])
            r_matrix_now_all.append(r_matrix_now[0][overlap_num:])
            dis_square_all.append(dis_square[overlap_num:])
            L_r_all.append(L_r[overlap_num:])
            L_rT_all.append(L_rT[overlap_num:])
            h_matrix_all.append(h_matrix[0][overlap_num:])
        else:
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(1,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(720000,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(18000,1e-3,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            if strain_plot=="C":
                loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(my_step,my_lr,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            else:
                loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(my_step,my_lr,use_point_num,use_begin,dic_group_loci_location_h_S,dic_loci_location_head_nearest_lad_distance_point_S)

            begin_5_r=r_matrix_now[0][-overlap_num:]
            begin_5_phi=phi_now[0][-overlap_num:]

            loss_now_all.append(loss_now)
            Res_matrix_now_all.append(Res_matrix_now)
            rc_now_all.append(rc_now)
            rs_now_all.append(rs_now)
            r_matrix_now_all.append(r_matrix_now)
            dis_square_all.append(dis_square)
            L_r_all.append(L_r)
            L_rT_all.append(L_rT)
            h_matrix_all.append(h_matrix)
        print("use_begin is %s, loss_now is %s"%(use_begin,loss_now))

# overlap with 56 points
if method=="overlap56":
    # run all the loci for my data (for method with overlap)
    loss_now_all=[]
    Res_matrix_now_all=[]
    rc_now_all=[]
    rs_now_all=[]
    r_matrix_now_all=[]
    dis_square_all=[]
    L_r_all=[]
    L_rT_all=[]
    h_matrix_all=[]

    calculate_point_num=200
    use_point_num=256
    overlap_num=use_point_num-calculate_point_num

    circle_num=int(len(dic_group_loci_location_h_C.keys())/calculate_point_num)+1
    begin_list = np.linspace(0, circle_num*calculate_point_num, circle_num+1)[:-1]

    for use_begin in begin_list:
        use_begin=int(use_begin)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi(use_point_num,use_begin,dic_group_loci_location_h_C)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(600000,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_C)
    #     del rc_now,rs_now,h_matrix
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(1800000,1e-5,use_point_num,use_begin,dic_group_loci_location_wSNP_C,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(810000,1e-5,use_point_num,use_begin,dic_group_loci_location_wSNP_C,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
        if use_begin != 0:
            if len(dic_group_loci_location_h_C.keys())-use_begin < use_point_num:
                use_point_num=len(dic_group_loci_location_h_C.keys())-use_begin
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(8,1e-5,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)    
            loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(810000,1e-5,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)    
            begin_5_r=r_matrix_now[0][-overlap_num:]
            begin_5_phi=phi_now[0][-overlap_num:]

            loss_now_all.append(loss_now)
            Res_matrix_now_all.append(Res_matrix_now[overlap_num:])
            rc_now_all.append(rc_now[0][overlap_num:])
            rs_now_all.append(rs_now[0][overlap_num:])
            r_matrix_now_all.append(r_matrix_now[0][overlap_num:])
            dis_square_all.append(dis_square[overlap_num:])
            L_r_all.append(L_r[overlap_num:])
            L_rT_all.append(L_rT[overlap_num:])
            h_matrix_all.append(h_matrix[0][overlap_num:])
        else:
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(7,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(720000,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            begin_5_r=r_matrix_now[0][-overlap_num:]
            begin_5_phi=phi_now[0][-overlap_num:]

            loss_now_all.append(loss_now)
            Res_matrix_now_all.append(Res_matrix_now)
            rc_now_all.append(rc_now)
            rs_now_all.append(rs_now)
            r_matrix_now_all.append(r_matrix_now)
            dis_square_all.append(dis_square)
            L_r_all.append(L_r)
            L_rT_all.append(L_rT)
            h_matrix_all.append(h_matrix)
        print("use_begin is %s, loss_now is %s"%(use_begin,loss_now))

# run method 2: hierarchically 
if method=="hierarchical":
    # run all the loci for my data
    loss_now_all=[]
    Res_matrix_now_all=[]
    rc_now_all=[]
    rs_now_all=[]
    r_matrix_now_all=[]
    dis_square_all=[]
    L_r_all=[]
    L_rT_all=[]
    h_matrix_all=[]

    var_list_location=[]
    if strain_plot=="C":
        for loci_location_i in dic_group_loci_location_h_C:
            var_list_location.append(loci_location_i)
        loci_num=len(dic_group_loci_location_h_C.keys())
    else:
        for loci_location_i in dic_group_loci_location_h_S:
            var_list_location.append(loci_location_i)
        loci_num=len(dic_group_loci_location_h_S.keys())
    sorted_var_location=np.sort(var_list_location)

    # run first level
    # loci_num=len(dic_group_loci_location_h_C.keys())
    #each use_point_num(200) with about 4 point fixed
    first_level_loci_index=[int(x) for x in np.linspace(0,loci_num-1,int(loci_num/50))]
    first_level_loci_location=[sorted_var_location[x] for x in first_level_loci_index]
    # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_given_point_fixed_point(810000,1e-5,len(first_level_loci_location),first_level_loci_location,[],[],[],dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    # loss_first_level,Res_matrix_first_level,rc_first_level,rs_first_level,phi_first_level,r_matrix_first_level,dis_square_first_level,L_r_first_level,L_rT_first_level,h_matrix_first_level,M_rc_res_square_first_level,M_rs_res_square_first_level,h_r_res_square_first_level,dis_square_first_level,limit_r_by_LAD_square_first_level,Res_matrix_original_first_level,limit_r_by_LAD_square_original_first_level=run_phi_with_LAD_with_given_point_fixed_point(18000,1e-3,len(first_level_loci_location),first_level_loci_location,[],[],[],dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    if strain_plot=="C":
        loss_first_level,Res_matrix_first_level,rc_first_level,rs_first_level,phi_first_level,r_matrix_first_level,dis_square_first_level,L_r_first_level,L_rT_first_level,h_matrix_first_level,M_rc_res_square_first_level,M_rs_res_square_first_level,h_r_res_square_first_level,dis_square_first_level,limit_r_by_LAD_square_first_level,Res_matrix_original_first_level,limit_r_by_LAD_square_original_first_level=run_phi_with_LAD_with_given_point_fixed_point(my_step,my_lr,len(first_level_loci_location),first_level_loci_location,[],[],[],dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    else:
        loss_first_level,Res_matrix_first_level,rc_first_level,rs_first_level,phi_first_level,r_matrix_first_level,dis_square_first_level,L_r_first_level,L_rT_first_level,h_matrix_first_level,M_rc_res_square_first_level,M_rs_res_square_first_level,h_r_res_square_first_level,dis_square_first_level,limit_r_by_LAD_square_first_level,Res_matrix_original_first_level,limit_r_by_LAD_square_original_first_level=run_phi_with_LAD_with_given_point_fixed_point(my_step,my_lr,len(first_level_loci_location),first_level_loci_location,[],[],[],dic_group_loci_location_h_S,dic_loci_location_head_nearest_lad_distance_point_S)

    # loss_first_level,Res_matrix_first_level,rc_first_level,rs_first_level,phi_first_level,r_matrix_first_level,dis_square_first_level,L_r_first_level,L_rT_first_level,h_matrix_first_level,M_rc_res_square_first_level,M_rs_res_square_first_level,h_r_res_square_first_level,dis_square_first_level,limit_r_by_LAD_square_first_level,Res_matrix_original_first_level,limit_r_by_LAD_square_original_first_level=run_phi_with_LAD_with_given_point_fixed_point(540000,1e-5,len(first_level_loci_location),first_level_loci_location,[],[],[],dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)

    # run the small level for method2
    use_point_num=200
    circle_num=int(loci_num/use_point_num)+1
    begin_list = np.linspace(0, circle_num*use_point_num, circle_num+1)[:-1]
    # begin_list =[0] + begin_list

    # run_phi_with_LAD_with_fixed_point(steps,lr,use_point_num,use_begin,fixed_point_loci_location,fixed_point_r,fixed_point_phi,dic_group_loci_location_h_cal,dic_loci_location_head_nearest_lad_distance_point_cal)

    # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(810000,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    # run_phi_with_LAD_with_given_point_fixed_point(steps,lr,use_loci_locaiton,fixed_point_loci_location,fixed_point_r,fixed_point_phi,dic_group_loci_location_h_cal,dic_loci_location_head_nearest_lad_distance_point_cal)

    for use_begin in begin_list:
        use_begin=int(use_begin)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi(use_point_num,use_begin,dic_group_loci_location_h_C)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(600000,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_C)
    #     del rc_now,rs_now,h_matrix
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(1800000,1e-5,use_point_num,use_begin,dic_group_loci_location_wSNP_C,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
    #     loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix=run_phi_with_LAD(810000,1e-5,use_point_num,use_begin,dic_group_loci_location_wSNP_C,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
        if strain_plot=="C":
            if len(dic_group_loci_location_h_C.keys())-use_begin < use_point_num:
                use_point_num=len(dic_group_loci_location_h_C.keys())-use_begin
        else:
            if len(dic_group_loci_location_h_S.keys())-use_begin < use_point_num:
                use_point_num=len(dic_group_loci_location_h_S.keys())-use_begin
        
        use_loci_locaiton=sorted_var_location[use_begin:use_begin+use_point_num]

        fixed_point_loci_location=[]
        fixed_point_r=[]
        fixed_point_phi=[]
        for first_level_loci_location_index in range(len(first_level_loci_location)):
            if first_level_loci_location[first_level_loci_location_index] in use_loci_locaiton:
                fixed_point_loci_location.append(first_level_loci_location[first_level_loci_location_index])
                fixed_point_r.append(r_matrix_first_level[0][first_level_loci_location_index])
                fixed_point_phi.append(phi_first_level[0][first_level_loci_location_index])

        print("fixed point number in this round is %s"%(len(fixed_point_loci_location)))

        # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_given_point_fixed_point(18000,1e-3,use_point_num,use_loci_locaiton,fixed_point_loci_location,fixed_point_r,fixed_point_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
        if strain_plot=="C":
            loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_given_point_fixed_point(my_step,my_lr,use_point_num,use_loci_locaiton,fixed_point_loci_location,fixed_point_r,fixed_point_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
        else:
            loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_given_point_fixed_point(my_step,my_lr,use_point_num,use_loci_locaiton,fixed_point_loci_location,fixed_point_r,fixed_point_phi,dic_group_loci_location_h_S,dic_loci_location_head_nearest_lad_distance_point_S)

        # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(810000,1e-5,use_point_num,use_begin,begin_5_r,begin_5_phi,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)

        loss_now_all.append(loss_now)
        Res_matrix_now_all.append(Res_matrix_now)
        rc_now_all.append(rc_now)
        rs_now_all.append(rs_now)
        r_matrix_now_all.append(r_matrix_now)
        dis_square_all.append(dis_square)
        L_r_all.append(L_r)
        L_rT_all.append(L_rT)
        h_matrix_all.append(h_matrix)
        
        print("use_begin is %s, loss_now is %s"%(use_begin,loss_now))
        

# #################################################################################
# whole genome in whole genome 
# 3d plot the genome_distance with spatial diatance
# the scatter plot
# #################################################################################
# plot the 3D fig for whole chr all loci
rc_all=np.array([])
rs_all=np.array([])
h_all=np.array([])
for i in range(len(rc_now_all)):
    rc_all=np.append(rc_all,rc_now_all[i])
    rs_all=np.append(rs_all,rs_now_all[i])
    h_all=np.append(h_all,h_matrix_all[i])


# #################################################################################
# calculate spatial distance
# #################################################################################
spatial_distance=[]
for i in range(len(rc_all)):
    spatial_distance.append(np.linalg.norm([rc_all[i],rs_all[i],h_all[i]]))

print("mean spatial distance is %s"%(np.mean(spatial_distance)))
print("std of spatial distance is %s"%(np.std(spatial_distance)))
# with open(png_dir+chr_name+"_10kb_spatial_distance_"+method+".txt",'w') as f:
#     f.write("mean spatial distance is "+str(np.mean(spatial_distance)))
#     f.write("/n")
#     f.write("std of spatial distance is "+str(np.std(spatial_distance)))

# #################################################################################
# whole genome in whole genome 
# sns plot the genome_distance with spatial diatance
# the scatter plot
# #################################################################################
# plot the spatial distance with different genome distance

var_list_location=[]
if strain_plot=="C":
    for loci_location_i in dic_group_loci_location_h_C:
        var_list_location.append(loci_location_i)
else:
    for loci_location_i in dic_group_loci_location_h_S:
        var_list_location.append(loci_location_i)

sorted_var_location=np.sort(var_list_location)

spatial_distance_all=[]
theoritical_spatial_distance_all=[]
genome_distance_plot_all=[]
for i in range(len(Res_matrix_now)):
    for j in range(i):
#         spatial_distance.append(np.sqrt(Res_matrix_now[i,j]))
        spatial_distance_all.append(np.linalg.norm([rc_all[i]-rc_all[j],rs_all[i]-rs_all[j],h_all[i]-h_all[j]]))
        # theoritical_spatial_distance.append(np.sqrt(dis_square_all[i][j]))
        theoritical_spatial_distance_all.append(get_spatial_distance(sorted_var_location[i]-sorted_var_location[j]))
        genome_distance_plot_all.append(np.abs(sorted_var_location[i]-sorted_var_location[j]))

distance_df_whole_genome=pd.DataFrame({"spatial_distance":np.array(spatial_distance_all+theoritical_spatial_distance_all),
                     "genome_distance_plot":np.array(genome_distance_plot_all+genome_distance_plot_all),
                     "label":np.array(["real"]*len(spatial_distance_all)+["theoritical"]*len(theoritical_spatial_distance_all))
                    })

print("plot gsline ...")
fig = plt.figure()
ax = plt.axes()
ax = sns.lineplot(x="genome_distance_plot", y="spatial_distance", hue="label",ci="sd", data=distance_df_whole_genome)

ax.set_xlabel('genome distance')
plt.title("spatial distance with different genome distance")
ax.set_ylabel('spatial distance')

fig.savefig(png_dir+chr_name+"_10kb_gsline_"+method+".eps",dpi=600,format='eps')
plt.close()
# plt.show()



# #################################################################################
# distribution plot: plot 3D genome spatial distance with DNA content
# #################################################################################
radial_distance_all=[]
for i in range(len(rc_all)):
    ra=np.linalg.norm([rc_all[i],rs_all[i],h_all[i]])
    # if ra<R:
    radial_distance_all.append(ra)

print("plot DNA_content ...")
fig = plt.figure()
ax = plt.axes()
sns.distplot(radial_distance_all, color="m")

plt.xlim([0,5])
plt.ylim([0,10])
ax.set_xlabel('radial distance')
plt.title("DNA content with different radial distance")
ax.set_ylabel('DNA content')
# plt.xlim([0,loci_location_range])
fig.savefig(png_dir+chr_name+"_10kb_DNA_content_"+method+".eps",dpi=600,format='eps')
plt.close()


# #################################################################################
# whole genome in whole genome 
# sns plot the genome_distance with spatial diatance
# the scatter plot
# #################################################################################
# plot the 3D fig for whole chr all loci 

print("plot loci_3D ...")
fig = plt.figure()
ax = plt.axes(projection='3d')
# ax.scatter3D(r_matrix*np.cos(phi_now),r_matrix*np.sin(phi_now),h_matrix,c='r', label='After calculation')
ax.scatter3D(rc_all,rs_all,h_all,c=sorted_var_location[0:len(rc_all)], label='After calculation')
ax.legend(loc='best')
ax.set_xlabel('x', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('y', fontdict={'size': 15, 'color': 'red'})
ax.set_zlabel('z', fontdict={'size': 15, 'color': 'red'})
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(-5, 5)

ax.view_init(10, 30)
fig.savefig(png_dir+chr_name+"_10kb_loci_3D_"+method+"_1.eps",dpi=600,format='eps')
ax.view_init(10, -30)
fig.savefig(png_dir+chr_name+"_10kb_loci_3D_"+method+"_2.eps",dpi=600,format='eps')
ax.view_init(80, 10)
fig.savefig(png_dir+chr_name+"_10kb_loci_3D_"+method+"_3.eps",dpi=600,format='eps')
plt.close()

# #################################################################################
# save whole genome
# save result for each [loci location] with [spatical coordinate:xyz]
# the scatter plot
# #################################################################################
# strain_plot="C"
dir_save=png_dir+chr_name+"_10kb_location_data_"+method+".txt"
print(dir_save)
with open(dir_save,"w") as f:
    for location_index in range(len(rc_all)):
        f.write(str(sorted_var_location[location_index])+"\t"+str(rc_all[location_index])+"\t"+str(rs_all[location_index])+"\t"+str(h_all[location_index]))
        f.write("\n")


# #################################################################################
# transfer 3D structure ie [distance matrix] to [contact probability matrix]
# #################################################################################
def cal_distance_matrix(resolution,sorted_var_location,rc_all,rs_all,h_all):
    min_loc_head=np.min(sorted_var_location)
    max_loc_head=np.max(sorted_var_location)
    matrix_length=int((max_loc_head-min_loc_head)/resolution+1)
    whole_distance_matrix=np.ones([matrix_length,matrix_length])*20
    # for location_head in range(min_loc_head,max_loc_head+1,resolution):
    #     if location_head in sorted_var_location:
    #         location_index=np.where(np.array(sorted_var_location)==location_head)
    for location_caputured_index in range(len(rc_all)):
        if (sorted_var_location[location_caputured_index]-min_loc_head)%resolution !=0:
            print("error in %s"%(sorted_var_location[location_caputured_index]))
            break
        else:
            location_in_whole_matrix=int((sorted_var_location[location_caputured_index]-min_loc_head)/resolution)
#             for location_caputured_index_other in range(location_caputured_index):
            for location_caputured_index_other in range(len(rc_all)):
                location_other_in_whole_matrix=int((sorted_var_location[location_caputured_index_other]-min_loc_head)/resolution)
                loci_1=np.array([rc_all[location_caputured_index],rs_all[location_caputured_index],h_all[location_caputured_index]])
                loci_2=np.array([rc_all[location_caputured_index_other],rs_all[location_caputured_index_other],h_all[location_caputured_index_other]])
                whole_distance_matrix[location_in_whole_matrix,location_other_in_whole_matrix]=np.linalg.norm(loci_1-loci_2)
                whole_distance_matrix[location_other_in_whole_matrix,location_in_whole_matrix]=np.linalg.norm(loci_1-loci_2)
#                   if location_caputured_index==location_caputured_index_other and whole_distance_matrix[location_caputured_index,location_other_in_whole_matrix] !=0:
#                     print("ERROR! location_caputured_index %s, loci_1 %s, loci_2 %s"%(location_caputured_index,loci_1,loci_2))
#                     break
    return whole_distance_matrix

whole_distance_matrix_cal=cal_distance_matrix(resolution,sorted_var_location,rc_all,rs_all,h_all)


# plot heat map for genome_distance_matrix
print("plot genome_distance_matrix ...")
f, ax = plt.subplots(figsize=(16, 16))
df = pd.DataFrame(whole_distance_matrix_cal)
sns.heatmap(df, annot=False, ax=ax,cmap = 'RdYlBu_r',cbar=False)
ax.set_xlabel('genome location')
plt.title("genome distance matrix")
ax.set_ylabel('genome location')
# plt.savefig(png_dir+chr_name+"_10kb_genome_distance_matrix_"+method+".eps",dpi=600,format='eps')
# plt.savefig(png_dir+chr_name+"_10kb_genome_distance_matrix_"+method+".eps",format='eps')
plt.savefig(png_dir+chr_name+"_10kb_genome_distance_matrix_"+method+"_notnorm.png")
plt.close()

# plot heat map for contact_probability_matrix
print("plot contact_probability_matrix ...")
f, ax = plt.subplots(figsize=(16, 16))
whole_distance_matrix_cal_paint=whole_distance_matrix_cal+9e-1
# whole_distance_matrix_cal_paint=whole_distance_matrix_cal+1e-9
df = pd.DataFrame((whole_distance_matrix_cal_paint)**(-3))
sns.heatmap(df, annot=False, ax=ax,cmap = 'RdYlBu_r',cbar=False)
ax.set_xlabel('genome location')
plt.title("contact probability matrix")
ax.set_ylabel('genome location')
# plt.savefig(png_dir+chr_name+"_10kb_contact_probability_matrix_"+method+".eps",dpi=600,format='eps')
# plt.savefig(png_dir+chr_name+"_10kb_contact_probability_matrix_"+method+".eps",format='eps')
plt.savefig(png_dir+chr_name+"_10kb_contact_probability_matrix_"+method+"_notnorm.png")
plt.close()



# #################################################################################
# use CpG to color the 3D model
# #################################################################################
# def load_cpg_data(chr_name,resolution):
#     dir_cpg=prefix+"/ygli/gam_paper_data/F123_cpg/GSM1027571_DNA_CpG_methcounts_E14_serum_LIF.bedGraph"
# #     dir_cpg=prefix+"/ygli/gam_paper_data/F123_cpg/esc_DNA_methylation.txt"
#     dir_cpg=prefix+"/ygli/gam_paper_data/mm9_cpg/cpg_calcualtion_result/"+chr_name+"_cpg.txt"
#     chr_needed=[chr_name]
#     chr_location_cpg={}
#     with open(dir_cpg,'r') as f:
#         f.readline()
#         all_lines=f.readlines()
#         for each_line in all_lines:
#             chr_name_this_l,pos_begin,pos_end,cpg_this_l=each_line.replace("\n","").split("\t")[0:4]
#             if chr_name_this_l in chr_needed:
#                 belong_loci=int(int(int(pos_begin)/resolution)*resolution)
# #                 print(belong_loci)
#                 if belong_loci in chr_location_cpg:
#                     chr_location_cpg[belong_loci].append(float(cpg_this_l))
#                 else:
#                     chr_location_cpg[belong_loci]=[float(cpg_this_l)]
#     return chr_location_cpg
# chr_location_cpg_cal=load_cpg_data(chr_name,resolution)

def load_cpg_data_strain(chr_name,resolution,strain_plot):
    if strain_plot=="C":
        dir_cpg=prefix+"/ygli/gam_paper_data/CAST_cpg/cpg_calcualtion_result_method_1/"+chr_name+"_cpg.txt"
    else:
        dir_cpg=prefix+"/ygli/gam_paper_data/129S1_cpg/cpg_calcualtion_result_method_1/"+chr_name+"_cpg.txt"
    chr_needed=[chr_name]
    chr_location_cpg={}
    with open(dir_cpg,'r') as f:
        f.readline()
        all_lines=f.readlines()
        for each_line in all_lines:
            chr_name_this_l,pos_begin,pos_end,cpg_this_l=each_line.replace("\n","").split("\t")[0:4]
            if chr_name_this_l in chr_needed:
                belong_loci=int(int(int(pos_begin)/resolution)*resolution)
#                 print(belong_loci)
                if belong_loci in chr_location_cpg:
                    chr_location_cpg[belong_loci].append(float(cpg_this_l))
                else:
                    chr_location_cpg[belong_loci]=[float(cpg_this_l)]
    return chr_location_cpg

chr_location_cpg_cal=load_cpg_data_strain(chr_name,resolution,strain_plot)

loci_location_cpg=[]
for loci_location_i in sorted_var_location:
    if loci_location_i in chr_location_cpg_cal:
#         percent=len(np.where(np.array(chr_location_cpg_cal[loci_location_i])>0.5)[0])/len(chr_location_cpg_cal[loci_location_i])
#         percent=len(np.where(np.array(chr_location_cpg_cal[loci_location_i]))[0])/resolution
        percent=np.mean(chr_location_cpg_cal[loci_location_i])
        loci_location_cpg.append(percent)
#         loci_location_cpg.append(np.mean(chr_location_cpg_cal[loci_location_i]))
    else:
        loci_location_cpg.append(0)
        


# with cpg as color: plot the 3D fig for whole chr all loci 
print("plot cpg_3D ...")
fig = plt.figure()
ax = plt.axes(projection='3d')
# ax.scatter3D(r_matrix*np.cos(phi_now),r_matrix*np.sin(phi_now),h_matrix,c='r', label='After calculation')
# sc=ax.scatter3D(rc_all,rs_all,h_all,c=loci_location_cpg[0:len(rc_all)],cmap='RdYlGn')
sc=ax.scatter3D(rc_all,rs_all,h_all,c=loci_location_cpg,cmap='RdYlGn_r')
plt.colorbar(sc)
ax.legend(loc='best')
ax.set_xlabel('x', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('y', fontdict={'size': 15, 'color': 'red'})
ax.set_zlabel('z', fontdict={'size': 15, 'color': 'red'})
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(-5, 5)

ax.view_init(10, 30)
fig.savefig(png_dir+chr_name+"_10kb_cpg_3D_"+method+"_1.eps",dpi=600,format='eps')
ax.view_init(10, -30)
fig.savefig(png_dir+chr_name+"_10kb_cpg_3D_"+method+"_2.eps",dpi=600,format='eps')
ax.view_init(80, 10)
fig.savefig(png_dir+chr_name+"_10kb_cpg_3D_"+method+"_3.eps",dpi=600,format='eps')
plt.close()

# # plot the relationship of spatial distance with cpg
spatial_distance_all=[]
cpg_all=[]
for i in range(len(rc_all)):
    spatial_distance_all.append(np.linalg.norm([rc_all[i],rs_all[i],h_all[i]]))
    cpg_all.append(loci_location_cpg[i])

cpg_df_whole_genome=pd.DataFrame({"spatial_distance":np.array(spatial_distance_all),
                     "cpg":np.array(cpg_all),
                     # "label":np.array(["real"]*len(spatial_distance_all)+["theoritical"]*len(theoritical_spatial_distance_all))
                    })

print("plot cpg_spatial_distance ...")
fig = plt.figure()
ax = plt.axes()
ax = sns.lineplot(x="spatial_distance", y="cpg", ci="sd",data=cpg_df_whole_genome)
ax.set_xlabel('spatial distance')
plt.title("cpg with different spatial distance")
ax.set_ylabel('cpg')
plt.savefig(png_dir+chr_name+"_10kb_cpg_spatial_distance_"+method+".eps",dpi=600,format='eps')
plt.close()


# # barplot: plot the relationship of spatial distance with cpg
# use 0.5 as bar[0,0.5,...]
bar_width=0.3
bar_num=int(4.5/bar_width)

spatial_distance_bar=[]
for i_loci in range(len(spatial_distance_all)):
    r_this=int(spatial_distance_all[i_loci]/bar_width)*bar_width
    spatial_distance_bar.append(r_this)

cpg_df_bar_whole_genome=pd.DataFrame({"spatial_distance":np.array(spatial_distance_bar),
                     "cpg":np.array(cpg_all)})

print("plot cpg_spatial_distance ...")
fig = plt.figure(figsize=(16, 4))
ax = plt.axes()
# ax = sns.lineplot(x="spatial_distance", y="cpg", ci=100,data=cpg_df_whole_genome)
ax = sns.lineplot(x="spatial_distance", y="cpg", ci="sd",data=cpg_df_bar_whole_genome)

ax.set_xlabel('spatial distance')
plt.title("cpg with different spatial distance")
ax.set_ylabel('cpg')
plt.savefig(png_dir+chr_name+"_10kb_cpg_spatial_distance_bar_"+method+".eps",dpi=600,format='eps')
plt.close()



# # plot the distribution of spatial distance of cpg
cpg_average=np.mean(cpg_all)
# spatial_distance_all=[]
cpg_high_r=[]
cpg_low_r=[]
for i in range(len(cpg_all)):
#     r_this=np.linalg.norm([rc_all[i],rs_all[i],h_all[i]])
    r_this=spatial_distance_all[i]
    if cpg_all[i]>2*cpg_average:
        cpg_high_r.append(r_this)
    else:
        cpg_low_r.append(r_this)

print("plot cpg_spatial_distance ...")
fig = plt.figure()
ax = plt.axes()
ax1 = sns.distplot(np.array(cpg_high_r), color="r")
ax2 = sns.distplot(np.array(cpg_low_r), color="b")
plt.legend(["high cpg","low cpg"])

ax.set_xlabel('spatial distance')
plt.title("cpg distribution")
ax.set_ylabel('probability')

plt.savefig(png_dir+chr_name+"_10kb_cpg_2_probability_"+method+".eps",dpi=600,format='eps')
plt.close()



# # plot the distribution of spatial distance of cpg
cpg_average=np.mean(cpg_all)
# spatial_distance_all=[]
cpg_high_r=[]
cpg_low_r=[]
for i in range(len(cpg_all)):
#     r_this=np.linalg.norm([rc_all[i],rs_all[i],h_all[i]])
    r_this=spatial_distance_all[i]
    if cpg_all[i]>cpg_average:
        cpg_high_r.append(r_this)
    else:
        cpg_low_r.append(r_this)

print("plot cpg_spatial_distance ...")
fig = plt.figure()
ax = plt.axes()
ax1 = sns.distplot(np.array(cpg_high_r), color="r")
ax2 = sns.distplot(np.array(cpg_low_r), color="b")
plt.legend(["high cpg","low cpg"])

ax.set_xlabel('spatial distance')
plt.title("cpg distribution")
ax.set_ylabel('probability')

plt.savefig(png_dir+chr_name+"_10kb_cpg_1_probability_"+method+".eps",dpi=600,format='eps')
plt.close()



#plot 2 pole result:  with cpg as color: plot the 3D fig for whole chr all loci 
print("plot cpg_3D ...")
fig = plt.figure()
ax = plt.axes(projection='3d')
# ax.scatter3D(r_matrix*np.cos(phi_now),r_matrix*np.sin(phi_now),h_matrix,c='r', label='After calculation')
sc=ax.scatter3D(rc_all,rs_all,h_all,c=np.sign(loci_location_cpg-cpg_average),cmap='RdYlGn_r')
plt.colorbar(sc)
ax.legend(loc='best')
ax.set_xlabel('x', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('y', fontdict={'size': 15, 'color': 'red'})
ax.set_zlabel('z', fontdict={'size': 15, 'color': 'red'})
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(-5, 5)

ax.view_init(10, 30)
fig.savefig(png_dir+chr_name+"_10kb_cpg_3D_binary_"+method+"_1.eps",dpi=600,format='eps')
ax.view_init(10, -30)
fig.savefig(png_dir+chr_name+"_10kb_cpg_3D_binary_"+method+"_2.eps",dpi=600,format='eps')
ax.view_init(80, 10)
fig.savefig(png_dir+chr_name+"_10kb_cpg_3D_binary_"+method+"_3.eps",dpi=600,format='eps')
plt.close()





# #plot 2 pole result: plot the relationship of spatial distance with cpg

# # barplot: plot the relationship of spatial distance with cpg
# use 0.5 as bar[0,0.5,...]
bar_width=0.3
bar_num=int(4.5/bar_width)

spatial_distance_bar_binary=[]
cpg_all_binary=[]
for i_loci in range(len(spatial_distance_all)):
    r_this=int(spatial_distance_all[i_loci]/bar_width)*bar_width
    spatial_distance_bar_binary.append(r_this)
    cpg_all_binary.append(np.sign(cpg_all[i_loci]-cpg_average))

cpg_df_bar_binary_whole_genome=pd.DataFrame({"spatial_distance":np.array(spatial_distance_bar_binary),
                     "cpg":np.array(cpg_all_binary)})

print("plot cpg_spatial_distance ...")
fig = plt.figure(figsize=(16, 4))
ax = plt.axes()
ax = sns.lineplot(x="spatial_distance", y="cpg", ci="sd",data=cpg_df_bar_binary_whole_genome)

ax.set_xlabel('spatial distance')
plt.title("cpg with different spatial distance")
ax.set_ylabel('cpg')

plt.savefig(png_dir+chr_name+"_10kb_cpg_spatial_distance_binary_"+method+".eps",dpi=600,format='eps')
plt.close()



# # plot the relationship of phi with cpg
spatial_phi_all=[]
cpg_all=[]
for i in range(len(rc_all)):
    r_this=np.linalg.norm([rc_all[i],rs_all[i],h_all[i]])
    phi=math.acos(rc_all[i]/r_this)
    if r_this<=R:
        spatial_phi_all.append(phi)
        cpg_all.append(loci_location_cpg[i])

cpg_phi_df_whole_genome=pd.DataFrame({"spatial_phi_all":np.array(spatial_phi_all),
                     "cpg":np.array(cpg_all)})

print("plot cpg_spatial_phi ...")
fig = plt.figure(figsize=(16, 4))
ax = plt.axes()
ax = sns.lineplot(x="spatial_phi_all", y="cpg", ci="sd",data=cpg_phi_df_whole_genome)
ax.set_xlabel('spatial phi')
plt.title("cpg with different spatial phi")
ax.set_ylabel('cpg')
plt.savefig(png_dir+chr_name+"_10kb_cpg_spatial_phi_"+method+".eps",dpi=600,format='eps')
plt.close()



# #################################################################################
# use Compartment AB to color the 3D model
# #################################################################################
def load_compartment_AB_data(chr_name,resolution):
    # dir_cpg=prefix+"/ygli/gam_paper_data/F123_cpg/GSM1027571_DNA_CpG_methcounts_E14_serum_LIF.bedGraph"
    dir_com=prefix+"/ygli/gam_paper_data/mes_hic/comp_AB_result_from_nij/"
    chr_location_com={}
    with open(dir_com+chr_name+"_40kb_AB_compartment.txt",'r') as f:
        all_lines=f.readlines()
        i_line=0
        for each_line in all_lines:
            pos_begin=i_line*4e+4
            comp_this_l=each_line.replace("\n","").split()[0]
            # belong_loci=int(int(pos_begin)/resolution)*resolution
            belong_loci_all=[pos_begin,pos_begin+1e+4,pos_begin+2e+4,pos_begin+3e+4]
            for belong_loci in belong_loci_all:
                belong_loci=int(belong_loci)
                if belong_loci in chr_location_com:
                    chr_location_com[belong_loci].append(float(comp_this_l))
                else:
                    chr_location_com[belong_loci]=[float(comp_this_l)]
            i_line+=1
    return chr_location_com

# chr_location_comp_cal=load_compartment_AB_data(chr_name,resolution)


def load_compartment_AB_data_cscore(chr_name,resolution):
    # dir_cpg=prefix+"/ygli/gam_paper_data/F123_cpg/GSM1027571_DNA_CpG_methcounts_E14_serum_LIF.bedGraph"
    dir_com=prefix+"/ygli/gam_paper_data/mes_hic/comp_AB_result_from_cscore/comp_result_"
    chr_location_com={}
    with open(dir_com+chr_name+"/"+chr_name+"_10kb_cscore.bedgraph",'r') as f:
        f.readline()
        all_lines=f.readlines()
        i_line=0
        for each_line in all_lines:
            pos_begin=each_line.replace("\n","").split()[1]
            comp_this_l=each_line.replace("\n","").split()[3]
            # belong_loci=int(int(pos_begin)/resolution)*resolution
            belong_loci_all=[pos_begin]
            for belong_loci in belong_loci_all:
                belong_loci=int(belong_loci)
                if belong_loci in chr_location_com:
                    chr_location_com[belong_loci].append(float(comp_this_l))
                else:
                    chr_location_com[belong_loci]=[float(comp_this_l)]
            i_line+=1
    return chr_location_com

chr_location_comp_cal=load_compartment_AB_data_cscore(chr_name,resolution)

loci_location_comp=[]
for loci_location_i in sorted_var_location:
    if loci_location_i in chr_location_comp_cal:
        # percent=len(np.where(np.array(chr_location_comp_cal[loci_location_i])>0.5)[0])/len(chr_location_comp_cal[loci_location_i])
        # percent=len(np.where(np.array(chr_location_comp_cal[loci_location_i])>0)[0])/resolution
        percent=np.sign(chr_location_comp_cal[loci_location_i])
        loci_location_comp.append(percent[0])
        # loci_location_comp.append(np.mean(chr_location_comp_cal[loci_location_i]))
    else:
        loci_location_comp.append(0)

# plot the +- for P1 result
# %matplotlib notebook
# %matplotlib qt5
fig = plt.figure(figsize=(16,2))
x_P1=chr_location_comp_cal.keys()

P1=[np.sign(chr_location_comp_cal[x]) for x in chr_location_comp_cal.keys()]
# P1=[np.sign(reduced_corr[x][0]) for x in range(len(reduced_corr))]
plt.scatter(x_P1,P1,s=0.9)
comp_dir=prefix+"/ygli/gam_paper_data/mes_hic/comp_AB_result_from_cscore/comp_result_"
plt.savefig(comp_dir+chr_name+"/"+chr_name+"_10kb_AB_compartment.png")
plt.close()


# with compartment AB as color: plot the 3D fig for whole chr all loci 
print("plot compartmentAB_3D ...")
fig = plt.figure()
ax = plt.axes(projection='3d')
sm=ax.scatter3D(rc_all,rs_all,h_all,c=loci_location_comp[0:len(rc_all)],cmap = 'RdBu',alpha = 0.6)
plt.colorbar(sm)

ax.legend(loc='best')
ax.set_xlabel('x', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('y', fontdict={'size': 15, 'color': 'red'})
ax.set_zlabel('z', fontdict={'size': 15, 'color': 'red'})
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(-5, 5)

ax.view_init(10, 30)
fig.savefig(png_dir+chr_name+"_10kb_compartmentAB_3D_"+method+"_1.eps",dpi=600,format='eps')
ax.view_init(10, -30)
fig.savefig(png_dir+chr_name+"_10kb_compartmentAB_3D_"+method+"_2.eps",dpi=600,format='eps')
ax.view_init(80, 10)
fig.savefig(png_dir+chr_name+"_10kb_compartmentAB_3D_"+method+"_3.eps",dpi=600,format='eps')
plt.close()


## plot the relationship of spatial distance with comp
spatial_distance_all=[]
comp_all=[]
for i in range(len(rc_all)):
    spatial_distance_all.append(np.linalg.norm([rc_all[i],rs_all[i],h_all[i]]))
    comp_all.append(loci_location_comp[i])

comp_df_whole_genome=pd.DataFrame({"spatial_distance":np.array(spatial_distance_all),
                     "comp":np.array(comp_all)})

fig = plt.figure()
ax = plt.axes()
comp_A_spatial_location=[]
for i in range(len(comp_all)):
    if comp_all[i]==1:
        comp_A_spatial_location.append(spatial_distance_all[i])

comp_B_spatial_location=[]
for i in range(len(comp_all)):
    if comp_all[i]==-1:
        comp_B_spatial_location.append(spatial_distance_all[i])

ax1 = sns.distplot(comp_A_spatial_location, color="r")
ax2 = sns.distplot(comp_B_spatial_location, color="b")

ax.set_xlabel('spatial distance')
plt.title("compartment AB with different spatial distance")
ax.set_ylabel('probability')
plt.legend(["Compartment A","Compartment B"])

fig.savefig(png_dir+chr_name+"_10kb_compartmentAB_distribution_"+method+".eps",dpi=600,format='eps')
plt.close()

