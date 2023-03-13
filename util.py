# util.py
from __future__ import division
from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import math
import random
import copy
import sys
import seaborn as sns; sns.set()
import pandas as pd
import os
from scipy.spatial.distance import pdist
import time
# sns.set_style("white")
# sns.set_context("paper") #"talk"
#read_group:infact_group(# 0 A1, 1 A2, 2 B1, 3 B2, 4 B3, 5 none)
# mapping_dic_odd={0:0,1:1,2:2,3:3,4:4,5:5}
# mapping_dic_even={0:0,1:1,2:2,3:3,4:4,5:5}

mapping_dic_odd={0:3,1:4,2:2,3:0,4:5,5:1}
mapping_dic_even={0:4,1:1,2:2,3:3,4:0,5:5}
R=4.3975
R_all = R+0.5
# cell_name_all=[1,2,3,4,5,6,7,8,9]

paired_color_list = [(0.6509803921568628, 0.807843137254902, 0.8901960784313725), (0.12156862745098039, 0.47058823529411764, 0.7058823529411765), (0.6980392156862745, 0.8745098039215686, 0.5411764705882353), (0.2, 0.6274509803921569, 0.17254901960784313), (0.984313725490196, 0.6039215686274509, 0.6), (0.8901960784313725, 0.10196078431372549, 0.10980392156862745), (0.9921568627450981, 0.7490196078431373, 0.43529411764705883), (1.0, 0.4980392156862745, 0.0), (0.792156862745098, 0.6980392156862745, 0.8392156862745098), (0.41568627450980394, 0.23921568627450981, 0.6039215686274509), (1.0, 1.0, 0.6), (0.6941176470588235, 0.34901960784313724, 0.1568627450980392)]

# modification_list=["test"]
# dic_modification_dir={"test":"/home/ygli/gam_paper_data/histone_m/test.bedgraph"}
modification_list=["CpG1","CpG2","WGBS_CpG1","WGBS_CpG2","WGBS_CpG3","LAD","RNA_seq","TSS","H3k27ac","H3K4me1","H3k4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3","early","late","super_enhancer","CTCF","RAD21"]
dic_modification_dir={"CpG1":"/home/ygli/gam_paper_data/mm10_cpg/cpg_calcualtion_result_1/merge_cpg.bedgraph",
                      "CpG2":"/home/ygli/gam_paper_data/mm10_cpg/cpg_calcualtion_result_2/merge_cpg.bedgraph",
                      "WGBS_CpG1":"/home/ygli/gam_paper_data/histone_m/GSM2425461_S_ESC_8_XX.bedGraph",
                      "WGBS_CpG2":"/home/ygli/gam_paper_data/histone_m/GSM2425462_S_ESC_12_XY_WGBS.bedGraph",
                      "WGBS_CpG3":"/home/ygli/gam_paper_data/histone_m/GSM2425463_W_S_ESC_31_XY.bedGraph",
                      "H3k27ac":"/home/ygli/gam_paper_data/histone_m/4DNFIR8QIFXO_H3K27ac.bedgraph",
                      "super_enhancer":"/home/ygli/gam_paper_data/histone_m/SEA00201_3col.bedgraph",
                      "H3k4me3":"/home/ygli/gam_paper_data/histone_m/4DNFIQFI5VUR_H3K4me3.bedgraph",
                      "H3K4me1":"/home/ygli/gam_paper_data/histone_m/sorted_ENCFF067VMS_H3K4me1.bedgraph",
                      "H3K9ac":"/home/ygli/gam_paper_data/histone_m/sorted_ENCFF960QSW_H3K9ac.bedgraph",
                      "H3K9me3":"/home/ygli/gam_paper_data/histone_m/sorted_ENCFF093JWL_H3K9me3.bedgraph",
                      "H3K27me3":"/home/ygli/gam_paper_data/histone_m/sorted_ENCFF423ZVR_H3K27me3.bedgraph",
                      "H3K36me3":"/home/ygli/gam_paper_data/histone_m/sorted_ENCFF759WKW_H3K36me3.bedgraph",
                      "early":"/home/ygli/gam_paper_data/histone_m/4DNFIBX82LLQ_early.bedGraph",
                      "late":"/home/ygli/gam_paper_data/histone_m/4DNFIDD5LIAX_late.bedGraph",
                      "RNA_seq":"/home/ygli/gam_paper_data/histone_m/4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "housekeeping_RNA_seq":"/home/ygli/gam_paper_data/histone_m/4DNFIC6NZQNU_transcript_housekeeping.bedgraph.deduplication",
                      "housekeeping_loc":"/home/ygli/gam_paper_data/histone_m/housekeeping_transcript_loc_sorted.bedgraph",
                      "compartment_A_RNA_seq":"/home/ygli/gam_paper_data/histone_m/compA_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "compartment_B_RNA_seq":"/home/ygli/gam_paper_data/histone_m/compB_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "LAD":"/home/ygli/gam_paper_data/histone_m/ES_non-allelic_LAD_coordinates.bed",
                      "TSS":"/home/ygli/gam_paper_data/TSS/MmRefseqTssTes.bedgraph",
                      "OSN_loc":"/home/ygli/gam_paper_data/histone_m/GSE103053_OSN_CHIPseq.bed",
                      "OSN_RNA_seq":"",
                      "subcompartment_A1":"/home/ygli/gam_paper_data/histone_m/subcomp0_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "subcompartment_A2":"/home/ygli/gam_paper_data/histone_m/subcomp1_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "subcompartment_B1":"/home/ygli/gam_paper_data/histone_m/subcomp2_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "subcompartment_B2":"/home/ygli/gam_paper_data/histone_m/subcomp3_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "subcompartment_B3":"/home/ygli/gam_paper_data/histone_m/subcomp4_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "subcompartment_none":"/home/ygli/gam_paper_data/histone_m/subcomp5_4DNFIC6NZQNU_transcript.bedgraph.deduplication",
                      "CTCF":"/home/ygli/gam_paper_data/histone_m/sorted_4DNFISPYLSSJ_CTCF.bedgraph",
                      "RAD21":"/home/ygli/gam_paper_data/histone_m/sorted_4DNFIIN8X1PC_RAD21.bedgraph",
                      "ATAC":"/home/ygli/gam_paper_data/histone_m/GSE119663_F123_ATACseq_replicated_peaks.bed",
                      "subcompartment_A1_loc":"/home/ygli/gam_paper_data/histone_m/subcomp_A1.bedgraph",
                      "subcompartment_A2_loc":"/home/ygli/gam_paper_data/histone_m/subcomp_A2.bedgraph",
                      "subcompartment_B1_loc":"/home/ygli/gam_paper_data/histone_m/subcomp_B1.bedgraph",
                      "subcompartment_B2_loc":"/home/ygli/gam_paper_data/histone_m/subcomp_B2.bedgraph",
                      "subcompartment_B3_loc":"/home/ygli/gam_paper_data/histone_m/subcomp_B3.bedgraph",
                      "test_feature":"/home/ygli/gam_paper_data/histone_m/test.bedgraph"
                      }

dir_com="/home/ygli/gam_paper_data/mes_hic_mm10/comp_AB_result_from_cscore/comp_result_"

cell_name_all=[17,18,19,20,25,26,27,28,38,39]
# cell_name_all=[17,18,19,20,25,26,27,28,29,38,39]
chr_index_list=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
chr_name_list=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"]
chr_reverse_ab=["chr2","chr3","chr6","chr7","chr8","chr11","chr13","chr15","chr17"]
dic_chr_size={"chr1":195471971,"chr2":182113224,"chr3":160039680,"chr4":156508116,"chr5":151834684,"chr6":149736546,"chr7":145441459,"chr8":129401213,"chr9":124595110,"chr10":130694993,"chr11":122082543,"chr12":120129022,"chr13":120421639,"chr14":124902244,"chr15":104043685,"chr16":98207768,"chr17":94987271,"chr18":90702639,"chr19":61431566,"chrX":171031299}
dic_chr_density = {'chr7': 14.7, 'chr6': 9.7, 'chr5': 10.2, 'chr4': 10.0, 'chr3': 8.6, 'chr2': 10.7, 'chr1': 8.4, 'chrY': 1.7, 'chr9': 11.7, 'chr8': 10.1, 'chr13': 8.9, 'chr12': 8.2, 'chr11': 15.9, 'chr10': 9.4, 'chr17': 13.1, 'chr16': 8.6, 'chr15': 9.4, 'chr14': 8.7, 'chr19': 14.1, 'chr18': 7.5, 'chrX': 7.9}
chr_tpm = {"chr1": 56032.31, "chr2": 55401.97, "chr3": 36505.20, "chr4": 52034.82, "chr5": 54309.34, "chr6": 46030.41, "chr7": 90891.35, "chr8": 45326.30, "chr9": 70674.76, "chr10": 52076.54, "chr11": 94751.69, "chr12": 31602.93, "chr13": 34898.72, "chr14": 31599.46, "chr15": 38827.28, "chr16": 23784.64, "chr17": 62046.45, "chr18": 21367.52, "chr19": 27640.01}
chr_logtpm = {key:np.log(chr_tpm[key]) for key in chr_tpm}
chr_logtpmpersize = {key:np.log(chr_tpm[key])/dic_chr_size[key]*1e6 for key in chr_tpm}


rnaseq = "/home/ygli/gam_paper_data/histone_m/4DNFIC6NZQNU_transcript.bedgraph.deduplication"
rnaseq_result = pd.read_csv(rnaseq, sep='\t', header=None)
rnaseq_result.columns = ['chr', 'head', 'tail', 'TPM', 'RPKM']
rnaseq_result = rnaseq_result.groupby(['chr'])[['TPM']].agg(['sum'])
rnaseq_result.columns = rnaseq_result.columns.droplevel([0])
rnaseq_result = rnaseq_result.reset_index()
rnaseq_result.head()
# rnaseq_result.loc[rnaseq_result['chr']=='chr1', 'sum'].values[0]


#read_group:infact_group(# 0 A1, 1 A2, 2 B1, 3 B2, 4 B3, 5 none)

def load_genome_location_coordinate(file_dir,chr_name,method,interpolation_flag=0):
    sorted_var_location=[]
    dic_location_index={}
    radial_distance=[]
    rc_all=[]
    rs_all=[]
    h_all=[]
    if interpolation_flag==1:
        dir_save=file_dir+chr_name+"_10kb_location_data_interpolation_"+method+".txt"
    else:
        dir_save=file_dir+chr_name+"_10kb_location_data_"+method+".txt"
    if "newmm10" in file_dir or "explore" in file_dir:
        if interpolation_flag==1:
            # print("use not interpolation data")
            dir_save=os.path.join(file_dir, chr_name+"_10kb_location_data_interpolation.txt")
        else:
            dir_save=os.path.join(file_dir, chr_name+"_10kb_location_data.txt")
            # print("use interpolation data")
    # print(dir_save)
    with open(dir_save,"r") as f:
        all_lines=f.readlines()
        i=0
        for each_line in all_lines:
            sorted_var_location_i,rc_all_i,rs_all_i,h_all_i=each_line.replace("\n","").split("\t")
            sorted_var_location.append(int(float(sorted_var_location_i)))
            dic_location_index[int(float(sorted_var_location_i))]=i
            rc_all.append(float(rc_all_i))
            rs_all.append(float(rs_all_i))
            h_all.append(float(h_all_i))
            radial_distance.append(np.linalg.norm([float(rc_all_i),float(rs_all_i),float(h_all_i)]))
            i+=1
            # 
    sorted_var_location=np.array(sorted_var_location)
    rc_all=np.array(rc_all)
    rs_all=np.array(rs_all)
    h_all=np.array(h_all)
    return dic_location_index,sorted_var_location,rc_all,rs_all,h_all,radial_distance


def load_RNAseq_data(chr_name,resolution):
    dir_RNAseq=dic_modification_dir["RNA_seq"]
    chr_location_data={}
    with open(dir_RNAseq,'r') as f:
        f.readline()
        all_lines=f.readlines()
        # i_line=0
        for each_line in all_lines:
            chr_name_this_line=each_line.replace("\n","").split()[0]
            if chr_name_this_line==chr_name:
                pos_begin=each_line.replace("\n","").split()[1]
                pos_end=each_line.replace("\n","").split()[2]
                data_this_l=each_line.replace("\n","").split()[3]
                # belong_loci=int(int(pos_begin)/resolution)*resolution
                belong_loci_begin=int(int(pos_begin)/resolution)*resolution
                belong_loci_end=int(int(pos_end)/resolution)*resolution
                for belong_loci in range(int(belong_loci_begin),int(belong_loci_end),int(resolution)):
                    belong_loci=int(belong_loci)
                    if belong_loci in chr_location_data:
                        chr_location_data[belong_loci].append(float(data_this_l))
                    else:
                        chr_location_data[belong_loci]=[float(data_this_l)]
    return chr_location_data


# resolution is the 
def load_subcompartment_data(chr_name,resolution,mapping_dic_even,mapping_dic_odd):
    if int(chr_name.replace("chr",""))%2==0:
        dir_subcom="/home/ygli/gam_paper_data/mes_hic_mm10/subcomp_eveninter_4DNFIASLBBXV_100kb.txt"
        # dir_subcom="/home/ygli/gam_paper_data/mes_hic/subcomp_eveninter_GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.txt"
        mapping_dic=mapping_dic_even
    else:
        dir_subcom="/home/ygli/gam_paper_data/mes_hic_mm10/subcomp_oddinter_4DNFIASLBBXV_100kb.txt"
        # dir_subcom="/home/ygli/gam_paper_data/mes_hic/subcomp_oddinter_GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.txt"
        mapping_dic=mapping_dic_odd
    chr_location_subcom={}
    # group_index_loci_head={0:[],1:[],2:[],3:[],4:[]}
    # group_index_loci_tail={0:[],1:[],2:[],3:[],4:[]}
    group_index_loci_head={0:[],1:[],2:[],3:[],4:[],5:[]}
    group_index_loci_tail={0:[],1:[],2:[],3:[],4:[],5:[]}
    with open(dir_subcom,'r') as f:
        # f.readline()
        all_lines=f.readlines()
        # i_line=0
        for each_line in all_lines:
            chr_name_this_line=each_line.replace("\n","").split()[0]
            if chr_name_this_line==chr_name:
                pos_begin=each_line.replace("\n","").split()[1]
                pos_end=each_line.replace("\n","").split()[2]
                comp_this_l=each_line.replace("\n","").split()[3]
                # belong_loci=int(int(pos_begin)/resolution)*resolution
                belong_loci_begin=int(float(pos_begin)/resolution)*resolution
                belong_loci_end=int(float(pos_end)/resolution)*resolution
                belong_loci_all=range(int(belong_loci_begin),int(belong_loci_end),int(resolution))
                for belong_loci in belong_loci_all:
                    belong_loci=int(belong_loci)
                    indeed_group=mapping_dic[int(float(comp_this_l))]
                    #changed here from original
                    chr_location_subcom[belong_loci]=[indeed_group]
                    group_index_loci_head[indeed_group].append(belong_loci)
                    group_index_loci_tail[indeed_group].append(belong_loci+int(resolution))
                # i_line+=1
    return chr_location_subcom,group_index_loci_head,group_index_loci_tail


def plot_part_distance_matrix():
    # #################################################################################
    # plot the distance matrix
    # #################################################################################
    show_distance_matrix=np.ones([len(inprint_gene_location_all),len(inprint_gene_location_all)])*2*R
    for loci_index_i in range(len(inprint_gene_location_all)):
        for loci_index_j in range(0,loci_index_i+1):
            if inprint_gene_location_all[loci_index_i] in dic_location_index and inprint_gene_location_all[loci_index_j] in dic_location_index:
                loci_1_index=dic_location_index[inprint_gene_location_all[loci_index_i]]
                loci_2_index=dic_location_index[inprint_gene_location_all[loci_index_j]]
                loci_1=np.array([rc_all[loci_1_index],rs_all[loci_1_index],h_all[loci_1_index]])
                loci_2=np.array([rc_all[loci_2_index],rs_all[loci_2_index],h_all[loci_2_index]])
                show_distance_matrix[loci_index_i,loci_index_j]=np.linalg.norm(loci_1-loci_2)
                show_distance_matrix[loci_index_j,loci_index_i]=np.linalg.norm(loci_1-loci_2)

    print("plot show_distance_matrix ...")
    f, ax = plt.subplots(figsize=(16, 16))
    df = pd.DataFrame(show_distance_matrix)
    sns.heatmap(df, annot=False, ax=ax,cmap = 'RdYlBu_r',cbar=False)
    # sns.heatmap(df, annot=False, ax=ax,cmap = 'RdBu',cbar=False)
    ax.set_xlabel('genome location')
    plt.title("genome distance matrix")
    ax.set_ylabel('genome location')
    print(png_dir+chr_name+"_10kb_"+str(inprint_gene)+"inprint_distance_matrix_"+method+"_notnorm.png")
    # plt.savefig(png_dir+chr_name+"_10kb_genome_distance_matrix_"+method+".eps",format='eps')
    plt.savefig(png_dir+chr_name+"_10kb_"+str(inprint_gene)+"inprint_distance_matrix_"+method+"_notnorm.png")
    plt.close()

def get_spatial_distance(x,resolution=1e+4):
    intercept=0.45
    para=0.00459
    inflexion=intercept+para*(5*resolution)**(1/3)
    if x<5*resolution:
        y=x/(5*resolution)*inflexion
    else:
        y=intercept+para*(x)**(1/3)
    return y



def sep_strlist_2_float(strlist,num_bool):
    # strlist = result.ground_TF[2]
    # num_bool = 'bool'
    # strlist = result.distance[2]
    # num_bool = 'num'
    if num_bool == 'bool':
        # strlist_t = list(strlist.replace(' ','').replace('[','').replace(']','').split(','))
        strlist_t = [True if i=='True' else False for i in list(strlist.replace(' ','').replace('[','').replace(']','').split(','))]
    elif num_bool == 'num':
        strlist_t = list(map(float,strlist.replace('[','').replace(']','').split(',')))
    return strlist_t



def read_all_cell(cell_name_all):
    chr_name_all_captured=[]
    radial_distance_all=[[] for i in range(len(chr_index_list))]
    radial_distance_all_without_blank=[]
    cell_name_all_each = []
    chr_name_all_each = []
    loci_location_all = []
    loci_strain_all = []
    for strain_plot in ['C','S']:
        for cell_name in cell_name_all:
    #         genome_location_dir=prefix+"ygli/gam_paper_data/single_cell_3D/set_1cell/set_1cell_"+str(cell_name)+"/"+file_prefix+strain_plot+"/"    
            genome_location_dir=prefix+"ygli/gam_paper_data/code/processed_data/set_1cell_"+str(cell_name)+"/"+file_prefix+strain_plot+"/"
    #         genome_location_dir=prefix+"ygli/gam_paper_data/single_cell_3D/set_1cell/set_1cell_"+str(cell_name)+"/"+file_prefix+strain_plot+"/"    
            for chr_index in chr_index_list:
                chr_name="chr"+str(chr_index)
    #             if os.path.exists(genome_location_dir+chr_name+"_10kb_location_data_"+method+".txt"):
                if os.path.exists(os.path.join(genome_location_dir, chr_name+"_10kb_location_data.txt")):
    #                 print("cell_name %s, chr_index %s"%(cell_name, chr_index))
                    _,sorted_var_location,_,_,_,radial_distance_this_chr=load_genome_location_coordinate(genome_location_dir,chr_name,method)
                    radial_distance_all[chr_index-1].extend(radial_distance_this_chr)
                    loci_location_all.extend(sorted_var_location)
                    loci_strain_all.extend([strain_plot]*len(radial_distance_this_chr))
                    cell_name_all_each.extend([cell_name]*len(radial_distance_this_chr))
                    chr_name_all_each.extend([chr_name]*len(radial_distance_this_chr))
    for chr_index in chr_index_list:
        if len(radial_distance_all[chr_index-1])!=0:
            chr_name_all_captured.append("chr"+str(chr_index))
            radial_distance_all_without_blank.append(radial_distance_all[chr_index-1])
#     print([len(cell_name_all_each),len(chr_name_all_each),len(np.array(radial_distance_all_without_blank).reshape(-1))])
    radial_df = pd.DataFrame({
        "cell_name": cell_name_all_each,
        "chr_name": chr_name_all_each,
        "loci_head": loci_location_all,
        "loci_strain": loci_strain_all,
        "radial_distance": np.array([j for i in radial_distance_all for j in i ])
    })
    return chr_name_all_captured, radial_distance_all_without_blank, radial_df



def get_TAD_location_by_TAD_headlist(R,resolution,prefix,chr_name,strain_plot,method,TAD_head_list,TAD_tail_list):
    # whole_TAD_distance_matrix_all=[]
    TAD_location_allcell=[]
    for cell_name in cell_name_all:
        genome_dir=prefix+"ygli/gam_paper_data/code/processed_data/set_1cell_"+str(cell_name)+"/"+file_prefix+strain_plot+"/"
        dic_location_index,sorted_var_location,rc_all,rs_all,h_all,_=load_genome_location_coordinate(genome_dir,chr_name,method) #without interpolation
        TAD_location=cal_TAD_location(resolution,dic_location_index,TAD_head_list,TAD_tail_list,rc_all,rs_all,h_all)
        # whole_TAD_distance_matrix_all.append(cal_TAD_distance_matrix(R,TAD_location))
        TAD_location_allcell.append(TAD_location)
    return TAD_location_allcell
