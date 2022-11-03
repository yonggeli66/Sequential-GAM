# python get_reads_zaxis.py

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import sys
import os
# from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.utils import shuffle
from sklearn.preprocessing import MinMaxScaler
import argparse
import random
random.seed(1995)

parser = argparse.ArgumentParser()
parser.add_argument("--experiment", type=str, default="gam_seq_mapped_219")
parser.add_argument("--np_zaxis_dir", type=str, default="/home/ygli/gam_paper_data/gam_seq_mapped_219/nuclear_size_B16_calculate_with_index.csv")
parser.add_argument("--cell_name", type=str, default="Cell18")
parser.add_argument("--dir_savefile", type=str, default="/home/ygli/gam_paper_data/gam_seq_mapped_219/z_location/")
args = parser.parse_args()
experiment = args.experiment
np_zaxis_dir = args.np_zaxis_dir
cell_name = args.cell_name
dir_savefile = args.dir_savefile

preprocess = MinMaxScaler()

# experiment = "gam_seq_mapped_219"
# np_zaxis_dir = "/home/ygli/gam_paper_data/gam_seq_mapped_219/nuclear_size_B16_calculate_with_index.csv"
# cell_list = ["Cell17", "Cell18", "Cell19", "Cell20"]

# experiment = "gam_seq_mapped_414"
# np_zaxis_dir = "/home/ygli/gam_paper_data/gam_seq_mapped_414/nuclear_size_B18_414_calculate_with_index.csv"
# cell_list = ["Cell25", "Cell26", "Cell27", "Cell28", "Cell29", "Cell38", "Cell39"]

np_zaxis = pd.read_csv(np_zaxis_dir, sep="\t", header=None)
np_zaxis.columns = ['Cell','NPs', 'mean_radius', 'NP_mean_nuclei_radius', 'zaxis']
cell_list = list(set(np_zaxis.Cell))
np_list = list(set(np_zaxis.NPs))
# np_zaxis.head()
# cell_list


def get_all_reads(cell_dir, experiment, cell_name, np_name):
    # read bowtie map results
    np_read_dir = os.path.join(cell_dir, np_name+"_mm10.sort.remove_dump.bam.sam.sort")
    all_reads = pd.read_csv(np_read_dir, sep='\t', header=None, usecols=[0,1,2,3,4,5,6,7])
    all_reads.columns = ['reads', 'flag', 'chr', 'start', 'mapq', 'CIGAR', 'rNext', 'posNext']
    all_reads['zaxis'] = [np_zaxis.loc[(np_zaxis.Cell==cell_name)&(np_zaxis.NPs==np_name)].zaxis.values[0]]*len(all_reads)
    all_reads['np_name'] = [np_name] *len(all_reads)
    # all_reads.head()

    # read dipc snp
    reads_dipc_dir = np_read_dir+".callsnp"
    reads_dipc_dir
    reads_dipc = pd.read_csv(reads_dipc_dir, sep='\t', header=None)
    reads_dipc.columns = ['reads', 'mate', 'rmate']
    reads_dipc['parents'] = [i.mate.split("!")[4] for _,i in reads_dipc.iterrows()]
    #. for unknown, 0 for paternal, and 1 for maternal
    # reads_dipc.head()

    # # get parents for each_reads in all_reads
    # all_reads_parents = []
    # for _,i in all_reads.iterrows():
    #     prt = reads_dipc[reads_dipc.reads==i.reads].parents
    #     if len(prt) != 0:
    #         all_reads_parents.append(prt.values[0])
    #     else:
    #         all_reads_parents.append('.')
    # all_reads['parents'] = all_reads_parents
    # # all_reads.head()

    # a much faster way to get parents(in reads_dipc) for each_reads in all_reads 
    df_merge = pd.merge(all_reads,reads_dipc.loc[:,["reads","parents"]],on='reads',how='left')
    df_merge.fillna(value='.', inplace=True)
    all_reads = df_merge
    
    return all_reads


# np_name = "B16F-NP638"
cell_dir = "/home/ygli/gam_paper_data/"+experiment+"/bowtie_result/"+cell_name
np_list_thiscell = list(set(np_zaxis[np_zaxis.Cell==cell_name].NPs))

all_reads_thiscell_list = []
for each_np_name in np_list_thiscell:
# for each_np_name in ['B16F-NP601','B16F-NP603']:
    print(each_np_name)
    all_reads_thisnp = get_all_reads(cell_dir, experiment, cell_name, each_np_name)
    all_reads_thiscell_list.append(all_reads_thisnp)

all_reads_thiscell = pd.concat(all_reads_thiscell_list, ignore_index=True)
all_reads_thiscell['parents_KNN'] = list(all_reads_thiscell.parents)
all_reads_thiscell['zaxis_minmax'] = [i[0] for i in preprocess.fit_transform(np.array(all_reads_thiscell['zaxis']).reshape(-1,1))]
all_reads_thiscell['start_minmax'] = [i[0] for i in preprocess.fit_transform(np.array(all_reads_thiscell['start']).reshape(-1,1))]
# all_reads_thiscell.head()


# #use KNN to classify the un
X = np.array([[i.zaxis_minmax, i.start_minmax] for _,i in all_reads_thiscell[all_reads_thiscell.parents!='.'].iterrows()])
y = np.array(map(int,list(all_reads_thiscell[all_reads_thiscell.parents!='.'].parents)))
indices = list(xrange(len(X)))
random.shuffle(indices)
X = X[indices]
y = y[indices]
# clf = LinearDiscriminantAnalysis()
clf = KNeighborsClassifier(n_neighbors=5)
clf.fit(X, y)
x_new = np.array([[i.zaxis, i.start_minmax] for _,i in all_reads_thiscell[all_reads_thiscell.parents=='.'].iterrows()])
# x_new = np.array([[i.zaxis_minmax, i.start_minmax] for _,i in all_reads_thiscell.iterrows()])
if len(x_new)>0:
    nonSNP_group=clf.predict(x_new)
else:
    nonSNP_group=[]
all_reads_thiscell.loc[all_reads_thiscell.parents=='.','parents_KNN']=nonSNP_group
parents_KNN_int = [ int(i.parents_KNN) for _,i in all_reads_thiscell.iterrows()]
all_reads_thiscell['parents_KNN'] = parents_KNN_int

# all_reads_thiscell.to_csv(experiment+"_"+cell_name+".csv", columns=['cell_name', 'np_name', 'zaxis', 'parents_KNN'])
all_reads_thiscell.to_csv(dir_savefile+experiment+"_"+cell_name+".csv",index=None)
print("%s, saved in %s"%(cell_name, dir_savefile+experiment+"_"+cell_name+".csv"))
print("before devide, 0 num %s, 1 num %s, unknown num %s"%(len(all_reads_thiscell[all_reads_thiscell.parents=='0']), \
    len(all_reads_thiscell[all_reads_thiscell.parents=='1']), \
    len(all_reads_thiscell[all_reads_thiscell.parents=='.'])))
# print("after devide, 0 num %s, 1 num %s, unknown num %s"%(len(all_reads_thiscell[all_reads_thiscell.parents_KNN=='0']), \
#     len(all_reads_thiscell[all_reads_thiscell.parents_KNN=='1']), \
#     len(all_reads_thiscell[all_reads_thiscell.parents_KNN=='.'])))
print("after devide, 0 num %s, 1 num %s, unknown num %s"%(len(all_reads_thiscell[all_reads_thiscell.parents_KNN==0]), \
    len(all_reads_thiscell[all_reads_thiscell.parents_KNN==1]), \
    0))


##########################################
# save the result for paternal and maternal
resolution=10000
all_reads_thiscell['start_resolution'] = [ int(i.start//resolution*resolution) for _,i in all_reads_thiscell.iterrows()]
all_reads_thiscell['end_resolution'] = [ int(int(i.start//resolution*resolution)+resolution) for _,i in all_reads_thiscell.iterrows()]
all_reads_thiscell_0 = all_reads_thiscell[all_reads_thiscell.parents_KNN==0]
all_reads_thiscell_0_zaxis = all_reads_thiscell_0.groupby(['chr', 'start_resolution', 'end_resolution'])[['zaxis']].agg(['mean'])
all_reads_thiscell_0_zaxis.columns = all_reads_thiscell_0_zaxis.columns.droplevel([0])
all_reads_thiscell_0_zaxis = all_reads_thiscell_0_zaxis.reset_index()
all_reads_thiscell_0_zaxis.rename(columns={'mean':'zaxis'},inplace=True)
all_reads_thiscell_0_zaxis.to_csv(dir_savefile+experiment+"_"+cell_name+"_zlocation_0.csv", index=None)

all_reads_thiscell_1 = all_reads_thiscell[all_reads_thiscell.parents_KNN==1]
all_reads_thiscell_1_zaxis = all_reads_thiscell_1.groupby(['chr', 'start_resolution', 'end_resolution'])[['zaxis']].agg(['mean'])
all_reads_thiscell_1_zaxis.columns = all_reads_thiscell_1_zaxis.columns.droplevel([0])
all_reads_thiscell_1_zaxis = all_reads_thiscell_1_zaxis.reset_index()
all_reads_thiscell_1_zaxis.rename(columns={'mean':'zaxis'},inplace=True)
all_reads_thiscell_1_zaxis.to_csv(dir_savefile+experiment+"_"+cell_name+"_zlocation_1.csv", index=None)

print("parents' zlocation saved in %s"%(dir_savefile+experiment+"_"+cell_name+"_zlocation_1.csv"))

