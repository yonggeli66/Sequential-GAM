# python calculate_radius.py

import seaborn as sns
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
import math
from scipy.stats import norm

# dir_cell_diameter = "/Users/april/Desktop/2017Autumn/microscope/gam/sequential_gam/第三批连续数据的chr分析/"
dir_cell_diameter = "/Users/april/Desktop/2017Autumn/microscope/gam/sequential_gam/第六批数据/"
# B16F信息.xlsx

# these are the area of the preknow data, eg. the area suspension cell and nucleus
area_1=[107.884,101.565,62.648,97.347,38.686,86.502,64.996,55.495,92.759,109.166,69.492,69.492,45.375,138.381,60.964,81.373,79.272,113.724,151.174,88.093,115.068,46.148,34.53,40.725,34.36,77,21.367,54.212,52.544,82.176,48.419,50.597,59.774,83.86,57.797,45.669,64.548,77.757,106.324,76.815,64.873,41.188,120.846,84.308,38.639,119.734,96.204,115.192,73.045,49.995,86.363,34.591,100.654,93.099,88.711,27.531,104.609,93.593,108.008,56.406,97.008,118.76,65.985,68.967]
# area_2_nu=[45.226,56.113,36.651,108.272,87.721,53.057,58.175,42.419,84.494,22.747,127.027,101.072,15.718,129.567,55.062,65.471,87.263,64.039,118.757,53.401,106.839,129.338,35.772,76.625,54.89,101.339,127.256,53.859,77.236,51.319,28.362,98.856,110.201,81.323,56.762,101.74,60.257,51.491,102.638,77.179,46.353,43.278,108.597,111.595,75.995,40.795,25.688,101.931,91.369,57.602,61.327,77.905,40.241,89.479,145.84,152.467,61.346,13.942,120.724,27.541,125.958,104.452,109.38,138.181,21.754,61.518,46.621,51.835,95.208,53.191,73.225,23.339,123.436,110.602,126.435,79.031,104.338,57.946,82.259,61.785,54.05,99.983,59.417,83.577,51.491,103.898,66.923,71.697,87.435,79.814,20.264,117.401,96.087,81.877,97.844,58.615,62.95,21.047,117.268,80.903,96.297,55.215,31.571,42.247,68.011,62.53,45.494,106.247,50.268,34.741,101.931,84.207,112.569,49.218,53.534,19.748,50.88]
area_2_nu=[53.057,58.175,42.419,84.494,22.747,127.027,101.072,15.718,129.567,55.062,65.471,87.263,64.039,118.757,53.401,106.839,129.338,35.772,76.625,54.89,101.339,127.256,53.859,77.236,51.319,28.362,98.856,110.201,81.323,56.762,101.74,60.257,51.491,102.638,77.179,46.353,43.278,108.597,111.595,75.995,40.795,25.688,101.931,91.369,57.602,61.327,77.905,40.241,89.479,145.84,152.467,61.346,13.942,120.724,27.541,125.958,104.452,109.38,138.181,21.754,61.518,46.621,51.835,95.208,53.191,73.225,23.339,123.436,110.602,126.435,79.031,104.338,57.946,82.259,61.785,54.05,99.983,59.417,83.577,51.491,103.898,66.923,71.697,87.435,79.814,20.264,117.401,96.087,81.877,97.844,58.615,62.95,21.047,117.268,80.903,96.297,55.215,31.571,42.247,68.011,62.53,45.494,106.247,50.268,34.741,101.931,84.207,112.569,49.218,53.534,19.748,50.88]
area_2_al=[127.352,97.634,110.64,144.98,68.909,206.651,158.674,162.952,204.378,125.499,88.963,157.356,154.721,181.077,177.295,137.78,243.168,101.664,213.889,131,141.007,203.041,122.844,134.858,129.968,82.068,194.16,242.48,171.699,225.73,184.133,108.463,84.226,216.085,167.001,87.76,85.124,155.16,218.645,202.888,85.219,145.056,176.455,206.135,150.595,181.86,157.967,79.566,210.069,181.478,295.002,157.719,97.767,168.262,143.49,175.214,208.618,227.869,169.35,92.114,86.194,114.46,93.107,183.598,124.029,164.843,134.571,367.12,180.943,248.324,138.926,144.98,148.628,150.958,114.002,101.129,209.84,144.541,159.342,95.055,227.984,177.868,160.584,212.457,131.668,76.51,189.29,187.437,221.949,129.72,161.252,150.423,141.485,219.848,190.416,227.64,135.736,78.363,72.576,133.559,148.38,202.086,198.61,98.092,104.242,180.37,187.571,158.884,105.006,99.085,64.497,162.322]

area_3_nu=[166.982,215.111,199.794,269.295,181.077,148.418,90.93,134.399,150.767,148.036,71.812,146.737,108.921,144.273,178.021,110.621,147.788,128.364,146.565,108.062,123.436,109.552,121.679,101.568,133.425,212.094,247.217,158.846,184.114,62.874,166.294,114.746,112.798,89.364,99.887,126.397,73.588,107.884,101.565,62.648,97.347,38.686,86.502,64.996,55.495,92.759,109.166,69.492,69.492,45.375,138.381,60.964,81.373,79.272,113.724,151.174,88.093,115.068,46.148,34.53,40.725,34.36,77,21.367,54.212,52.544,82.176,48.419,50.597,59.774,83.86,57.797,45.669,64.548,77.757,106.324,76.815,64.873,41.188,120.846,84.308,38.639,119.734,96.204,115.192,73.045,49.995,86.363,34.591,100.654,93.099,88.711,27.531,104.609,93.593,108.008,56.406,97.008,118.76,65.985,68.967]
area_3_al=[278.558,280.888,332.283,327.318,240.055,232.511,164.805,193.358,193.568,194.886,93.26,182.223,147.387,182.872,221.013,139.881,179.205,170.706,182.891,179.454,167.46,142.039,163.812,157.28,206.727,252.832,330.775,221.299,274.891,142.402,244.868,167.746,141.409,146.966,154.644,164.117,110.774,176.325,137.454,90.89,144.978,71.686,128.216,124.044,109.012,111.561,149.242,110.371,142.043,134.89,205.309,78.252,126.748,143.217,154.372,215.32,127.675,175.692,61.844,51.431,55.726,43.738,91.477,33.665,67.02,70.017,107.884,83.613,65.166,81.635,115.655,84.571,57.905,94.721,112.751,136.836,107.034,108.919,67.267,147.759,115.315,67.344,178.535,140.204,156.72,88.232,118.56,130.456,60.995,134.164,128.463,126.223,46.395,138.551,145.071,134.766,113.446,148.254,200.025,84.555,108.409,]

radius_2_nu=(np.array(area_2_nu)/math.pi)**(0.5)
radius_2_al=(np.array(area_2_al)/math.pi)**(0.5)

radius_3_nu=(np.array(area_3_nu)/math.pi)**(0.5)
radius_3_al=(np.array(area_3_al)/math.pi)**(0.5)

radius=radius_2_al

# plt.figure(figsize=(9,6))
# sns.distplot(radius, bins=9, hist=True, kde=True, norm_hist=False,
#             rug=True, vertical=False,label='distplot',
#             axlabel='radius/um',fit=norm)

# plt.legend()
# plt.grid(linestyle='--')
# plt.show()

name=["slice_nuclei","slice_cell","suspension_nuclei","suspension_cell"]
i=0
name_all=[]
for area in [area_2_nu,area_2_al,area_3_nu,area_3_al]:
    print(name[i])
    print(len(area))
    radius=(np.array(area)/math.pi)**(0.5)
    name_all.append(norm.fit(radius))
    print("the parameter of fit norm distribution:",norm.fit(radius))
    print("mean radius is ", np.mean(radius))
    print("calculated radius of sphere is", 4/math.pi*np.mean(radius))
    i+=1

# combine this two kinds of distribution, i.e. sliced and suspensioned

w=0.9
mean_fitted_nuclei=w*name_all[0][0]+(1-w)*name_all[2][0]
variance_fitted_nuclei=w**2*name_all[0][1]+(1-w)**2*name_all[2][1]
mean_fitted_cell=w*name_all[1][0]+(1-w)*name_all[3][0]
variance_fitted_cell=w**2*name_all[1][1]+(1-w)**2*name_all[3][1]

nuclei_distribution=[mean_fitted_nuclei,variance_fitted_nuclei]
cell_distribution=[mean_fitted_cell,variance_fitted_cell]

print("after weighted average, the mean and the variance of nuclei is:", nuclei_distribution)
print("after weighted average, the mean and the variance of cell is:", cell_distribution)
print("then the radius of the sphere of nuclei is:",4/math.pi*mean_fitted_nuclei)
print("then the radius of the sphere of cell is:",4/math.pi*mean_fitted_cell)

##################################################################################
######## transfer the raduis of the data we measured from cell to nucluis ########
##################################################################################

mu_1,sigma_1=cell_distribution
mu_2,sigma_2=nuclei_distribution
def transfer_r_from_cell_to_nucluis(r_cell):
    a=(2*sigma_2**2*(math.log(sigma_1/sigma_2)+(r_cell-mu_1)**2/(2*sigma_1**2)))**(0.5)
    if a<mu_2:
        r_nuclei=mu_2-a
    else:
        r_nuclei=mu_2+a
    return r_nuclei

# transfer_r_from_cell_to_nucluis(10)

Cell_name=[]
NP_name=[]
NP_mean_cell_radius=[]
NP_mean_nuclei_radius=[]
split_symbol="\t"
# split_symbol=","
# rewrite the radius in each cells
# output: Cell_name"\t"NP_name"\t"mean_cell_r"\t"mean_nuclear_r"\n"
with open(dir_cell_diameter+"nuclear_size_B16.csv",'r') as f:
    # all=readlines()
    # head=f.readline().replace("\n","").split(",")
    head=f.readline().replace("\n","").split(split_symbol)
    #head: ['Cell', 'experiments', 'NPs', 'nuclear size', '', '', '', '', 'mean_radius\n']
    Cell_name_index=[x for x in range(len(head)) if head[x] == 'Cell']
    NPs_index=[x for x in range(len(head)) if head[x] == 'NPs']
    mean_cell_radius_index=[x for x in range(len(head)) if head[x] == 'mean_radius']
    for each_line in f:
        NP_mean_cell_radius_each=each_line.replace("\n","").split(split_symbol)[mean_cell_radius_index[0]]
        if NP_mean_cell_radius_each!='' and NP_mean_cell_radius_each!=" ":
            Cell_name_each=each_line.replace("\n","").split(split_symbol)[Cell_name_index[0]]
            NP_name_each=each_line.replace("\n","").split(split_symbol)[NPs_index[0]]
            if Cell_name_each=="":
                print(each_line)
                print("here no name")
                break
                Cell_name.append("Cell_no_name")
            else:
                Cell_name.append(Cell_name_each)
            if NP_name_each=="":
                NP_name.append("NP_no_name")
            else:
                NP_name.append(NP_name_each)
            NP_mean_cell_radius.append(NP_mean_cell_radius_each)
            if NP_mean_cell_radius_each!='None' and NP_mean_cell_radius_each!=" ":
                NP_mean_nuclei_radius.append(transfer_r_from_cell_to_nucluis(float(NP_mean_cell_radius_each)))
                # print("cell:",float(NP_mean_cell_radius_each))
                # print("nuclei:",transfer_r_from_cell_to_nucluis(float(NP_mean_cell_radius_each)))
            else:
                NP_mean_nuclei_radius.append("None")

aa=[float(x) for x in NP_mean_cell_radius if x!='None']
# aa=[float(x) for x in NP_mean_cell_radius if x!=' ']
norm.fit(aa)
print("radiu of sphere of cell:",4/math.pi*norm.fit(aa)[0])

aa=[float(x) for x in NP_mean_nuclei_radius if x!='None']
norm.fit(aa)
print("radiu of sphere of nuclei:",4/math.pi*norm.fit(aa)[0])

with open(dir_cell_diameter+"nuclear_size_B16_calculate.csv",'w') as f:
    for each_index in range(len(NP_name)):
        if NP_mean_cell_radius[each_index]!='':
            f.write(Cell_name[each_index]+"\t"+NP_name[each_index]+"\t"+NP_mean_cell_radius[each_index]+"\t"+str(NP_mean_nuclei_radius[each_index])+"\n")


#################################################################
### after calculated the index of NP in calculate_NP_index.py ###
#################################################################
# after_index_NP_name=["B12F-NP409","B12F-NP410","B12F-NP413","B12F-NP414","B12F-NP415","B12F-NP416","B12F-NP417","B12F-NP418","B12F-NP419","B12F-NP420","B12F-NP421","B12F-NP422","B12F-NP423","B12F-NP424","B12F-NP425","B12F-NP426","B12F-NP427","B12F-NP428","B12F-NP429","B12F-NP430","B12F-NP432","B12F-NP433","B12F-NP434","B12F-NP435","B12F-NP436","B12F-NP437","B12F-NP438","B12F-NP439","B12F-NP440","B12F-NP441","B12F-NP442","B12F-NP443","B12F-NP445","B12F-NP446","B12F-NP447","B12F-NP448","B12F-NP449","B12F-NP450","B12F-NP451","B12F-NP452","B12F-NP453","B12F-NP456","B12F-NP457","B12F-NP458","B12F-NP459","B12F-NP460","B12F-NP461","B12F-NP462","B12F-NP464","B12F-NP467","B12F-NP468","B12F-NP469","B12F-NP470","B12F-NP471","B12F-NP472","B12F-NP473","B12F-NP474","B12F-NP475","B12F-NP476","B12F-NP477","B12F-NP478","B12F-NP479","B12F-NP480","B12F-NP481","B12F-NP482","B12F-NP483"]
# after_index_NP_distance=[0.8499853258319932,1.049985325831993,1.6499853258319932,1.8499853258319932,2.0499853258319933,2.249985325831993,2.4499853258319932,2.6499853258319934,2.849985325831993,3.0499853258319933,3.249985325831993,3.4499853258319932,3.6499853258319934,3.849985325831993,4.049985325831993,4.249985325831993,4.649985325831993,4.849985325831994,5.049985325831994,5.249985325831993,6.059566181,6.259566181,6.459566181,6.659566181,6.859566181,7.059566181,7.259566181,7.459566181,7.859566181,8.059566181,8.259566181,8.459566181,0.504946692,0.704946692,0.904946692,1.104946692,1.304946692,1.504946692,1.704946692,1.904946692,2.104946692,0.00083358,0.20083358,0.40083358,0.60083358,0.80083358,1.00083358,1.20083358,2.00083358,2.20083358,0.037359362,0.237359362,0.437359362,0.637359362,0.837359362,1.037359362,1.237359362,1.437359362,1.637359362,1.837359362,2.037359362,2.237359362,2.437359362,2.637359362,2.837359362,3.037359362]
after_index_NP_name=["B16F-NP620","B16F-NP618","B16F-NP616","B16F-NP614","B16F-NP613","B16F-NP600","B16F-NP602","B16F-NP604","B16F-NP606","B16F-NP608","B16F-NP609","B16F-NP611","B16F-NP623","B16F-NP625","B16F-NP627","B16F-NP629","B16F-NP631","B16F-NP633","B16F-NP635","B16F-NP637","B16F-NP639","B16F-NP640","B16F-NP641","B16F-NP642","B16F-NP621","B16F-NP619","B16F-NP617","B16F-NP615","B16F-NP613","B16F-NP601","B16F-NP603","B16F-NP605","B16F-NP607","B16F-NP610","B16F-NP612","B16F-NP624","B16F-NP626","B16F-NP628","B16F-NP630","B16F-NP632","B16F-NP634","B16F-NP636","B16F-NP638","B16F-NP658","B16F-NP660","B16F-NP657","B16F-NP656","B16F-NP655","B16F-NP654","B16F-NP653","B16F-NP652","B16F-NP651","B16F-NP650","B16F-NP648","B16F-NP646","B16F-NP645","B16F-NP644","B16F-NP643","B16F-NP672","B16F-NP671","B16F-NP670","B16F-NP669","B16F-NP668","B16F-NP667","B16F-NP666","B16F-NP665","B16F-NP664","B16F-NP663","B16F-NP662","B16F-NP661"]
after_index_NP_distance=[3.608297279,3.808297279,4.008297279,4.208297279,4.408297279,4.608297279,4.808297279,5.008297279,5.208297279,5.408297279,5.608297279,5.808297279,6.008297279,6.208297279,6.608297279,6.808297279,7.008297279,7.208297279,7.408297279,7.608297279,8.008297279,8.208297279,8.408297279,8.608297279,4.782533685,4.982533685,5.182533685,5.382533685,5.582533685,5.782533685,5.982533685,6.182533685,6.382533685,6.782533685,6.982533685,7.182533685,7.382533685,7.782533685,7.982533685,8.182533685,8.382533685,8.582533685,8.782533685,4.49316056,4.79316056,5.09316056,5.29316056,5.49316056,5.69316056,5.89316056,6.09316056,6.29316056,6.49316056,6.89316056,7.09316056,7.29316056,7.49316056,7.69316056,0.732993659,0.932993659,1.132993659,1.332993659,1.532993659,1.732993659,2.332993659,2.732993659,2.932993659,3.132993659,3.332993659,3.532993659]



after_index_NP_name_dir={}
for i in range(len(after_index_NP_name)):
    after_index_NP_name_dir[after_index_NP_name[i]]=after_index_NP_distance[i]

with open(dir_cell_diameter+"nuclear_size_B16_calculate_with_index.csv",'w') as f:
    f.write("Cell_name\tNP_name\tCell_radius\tNuclear_radius\tdistance\n")
    for each_index in range(len(NP_name)):
        if NP_mean_cell_radius[each_index]!='':            
            if NP_name[each_index] in after_index_NP_name_dir:
                f.write(Cell_name[each_index]+"\t"+NP_name[each_index]+"\t"+NP_mean_cell_radius[each_index]+"\t"+str(NP_mean_nuclei_radius[each_index])+"\t")
                f.write(str(after_index_NP_name_dir[NP_name[each_index]]))
                f.write("\n")









