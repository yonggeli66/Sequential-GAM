# build_3D_structure_snpstrain_with_parameter_withchangedir.py

from __future__ import division
from __future__ import print_function
# from mpl_toolkits import mplot3d
import tensorflow as tf
# import tensorflow.compat.v1 as tf 
import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
import math
import random
import copy
from scipy.optimize import minimize
import sys
import pandas as pd
import os
import argparse
from scipy.stats import beta, expon, lognorm
from itertools import product 
from scipy.spatial.distance import pdist
import time
import warnings
warnings.filterwarnings('ignore')

# tf.disable_v2_behavior()34e
# matplotlib.use('Agg')


def sample_paramter_result(parameter_fit_bestresult, sample_chr_name, sampler_parameter, seed_sample=1995):
    sampler_parameter_all = ['i', 'p', 'b']
#     sampler_distribution = ['beta', 'expon', 'lognorm']
    sampler_distribution = ['beta', 'lognorm', 'lognorm']
    sampler_parameter_to_distribution = dict(zip(sampler_parameter_all, sampler_distribution))
    np.random.seed(seed_sample)

    distribution_here = sampler_parameter_to_distribution[sampler_parameter]
    str_para = parameter_fit_bestresult.loc[(parameter_fit_bestresult.para_name==sampler_parameter) & (parameter_fit_bestresult.chr_name==sample_chr_name)].Parameters.values[0]
    float_para = map(float,str_para.replace('(','').replace(')','').split(","))
    if distribution_here=='beta':    
        a, b, loc, scale = float_para
        sampler_result = beta.rvs(a, b, loc=loc, scale=scale, size=10000)
    elif distribution_here=='expon':
        loc, scale = float_para
        sampler_result = expon.rvs(loc=loc, scale=scale, size=10000)
    elif distribution_here=='lognorm':
        s, loc, scale = float_para
        sampler_result = lognorm.rvs(s, loc=loc, scale=scale, size=10000)
#     print(sampler_result)
    sampler_result_output = np.median(sampler_result)
    iter_round_upperbound = 300
    iter_round = 0
    p_cutoff = 10
    while sampler_parameter=='p' and sampler_result_output>p_cutoff and iter_round<iter_round_upperbound:
        seed_sample+=1
        np.random.seed(seed_sample)
        s, loc, scale = float_para
        sampler_result = lognorm.rvs(s, loc=loc, scale=scale, size=10000)
        sampler_result_output = np.median(sampler_result)
        iter_round+=1
    if sampler_parameter=='p' and sampler_result_output>p_cutoff: #if after iter, still big
        sampler_result_output = np.median(sampler_result[np.where(sampler_result<p_cutoff)])

    return sampler_result_output

# def sample_paramter_result(parameter_fit_bestresult, sample_chr_name, sampler_parameter, seed_sample=1995):
#     sampler_parameter_all = ['i', 'p', 'b']
#     sampler_parameter_dic = {'i': 0.45*1e3, 'p': 0.00459*1e3, 'b':3.0}
# # #     sampler_distribution = ['beta', 'expon', 'lognorm']
# #     sampler_distribution = ['beta', 'lognorm', 'lognorm']
# #     sampler_parameter_to_distribution = dict(zip(sampler_parameter_all, sampler_distribution))
# #     np.random.seed(seed_sample)

# #     distribution_here = sampler_parameter_to_distribution[sampler_parameter]
# #     str_para = parameter_fit_bestresult.loc[(parameter_fit_bestresult.para_name==sampler_parameter) & (parameter_fit_bestresult.chr_name==sample_chr_name)].Parameters.values[0]
# #     float_para = map(float,str_para.replace('(','').replace(')','').split(","))
# #     if distribution_here=='beta':    
# #         a, b, loc, scale = float_para
# #         sampler_result = beta.rvs(a, b, loc=loc, scale=scale, size=10000)
# #     elif distribution_here=='expon':
# #         loc, scale = float_para
# #         sampler_result = expon.rvs(loc=loc, scale=scale, size=10000)
# #     elif distribution_here=='lognorm':
# #         s, loc, scale = float_para
# #         sampler_result = lognorm.rvs(s, loc=loc, scale=scale, size=10000)
# # #     print(sampler_result)
#     return sampler_parameter_dic[sampler_parameter]


#################################################################################
#calculate point_r_set_kb_linear -> point_r_set_linear and point_r_set_kb_np -> point_r_set_np 
# and compare point_r_set_linear and point_r_set_np to get max
#################################################################################
# def the functions to calculate distances
def get_spatial_distance(x, sample_i, sample_p, sample_b):
    i, p, b = sample_i, sample_p, sample_b
    intercept=i
    para=p
    scale=1000 #for nm
    inflexion=intercept+para*(5*resolution)**(1/b)
    if x<5*resolution:
        y=x/(5*resolution)*inflexion
    else:
        y=intercept+para*(x)**(1/b)
    return y/scale


def get_spatial_distance_list(x,sample_i, sample_p, sample_b):
    i, p, b = sample_i, sample_p, sample_b
    intercept=i
    para=p
    scale=1000 #for nm
    inflexion=intercept+para*(5*resolution)**(1/b)
    y=np.zeros_like(x)
    for x_i in range(len(x)):
        if x[x_i]<5*resolution:
            y[x_i]=x[x_i]/(5*resolution)*inflexion
        else:
            y[x_i]=intercept+para*(x[x_i])**(1/b)
    return y/scale

def get_spatial_distance_array(x, sample_i, sample_p, sample_b):
    i, p, b = sample_i, sample_p, sample_b
    intercept=i
    para=p
    scale=1000 #for nm
    inflexion=intercept+para*(5*resolution)**(1/b)
    y=np.zeros_like(x)
    for x_i in range(len(x)):
        for x_j in range(len(x[0])):
            if x[x_i][x_j]<5*resolution:
                y[x_i][x_j]=(x[x_i][x_j]/(5*resolution)*inflexion)
            else:
                y[x_i][x_j]=(intercept+para*((x[x_i][x_j])**(1/b)))
    return y/scale


# #################################################################################
# #################################################################################
# build the result for whole chr
# #################################################################################
# #################################################################################
# begin_5_r,begin_5_phi is the 5 overlap point from end of last use_begin.
def run_phi_with_LAD_with_overlap(steps, lr, use_point_num, use_begin, begin_5_r, begin_5_phi, sorted_var_location, sorted_h, sorted_nearest_lad_distance_all, parameters, with_overlap=True):
    random_seed, limit_LAD, limit_res, limit_power_LAD, limit_power_res,constant_LAD, constant_res = parameters
    tf.set_random_seed(random_seed)#tf<2.0
    # tf.random.set_seed(random_seed)#tf>2.0
    
    location_matrix=np.array(sorted_var_location[use_begin:use_begin+use_point_num],dtype='float32')
    h_matrix=np.array(sorted_h[use_begin:use_begin+use_point_num],dtype='float32')
    sorted_nearest_lad_distance=np.array(sorted_nearest_lad_distance_all[use_begin:use_begin+use_point_num])
    # r_matrix=np.array(sorted_r[use_begin:use_begin+use_point_num],dtype='float32')

    location_matrix=location_matrix.reshape([1,use_point_num])
    h_matrix=h_matrix.reshape([1,use_point_num])
    # r_matrix=r_matrix.reshape([1,use_point_num])

    # cos_phi_sq=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=1.0, dtype=tf.float32, seed=1234))
    # sin_phi_sq=tf.ones([1,use_point_num])-cos_phi_sq

    #tf<2.0
    phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=random_seed)) 
    #tf>2.0
    # phi=tf.Variable(tf.random.uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=random_seed))

    
    #     phi=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=2*math.pi, dtype=tf.float32, seed=2020))
    phi_place=tf.placeholder(tf.float32,[1,use_point_num])
    phi_assign_opt=tf.assign(phi,phi_place)

    #tf<2.0
    r_matrix=tf.Variable(tf.random_uniform([1,use_point_num], minval=0.0, maxval=R-0.5, dtype=tf.float32, seed=random_seed))
    #tf>2.0
    # r_matrix=tf.Variable(tf.random.uniform([1,use_point_num], minval=0.0, maxval=R-0.5, dtype=tf.float32, seed=random_seed))
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
    dis=get_spatial_distance_array(np.abs(L_r-L_rT), sample_i, sample_p, sample_b)
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

    # r_z =(h**2+r**2)**2
    #calculate abs( (R-r)**2 - ( min LAD distance )**2  )
    r_sphere_matrix=tf.math.sqrt(tf.math.square(r_matrix)+tf.math.square(h_matrix))
    R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_sphere_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    # R_minus_matrix=tf.fill([use_point_num,use_point_num],R)-tf.convert_to_tensor(tf.tile(r_matrix,(use_point_num,1)))#[use_point_num,use_point_num]
    R_minus_matrix_square=tf.math.square(R_minus_matrix)

    # LAD_distance_array=get_spatial_distance(np.abs(sorted_nearest_lad_distance))
    LAD_distance_array=get_spatial_distance_list(np.abs(sorted_nearest_lad_distance),sample_i, sample_p, sample_b)
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

    if with_overlap == True:  # run results for the loci from beginning, like 0, 10kb, ...
        with tf.Session() as sess:
            sess.run(init)
            phi_now=sess.run(phi)
            phi_now[0][0:len(begin_5_phi)]=begin_5_phi
            r_now=sess.run(r_matrix)
            r_now[0][0:len(begin_5_r)]=begin_5_r
            loss_now = np.inf
            loss_before = [np.inf] *10
            # _,_=sess.run([parameter_A_assign_opt,parameter_C_assign_opt],feed_dict={parameter_A_place:a_for_parameter,parameter_C_place:c_for_parameter})
            for i in range(steps):
                if loss_now>=(0.1/250)*use_point_num:
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
                    loss_before.append(loss_now)
                    loss_before = loss_before[1:]
                    if np.std(loss_before)<loss_now*1e-7: #early stop with patience
                        break
                else:
                    break
            loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now=sess.run([cost_func,Res_matrix,rc,rs,r_matrix])
                    # print("Step %d,loss_: %f, phi_now:%s"%(i,loss_now,phi_now))
            M_rc_res_square_now,M_rs_res_square_now,limit_r_by_LAD_square_now=sess.run([M_rc_res_square,M_rs_res_square,limit_r_by_LAD_square])
            Res_matrix_original_now,limit_r_by_LAD_square_original_now=sess.run([Res_matrix_original,limit_r_by_LAD_square_original])
    
    else: # run results with restrain of the loci location has been calculated (eg 56 loci locations in last optimizing)
        with tf.Session() as sess:
            sess.run(init)
            phi_now=sess.run(phi)
            r_now=sess.run(r_matrix)
            loss_now = np.inf
            loss_before = [np.inf] *10
            # _,_=sess.run([parameter_A_assign_opt,parameter_C_assign_opt],feed_dict={parameter_A_place:a_for_parameter,parameter_C_place:c_for_parameter})
            for i in range(steps):
                if loss_now>=(0.1/250)*use_point_num:
            #         phi_now%(2*math.pi)
            #         _=sess.run(phi_assign_opt,feed_dice={phi_place:phi_now})
                    for i_r in range(len(r_now[0])):
                        if r_now[0][i_r]**2 + h_matrix[0][i_r]**2 > R**2:
                            r_now[0][i_r]=(R**2-h_matrix[0][i_r]**2)**(0.5)
                    _,_,_,loss_now,phi_now,r_now=sess.run([phi_assign_opt,r_assign_opt,opt_op,cost_func,phi,r_matrix],feed_dict={phi_place:phi_now%(2*math.pi),r_place:r_now})
                    if i%10000==0:
                        print("Step %d,loss_: %f"%(i,loss_now))
                    loss_before.append(loss_now)
                    loss_before = loss_before[1:]
                    if np.std(loss_before)<loss_now*1e-7:
                        break
                else:
                    break
            loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now=sess.run([cost_func,Res_matrix,rc,rs,r_matrix])
                    # print("Step %d,loss_: %f, phi_now:%s"%(i,loss_now,phi_now))
            M_rc_res_square_now,M_rs_res_square_now,limit_r_by_LAD_square_now=sess.run([M_rc_res_square,M_rs_res_square,limit_r_by_LAD_square])
            Res_matrix_original_now,limit_r_by_LAD_square_original_now=sess.run([Res_matrix_original,limit_r_by_LAD_square_original])

    return loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square,dis_square,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now
    # return loss_now,Res_matrix_now,rc_now,rs_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix


def plot_3d_fig(png_save_dir,chr_name, rc_all, rs_all, h_all, sorted_var_location, title_para):
    random_seed, limit_LAD, limit_res, limit_power_LAD, limit_power_res,constant_LAD, constant_res = title_para
    print("plot loci_3D ...")
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # ax.scatter3D(r_matrix*np.cos(phi_now),r_matrix*np.sin(phi_now),h_matrix,c='r', label='After calculation')
    ax.scatter3D(rc_all,rs_all,h_all,c=sorted_var_location, label='Sequential-GAM')
    ax.legend(loc='best')
    ax.set_xlabel('x', fontdict={'size': 15, 'color': 'red'})
    ax.set_ylabel('y', fontdict={'size': 15, 'color': 'red'})
    ax.set_zlabel('z', fontdict={'size': 15, 'color': 'red'})
    ax.set_xlim3d(-5, 5)
    ax.set_ylim3d(-5, 5)
    ax.set_zlim3d(-5, 5)
#     ax.title.set_text("limit_LAD %s, limit_res %s"%(limit_LAD, limit_res))
    # ax.title.set_text("limit_LAD %s, limit_res %s, limit_power_LAD %.3f, \nlimit_power_res %.3f, constant_LAD %s, constant_res %s"%(limit_LAD, limit_res, limit_power_LAD, limit_power_res,constant_LAD, constant_res ))
    ax.title.set_text(chr_name)
#     plt.show()
    ax.view_init(10, 30)
    fig.savefig(png_save_dir+chr_name+"_10kb_loci_3D_1.eps",dpi=600,format='eps')
    ax.view_init(10, -30)
    fig.savefig(png_save_dir+chr_name+"_10kb_loci_3D_2.eps",dpi=600,format='eps')
    ax.view_init(80, 10)
    fig.savefig(png_save_dir+chr_name+"_10kb_loci_3D_3.eps",dpi=600,format='eps')
    plt.close()
    


if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(description='parameter for sequential-GAM')
    parser.add_argument('--chr_name', type=str, default = "chr1")
    parser.add_argument('--cell_name', type=str, default= "6")
    # parser.add_argument('--method', type=str, default="overlap")
    parser.add_argument('--my_step', type=int, default=180000)
    parser.add_argument('--my_lr', type=float, default=1e-4)
    parser.add_argument('--strain_plot', type=str, default="C") #C means CAST, and S means 129S1
    parser.add_argument('--limit_LAD', type=float, default=1e-3)
    parser.add_argument('--limit_res', type=float, default=1e-3)
    parser.add_argument('--constant_LAD', type=float, default=1e-8)
    parser.add_argument('--constant_res', type=float, default=1e-3)
    parser.add_argument('--limit_power_LAD', type=int, default=3)
    parser.add_argument('--limit_power_res', type=int, default=1)
    parser.add_argument('--save_file_prefix', type=str, default="newmm10_depart_parameter_all_chr_10kb_")
    # parser.add_argument('--location_save_dir', type=str, default="./")
    parser.add_argument('--input_z_file_prefix', type=str, default="/home/ygli/gam_paper_data/gam_seq_mapped_414/z_location/gam_seq_mapped_414_")
    parser.add_argument('--input_LAD_file_prefix', type=str, default="/home/ygli/gam_paper_data/lad_data/GSE17051_cLAD_regions_mm10_10000/")
    parser.add_argument('--powerlaw_paramter_file', type=str, default="/data01/ygli/cailongseqFISH/DNAseqFISH+/fit_parameter_each_chr_sampler.csv")
    parser.add_argument('--run_file_prefix', type=str, default="./")
    parser.add_argument('--R', type=float, default=4.3975)
    parser.add_argument('--resolution', type=float, default=1e+4)
    parser.add_argument('--cuda', type=str, default="None")
    parser.add_argument('--cLAD_dropout', type=float, default=0)
    parser.add_argument('--power_law_i', type=str, default="450") #or "sampler"
    parser.add_argument('--power_law_p', type=str, default="4.59") #or "sampler"
    parser.add_argument('--power_law_b', type=str, default="3.0") #or "sampler"
    args = parser.parse_args()
    # args = parser.parse_args(args=['--cell_name','38', '--save_file_prefix', 'newcell_mm10_depart_parameter_all_chr_10kb_'])

    CUDA_VISIBLE_DEVICES = int(args.cuda) if args.cuda !="None" else None
    chr_name = args.chr_name
    cell_name = args.cell_name
    strain_plot = args.strain_plot
    # method = args.method
    steps = args.my_step
    lr = args.my_lr
    limit_LAD = args.limit_LAD
    constant_LAD = args.constant_LAD
    limit_res = args.limit_res
    constant_res = args.constant_res
    limit_power_LAD = 1/float(args.limit_power_LAD)
    limit_power_res = 1/float(args.limit_power_res)
    save_file_prefix = args.save_file_prefix
    R = args.R 
    resolution = args.resolution # consist with input mapped data

    prefix = args.run_file_prefix
    input_z_file_prefix = args.input_z_file_prefix
    input_LAD_file_prefix = args.input_LAD_file_prefix
    powerlaw_paramter_file = args.powerlaw_paramter_file
    cLAD_dropout = args.cLAD_dropout
    seed_sample = int(cell_name)

    power_law_i = args.power_law_i
    power_law_p = args.power_law_p
    power_law_b = args.power_law_b

    location_save_dir = prefix + "processed_data/set_1cell_"+str(cell_name)+"/"+save_file_prefix+strain_plot+"/"
    if not os.path.exists(location_save_dir):
        os.makedirs(location_save_dir)

    # z_axis_dir = prefix + "processed_data/set_1cell_"+str(cell_name)+"/all_chr_10kb_"+strain_plot+"/"


    # load NP data
    if strain_plot == "C":
        np_info = pd.read_csv(input_z_file_prefix+"Cell"+str(cell_name)+"_zlocation_1.csv")
    elif strain_plot == "S":
        np_info = pd.read_csv(input_z_file_prefix+"Cell"+str(cell_name)+"_zlocation_0.csv")

    np_info = np_info[np_info.chr==chr_name]
    np_info.loc[:,"zaxis"]-=R
    # np_info.head()

    # load cLAD data
    if cLAD_dropout==0:
        LAD_dir = os.path.join(input_LAD_file_prefix,chr_name+"_shortest_distance.csv")
    else:
        LAD_dir = os.path.join(input_LAD_file_prefix,chr_name+"_shortest_distance_dpot"+str(cLAD_dropout)+"_seed"+str(cell_name)+".csv")
    LAD_info = pd.read_csv(LAD_dir)
    # LAD_info.head()


    # merge info in np_info and LAD_info
    np_info["shortest_distance"] = [LAD_info[LAD_info.head_loci == i.start_resolution].shortest_distance.values[0] for _,i in np_info.iterrows()]

    # sample paramters of power-law for cell, extimated 3d distance = i +p*(x)**(1/b), x is 2d distance
    parameter_fit_bestresult = pd.read_csv(powerlaw_paramter_file, sep='\t')

    if power_law_i=="sampler":
        sample_i = sample_paramter_result(parameter_fit_bestresult, chr_name, 'i', seed_sample)
    else:
        sample_i = float(power_law_i)
    if power_law_p=="sampler":
        sample_p = sample_paramter_result(parameter_fit_bestresult, chr_name, 'p', seed_sample)
    else:
        sample_p = float(power_law_p)
    if power_law_b=="sampler":
        sample_b = sample_paramter_result(parameter_fit_bestresult, chr_name, 'b', seed_sample)
    else:
        sample_b = float(power_law_b)    
    
    print("cell%s, (sample_i, sample_p, sample_b)", (cell_name, sample_i, sample_p, sample_b))



    # overlap with 56 points
    # run all the loci for my data (for method with overlap)
    sorted_var_location = list(np_info.start_resolution)
    sorted_h = list(np_info.zaxis)
    sorted_nearest_lad_distance_all = list(np_info.shortest_distance)

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

    circle_num=int(len(sorted_var_location)/calculate_point_num)+1
    begin_list = np.linspace(0, circle_num*calculate_point_num, circle_num+1)[:-1]
    if len(sorted_var_location)-begin_list[-2] < use_point_num: #if the num overlaped points(use_point_num-calculate_point_num) is large than the left points 
        begin_list = begin_list[:-1]

    
    random_seed = int(cell_name)
    parameters = [random_seed, limit_LAD, limit_res, limit_power_LAD, limit_power_res,constant_LAD, constant_res]
    print("run with paramters as: random_seed %s, limit_LAD %s, limit_res %s, limit_power_LAD %s, limit_power_res %s,constant_LAD %s, constant_res %s."%(random_seed, limit_LAD, limit_res, limit_power_LAD, limit_power_res, constant_LAD, constant_res))
    # sorted_var_location, sorted_h, sorted_nearest_lad_distance_all


    for use_begin in begin_list:
        use_begin=int(use_begin)
        if use_begin != 0:
            if len(sorted_var_location)-use_begin < use_point_num: #when the last epoch do not have enough loci
                use_point_num=len(sorted_var_location)-use_begin
                lr *= 1e-1
            loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(steps, lr, use_point_num,use_begin, begin_5_r,begin_5_phi,sorted_var_location, sorted_h, sorted_nearest_lad_distance_all, parameters, with_overlap=True)
            loop_num=0
            lr_small = copy.deepcopy(lr)
            while np.isnan(loss_now) and loop_num<3:
                lr_small*=0.1
                loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(steps, lr_small, use_point_num,use_begin, begin_5_r,begin_5_phi,sorted_var_location, sorted_h, sorted_nearest_lad_distance_all, parameters, with_overlap=True)
                loop_num+=1
            if np.isnan(loss_now):
                loss_now_all.append(loss_now)
                Res_matrix_now_all.append(Res_matrix_now[overlap_num:])
                rc_now_all.append(np.nan*np.ones_like(rc_now[0][overlap_num:]))
                rs_now_all.append(np.nan*np.ones_like(rs_now[0][overlap_num:]))
                r_matrix_now_all.append(np.nan*np.ones_like(r_matrix_now[0][overlap_num:]))
                dis_square_all.append(dis_square[overlap_num:])
                L_r_all.append(L_r[overlap_num:])
                L_rT_all.append(L_rT[overlap_num:])
                h_matrix_all.append(h_matrix[0][overlap_num:])
                continue
                
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
            begin_5_r,begin_5_phi = [], []
            # loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD(7,1e-5,use_point_num,use_begin,dic_group_loci_location_h_C,dic_loci_location_head_nearest_lad_distance_point_C)
            loss_now,Res_matrix_now,rc_now,rs_now,phi_now,r_matrix_now,dis_square,L_r,L_rT,h_matrix,M_rc_res_square_now,M_rs_res_square_now,h_r_res_square_now,dis_square_now,limit_r_by_LAD_square_now,Res_matrix_original_now,limit_r_by_LAD_square_original_now=run_phi_with_LAD_with_overlap(steps, lr, use_point_num, use_begin, begin_5_r,begin_5_phi, sorted_var_location, sorted_h, sorted_nearest_lad_distance_all, parameters, with_overlap=False)
            begin_5_r=r_matrix_now[0][-overlap_num:]
            begin_5_phi=phi_now[0][-overlap_num:]

            loss_now_all.append(loss_now)
            Res_matrix_now_all.append(Res_matrix_now)
            rc_now_all.append(rc_now[0])
            rs_now_all.append(rs_now[0])
            r_matrix_now_all.append(r_matrix_now[0])
            dis_square_all.append(dis_square)
            L_r_all.append(L_r)
            L_rT_all.append(L_rT)
            h_matrix_all.append(h_matrix[0])
        print("calculated loci %s to %s / %s, loss_now is %s"%(use_begin, use_begin+use_point_num, len(sorted_var_location), loss_now))

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
    # save whole genome
    # save result for each [loci location] with [spatical coordinate:xyz]
    # the scatter plot
    # #################################################################################
    # strain_plot="C"
    dir_save=location_save_dir+chr_name+"_10kb_location_data"
    # print(dir_save)
    location_result_df = pd.DataFrame({
        'chr_name': [chr_name]*len(sorted_var_location),
        'loci_head': sorted_var_location,
        'rc_all': rc_all,
        'rs_all': rs_all,
        'h_all': h_all,
        })
    location_result_df.to_csv(dir_save+".csv", index=0)
    with open(dir_save+".txt","w") as f:
        for location_index in range(len(rc_all)):
            if ~np.isnan(rc_all[location_index]):
                f.write(str(sorted_var_location[location_index])+"\t"+str(rc_all[location_index])+"\t"+str(rs_all[location_index])+"\t"+str(h_all[location_index]))
                f.write("\n")
    end_time = time.time()

    print("saved in %s"%(dir_save+".csv"))
    print("time consumed: %s s"%(end_time-start_time))
    # plot_3d_fig(location_save_dir, chr_name, rc_all, rs_all, h_all, sorted_var_location, parameters)

