# calculate_NP_index.py

import seaborn as sns
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import optimize
from scipy import stats as st
import sys

dir_fig="/Users/april/Desktop/2017Autumn/microscope/gam/sequential_gam/第六批数据/radius/"
# dir_fig="/Users/april/Desktop/2017Autumn/microscope/gam/sequential_gam/第三批连续数据的chr分析/radius/"
cell_here=sys.argv[1]
# cell_here=Cell7

###############################################
#######calculate the minimize likelyhood#######
###############################################
# input_data=[3.4601216993561756, 3.7272141887783925, 'None', 'None', 4.068139776500302,4.386564185328192]
# [important] this is from the middle result of calculate_radius.py
# all_input_data is the NP_mean_nuclei_radius
# Cell_name=['Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell7', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell8', 'Cell9', 'Cell9', 'Cell9', 'Cell9', 'Cell9', 'Cell9', 'Cell9', 'Cell9', 'Cell9', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell10', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11', 'Cell11']
# NP_name=['B12F-NP409', 'B12F-NP410', 'B12F-NP411', 'B12F-NP412', 'B12F-NP413', 'B12F-NP414', 'B12F-NP415', 'B12F-NP416', 'B12F-NP417', 'B12F-NP418', 'B12F-NP419', 'B12F-NP420', 'B12F-NP421', 'B12F-NP422', 'B12F-NP423', 'B12F-NP424', 'B12F-NP425', 'B12F-NP426', 'NP_no_name', 'B12F-NP427', 'B12F-NP428', 'B12F-NP429', 'B12F-NP430', 'B12F-NP432', 'B12F-NP433', 'B12F-NP434', 'B12F-NP435', 'B12F-NP436', 'B12F-NP437', 'B12F-NP438', 'B12F-NP439', 'NP_no_name', 'B12F-NP440', 'B12F-NP441', 'B12F-NP442', 'B12F-NP443', 'B12F-NP445', 'B12F-NP446', 'B12F-NP447', 'B12F-NP448', 'B12F-NP449', 'B12F-NP450', 'B12F-NP451', 'B12F-NP452', 'B12F-NP453', 'B12F-NP456', 'B12F-NP457', 'B12F-NP458', 'B12F-NP459', 'B12F-NP460', 'B12F-NP461', 'B12F-NP462', 'NP_no_name', 'NP_no_name', 'B12F-NP463', 'B12F-NP464', 'B12F-NP467', 'B12F-NP468', 'B12F-NP469', 'B12F-NP470', 'B12F-NP471', 'B12F-NP472', 'B12F-NP473', 'B12F-NP474', 'B12F-NP475', 'B12F-NP476', 'B12F-NP477', 'B12F-NP478', 'B12F-NP479', 'B12F-NP480', 'B12F-NP481', 'B12F-NP482', 'B12F-NP483']
# all_input_data=[2.9781395853696306, 3.253727432903659, 'None', 'None', 3.6056701430409994, 3.9352563555408784, 3.7889065686387413, 3.712094525881632, 3.4060551557376426, 3.0813204736412017, 3.883899532208773, 2.5068045298539987, 3.635299361049643, 3.9345148845553792, 4.2227022343687866, 4.255230439680869, 4.3059784795919755, 4.233811175284054, 'None', 4.007140046081906, 3.730664495303529, 3.952235765836386, 3.8321756985826334, 3.779658243373227, 3.776072944489558, 3.2570288234679925, 3.416159299838945, 3.4960678142038617, 3.158963720954577, 3.6059035920487004, 2.5638535706969456, 'None', 2.648318834323953, 2.885026268533311, 1.9508146049912853, 2.1612914481769097, 2.323053753599689, 2.0488655792755464, 2.3710003993816446, 3.661163244343678, 3.1937955408036816, 3.407700201356807, 3.497356961973009, 3.137346265673301, 3.2877921429137364, 0.2349717125613644, 0.5778267092629807, 1.1094356236164136, 1.659833286187777, 1.645403805370571, 2.3983594032984206, 2.4019276876396627, 'None', 'None', 'None', 3.545601278989679, 3.5242993161433223, 1.8971895182800074, 2.2704547135515183, 1.8695394942495622, 1.8538066108962488, 1.506213392940389, 2.037787289479426, 0.8756735970064651, 1.4054040241190164, 1.6851135271787445, 2.0801922110592836, 2.52011774832299, 2.5980731331218774, 3.0492687524910895, 3.080681942793446, 3.472502287087821, 2.96771948715669]

Cell_name=['Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell17', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell18', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell19', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20', 'Cell20']
NP_name=['B16F-NP620', 'B16F-NP618', 'B16F-NP616', 'B16F-NP614', 'B16F-NP613', 'B16F-NP600', 'B16F-NP602', 'B16F-NP604', 'B16F-NP606', 'B16F-NP608', 'B16F-NP609', 'B16F-NP611', 'B16F-NP623', 'B16F-NP625', 'NP_no_name', 'B16F-NP627', 'B16F-NP629', 'B16F-NP631', 'B16F-NP633', 'B16F-NP635', 'B16F-NP637', 'NP_no_name', 'B16F-NP639', 'B16F-NP640', 'B16F-NP641', 'B16F-NP642', 'B16F-NP621', 'B16F-NP619', 'B16F-NP617', 'B16F-NP615', 'B16F-NP613', 'B16F-NP601', 'B16F-NP603', 'B16F-NP605', 'B16F-NP607', 'NP_no_name', 'B16F-NP610', 'B16F-NP612', 'B16F-NP624', 'B16F-NP626', 'NP_no_name', 'B16F-NP628', 'B16F-NP630', 'B16F-NP632', 'B16F-NP634', 'B16F-NP636', 'B16F-NP638', 'B16F-NP658', 'B16F-NP660', 'NP_no_name', 'B16F-NP657', 'B16F-NP656', 'B16F-NP655', 'B16F-NP654', 'B16F-NP653', 'B16F-NP652', 'B16F-NP651', 'B16F-NP650', 'NP_no_name', 'B16F-NP648', 'B16F-NP646', 'B16F-NP645', 'B16F-NP644', 'B16F-NP643', 'B16F-NP672', 'B16F-NP671', 'B16F-NP670', 'B16F-NP669', 'B16F-NP668', 'B16F-NP667', 'NP_no_name', 'NP_no_name', 'B16F-NP666', 'NP_no_name', 'B16F-NP665', 'B16F-NP664', 'B16F-NP663', 'B16F-NP662', 'B16F-NP661']
all_input_data=[3.5075510424425747, 3.4771112386270273, 3.564899308411523, 3.186122478779855, 3.5461862694141164, 3.6380969588168397, 2.964166892416335, 3.8348248722034013, 3.7304324938077515, 3.229432254048545, 3.375374005824706, 3.2223539226723314, 3.34855299526693, 3.216100639584885, 'None', 3.6514971029504153, 3.19922509678647, 3.213150733749379, 3.09822761470771, 3.029273821866985, 2.96712739816121, 'None', 2.278666699174271, 2.019679944089721, 2.1456927756123836, 1.9846534519189225, 3.8553076877121066, 4.083502191161806, 3.96236142259721, 'None', 3.5461862694141164, 3.8938715968680033, 2.148550617678009, 3.1845877461207843, 2.8581177859419866, 'None', 2.6951035683246403, 1.7124192571994037, 2.2265334646558608, 2.0963901751311553, 'None', 2.51774046556683, 1.5492749650511408, 1.9156613539708043, 1.8057694067126433, 1.299922017195775, 0.9837286925444757, 3.353494997356173, 'None', 'None', 4.115729159708576, 4.019897706288727, 3.841272659097226, 3.8365522895925244, 3.884128882095469, 3.6497497341388385, 3.4974741545618717, 4.278714070346387, 'None', 3.645089381598788, 3.439529174200149, 3.4093451753962567, 3.401707192451183, 3.1304933058201936, 2.5813207605228197, 3.0211087259867346, 3.2603300025117417, 3.126003045049468, 3.3803128115525194, 2.748280033889766, 'None', 'None', 3.9714522210218024, 'None', 3.846796613691462, 3.691296544532363, 3.7124428933589884, 4.085843400976571, 4.242001311805865]

mean_sphere=6.145
sigma_sphere=1.165
h=0.2
sigma_slice=0.5/3
# sigma_slice=0.5/1.5

input_data=[all_input_data[x] for x in range(len(Cell_name)) if Cell_name[x]==cell_here]
input_NP_name=[NP_name[x] for x in range(len(Cell_name)) if Cell_name[x]==cell_here]

xdata = range(0,int(mean_sphere*2/h)+1)
x = np.array(xdata[:len(input_data)])
y = input_data
index_not_none = [x for x in range(len(input_data)) if input_data[x]!="None" ]
x_data=[x[i] for i in index_not_none]
y_data=[y[i] for i in index_not_none]
NP_name_data=[input_NP_name[i] for i in index_not_none]
##########
# try directly use optimization
##########

# theta is [x0,R]
def loss_likely(theta):
    return 1/(2*sigma_slice**2)*np.sum([(y_data[index]-((theta[1]**2-(theta[1]-(x_data[index]+theta[0])*h)**2)**(0.5)))**2 for index in range(len(x_data))])-math.log(1/((2*math.pi)**(0.5)*sigma_sphere)*math.exp(-(theta[1]-mean_sphere)**2/(2*sigma_sphere**2)))

# def loss_likely(x_data_here,x0,R):
#     return 1/(2*sigma_slice**2)*np.sum([(y_data[index]-((R**2-(R-(x_data[index]-x0)*h)**2)**(0.5)))**2 for index in range(len(x_data))])-math.log(1/((2*math.pi)**(0.5)*sigma_sphere)*math.exp(-(R-mean_sphere)**2/(2*sigma_sphere**2)))

range_r_sphere=1.5
cons=({'type': 'ineq',
       'fun': lambda theta: theta[0]},
      {'type': 'ineq',
       'fun': lambda theta: (mean_sphere+range_r_sphere*sigma_sphere)/h*2-(np.max(x_data)+theta[0])},
      {'type': 'ineq',
       'fun': lambda theta: mean_sphere-range_r_sphere*sigma_sphere-theta[1]},
      {'type': 'ineq',
       'fun': lambda theta: theta[1]-mean_sphere+range_r_sphere*sigma_sphere})

min_loss=10000
theta_est=0
for begin in range(0,int((mean_sphere+range_r_sphere*sigma_sphere)/h*2-np.max(x_data)),10):
# for begin in range(0,20,5):
    first_guess=[begin,mean_sphere-1.4*sigma_sphere]
    theta_est_ = optimize.minimize(fun=loss_likely, x0=first_guess,constraints=cons,method='COBYLA')
    print(begin,theta_est_)
    if theta_est_['fun']<min_loss:
        min_loss=theta_est_['fun']
        theta_est=theta_est_

print(theta_est)


x0_est,R_est=theta_est['x']

R=R_est
a=(x_data+x0_est)*h
print("output NP name and fitted result with sphere radius is:",R)
for i in range(len(NP_name_data)):
    print(NP_name_data[i])
for i in range(len(x_data)):
    print(a[i])

def return_R(h):
    return (R**2-(R-h)**2)**(0.5)

hh=np.linspace(0, 2*R, num=50, endpoint=True)
rr=return_R(hh)
plt.figure(figsize=(9,6))
plt.plot(hh,rr,'b-')
plt.plot((x_data+x0_est)*h,y_data,'r')
plt.xlabel("distance from buttom or top/um ")
plt.ylabel("radius/um")
plt.title("raduis with R=%f um,h=%f um, x0=%f um"%(R,h,x0_est))
plt.legend(["fitted value","NP radius"])
plt.grid(linestyle='--')
# plt.show()
plt.savefig(dir_fig+cell_here+".png")
plt.close()


