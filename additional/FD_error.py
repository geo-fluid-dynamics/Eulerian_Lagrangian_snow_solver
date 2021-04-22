import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm 
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
from matplotlib.patches import Rectangle
from model_geometry import node_distance
import os
os.chdir(r"C:\Users\Anna Lara\Documents\05_GIT\Snowmodel\FD_error")


all_T_noFD = np.loadtxt('all_T_vTc_eta_n1_noFDterm')
all_T_FD = np.loadtxt('all_T_vTc_eta_n1_FDterm')
all_t_passed_noFD = np.loadtxt('all_t_passed_vTc_eta_n1_noFDterm')
all_t_passed_FD = np.loadtxt('all_t_passed_vTc_eta_n1_FDterm')
all_coord_noFD = np.loadtxt('all_coord_vTc_eta_n1_noFDterm')
all_coord_FD = np.loadtxt('all_coord_vTc_eta_n1_FDterm')
all_phi_FD = np.loadtxt('all_phi_vTc_eta_n1_FDterm')
all_phi_noFD = np.loadtxt('all_phi_vTc_eta_n1_noFDterm')
all_phi_noFD = np.loadtxt('all_phi_vTc_eta_n1_noFDterm')


t3 = 24*2*3600  # total time [s]
t1 = 3600*24/2               # 1/3 time
t2 = 3600*24              # 2/3 time 
t4 = 3600 *24 *1.5            # 2/3 time 

length_noFD = len(all_t_passed_noFD)
t1_array_noFD = np.ones_like(length_noFD) * t1
t2_array_noFD = np.ones_like(length_noFD) * t2
t3_array_noFD = np.ones_like(length_noFD) * t3
t4_array_noFD = np.ones_like(length_noFD) * t4

t1_diff_noFD = np.absolute(t1_array_noFD - all_t_passed_noFD)
t2_diff_noFD = np.absolute(t2_array_noFD - all_t_passed_noFD)
t3_diff_noFD = np.absolute(t3_array_noFD - all_t_passed_noFD)
t4_diff_noFD = np.absolute(t4_array_noFD - all_t_passed_noFD)

t1_diff_list_noFD = list( t1_diff_noFD)
t2_diff_list_noFD = list( t2_diff_noFD)
t3_diff_list_noFD = list( t3_diff_noFD)
t4_diff_list_noFD = list( t4_diff_noFD)

# find indices that correspond to t1 t2 t3
t1_index_noFD = t1_diff_list_noFD.index(min(t1_diff_list_noFD[:-2]))
t2_index_noFD = t2_diff_list_noFD.index(min(t2_diff_list_noFD[:-2]))
t3_index_noFD = t3_diff_list_noFD.index(min(t3_diff_list_noFD[:-2]))
t4_index_noFD = t4_diff_list_noFD.index(min(t4_diff_list_noFD[:-2]))

length_FD = len(all_t_passed_FD)
t1_array_FD = np.ones_like(length_FD) * t1
t2_array_FD = np.ones_like(length_FD) * t2
t3_array_FD = np.ones_like(length_FD) * t3
t4_array_FD = np.ones_like(length_FD) * t4

t1_diff_FD = np.absolute(t1_array_FD - all_t_passed_FD)
t2_diff_FD = np.absolute(t2_array_FD - all_t_passed_FD)
t3_diff_FD = np.absolute(t3_array_FD - all_t_passed_FD)
t4_diff_FD = np.absolute(t4_array_FD - all_t_passed_FD)

t1_diff_list_FD = list( t1_diff_FD)
t2_diff_list_FD = list( t2_diff_FD)
t3_diff_list_FD = list( t3_diff_FD)
t4_diff_list_FD = list( t4_diff_FD)

# find indices that correspond to t1 t2 t3
t1_index_FD = t1_diff_list_FD.index(min(t1_diff_list_FD[:-2]))
t2_index_FD = t2_diff_list_FD.index(min(t2_diff_list_FD[:-2]))
t3_index_FD = t3_diff_list_FD.index(min(t3_diff_list_FD[:-2]))
t4_index_FD = t4_diff_list_FD.index(min(t4_diff_list_FD[:-2]))

coord_12h = abs(all_coord_noFD[t1_index_noFD,:])
coord_24h = abs(all_coord_noFD[t2_index_noFD,:])
coord_48h = abs(all_coord_noFD[t3_index_noFD,:])
coord_36h = abs(all_coord_noFD[t4_index_noFD,:])

dz_12h = node_distance(coord_12h, 101)
dz_24h = node_distance(coord_24h, 101)
dz_36h = node_distance(coord_36h, 101)
dz_48h = node_distance(coord_48h, 101)

T_diff_12h = abs(all_T_noFD[t1_index_noFD,:] - all_T_FD[t1_index_FD,:])
T_diff_24h = abs(all_T_noFD[t2_index_noFD,:] - all_T_FD[t2_index_FD,:])
T_diff_36h = abs(all_T_noFD[t4_index_noFD,:] - all_T_FD[t4_index_FD,:])
T_diff_48h = abs(all_T_noFD[t3_index_noFD,:] - all_T_FD[t3_index_FD,:])

T_diff_24h_dz = (T_diff_24h[1:] * dz_24h)
T_diff_36h_dz = (T_diff_36h[1:] * dz_36h)
T_diff_48h_dz = (T_diff_48h[1:] * dz_48h)


T_diff_24h_dz_sum = np.sum(T_diff_24h[1:] * dz_24h)
T_diff_36h_dz_sum = np.sum(T_diff_36h[1:] * dz_36h)
T_diff_48h_dz_sum = np.sum(T_diff_48h[1:] * dz_48h)


height_24h = np.sum(dz_24h)
height_36h = np.sum(dz_36h)
height_48h = np.sum(dz_48h)

T_diff_24h_l1 = T_diff_24h_dz_sum/height_24h
T_diff_36h_l1 = T_diff_36h_dz_sum/height_36h
T_diff_48h_l1 = T_diff_48h_dz_sum/height_48h


plt.plot(all_coord_noFD[t3_index_noFD],all_T_noFD[t3_index_noFD])
plt.plot(all_coord_FD[t3_index_FD],all_T_FD[t3_index_FD])
plt.show()
plt.close()

plt.plot(T_diff_12h, label= '12 h')
plt.plot(T_diff_24h, label= '24 h')
plt.legend()
plt.title('absolute error temperature')
plt.show()
plt.close()
plt.plot(T_diff_48h, label = 'T_diff, 48 h')
plt.plot(T_diff_36h, label = 'T_diff, 36 h')
#plt.plot(all_coord_FD[t3_index_FD], all_phi_FD[t3_index_FD], label = 'ice volume fraction FD')
#plt.plot(all_coord_noFD[t3_index_noFD], all_phi_noFD[t3_index_noFD], label = 'ice volume fraction noFD')

plt.legend()
plt.title('absolute error temperature')
plt.show()
plt.close()
