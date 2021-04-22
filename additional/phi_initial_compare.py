import numpy as np
import matplotlib.pyplot as plt
import collections

nz = 101
rho_eff_low = np.ones(nz)
x1 = 0.5            
nz1 = int(x1 * nz)
nz2 = nz
rho_eff_low[0] = 150
for i in range(nz1-1):
    rho_eff_low[i]= 150
rho_eff_low[nz1-1] = 131.25
rho_eff_low[nz1] = 112.5
rho_eff_low[nz1+1] = 93.75
rho_eff_low[nz1+2:nz2] = 75
counter=collections.Counter(rho_eff_low)
print(counter)


nz = 251
x1 = 0.5            
nz1 = int(x1 * nz)
nz2 = nz

rho_eff_high = np.ones(nz)#*75
rho_eff_high[0] = 150
for i in range(nz1-5):
    rho_eff_high[i]= 150
fill_list1 = np.linspace(150,112.5,6)
fill_list2 = np.linspace(112.5,75,6)

rho_eff_high[nz1-5:nz1] = fill_list1[:-1]
rho_eff_high[nz1] = 112.5
rho_eff_high[nz1+1:nz1+6] = fill_list2[1:]
rho_eff_high[nz1+6:] = 75
counter=collections.Counter(rho_eff_high)
print(counter)
#plt.plot(rho_eff_low, np.linspace(0,0.5,101), label = '101 nodes')
plt.plot(rho_eff_high, np.linspace(0,0.5,251), label ='251 nodes')
plt.show()