#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 14:00:27 2020

@author: ckondur
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# plt.close('all')

if(len(sys.argv) == 1):
    dt = 0.2
    T  = 300
    Pot_suff = "./Merged_Pot.dat"

elif(len(sys.argv) == 4):
    dt = float(sys.argv[1])
    T  = int(sys.argv[2])
    Pot_suff = str(sys.argv[3])
else:
    print(" Wrong number of inputs.")
    sys.exit()

Output_suff="./" + str(dt) + "fs_V_z_" + str(T) + "K"
write_fluc = False

Data  = np.loadtxt( Pot_suff,usecols=(1,2) )
num_z = len( np.unique(Data[:,0]) )
num_t = len(Data[:,0])//num_z

print("Total Number of time steps = ",num_t)

V      = np.zeros((num_t,num_z))
V_avg  = np.zeros((num_z))
V_fluc = np.zeros_like(V)
z      = np.zeros_like(V_avg)

V = np.reshape(Data[:,1],(num_t,num_z))
z = Data[:num_z,0]

for i in range(num_t):
    
    V[i] = V[i] - V[i][0]
    
V_avg = np.mean(V,axis=0)

for i in range(num_t):
    V_fluc[i] = V[i] - V_avg

V_rms = (np.mean(V_fluc**2,axis=0))**0.5

#### Writing the output files
f_avg  = open(Output_suff+".avg","w")
f_rms  = open(Output_suff+".rms","w")

for i in range(num_z):
    f_avg.write("%4.6f  %6.8f\n"%(z[i],V_avg[i]))
    f_rms.write("%4.6f  %6.8f\n"%(z[i],V_rms[i]))

if(write_fluc):
    f_fluc = open(Output_suff+".fluc","w")
    for i in range(num_t):
        for j in range(num_z):
            f_fluc.write("%6.8f "%(V_fluc[i][j]))
        f_fluc.write("\n")
    f_fluc.close()

f_avg.close()
f_rms.close()

######### Plotting ###
#plt.plot(z,V_avg)    
#plt.show()
