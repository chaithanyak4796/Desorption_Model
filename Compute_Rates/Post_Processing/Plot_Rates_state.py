#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 15:11:22 2021

@author: chaithanya
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
plt.close('all')

if(len(sys.argv) == 1):
    T = 800
else:
    T = int(sys.argv[1])

dt    = 0.2
t_max = 40
label = str(dt) + "fs_" + str(t_max) + "ps"

# Dir = "/media/chaithanya/Chaithanya_New/Surface_Chem/Desorption/Density_Matrix/Results/Oxygen/Edge_site/Model_3/"
# Dir = Dir + label + "/"

label = "Method_1"
Dir = "/media/chaithanya/Chaithanya_New/Surface_Chem/Desorption/Density_Matrix/Results/Oxygen/Bridge_site/Test_dz/"
Dir = Dir + label + "/"

# T = 900
ini_state = [0,5,10,20,25]
colors    = ['b','g','r','c','m','y','k']
lst       = '-'

pref     = "T_" + str(T) + "K"
fname_bb = Dir + pref + ".Wbb"
fname_bc = Dir + pref + ".Wbc"
fname_be = Dir + pref + ".Wbe"

data_bb = np.loadtxt(fname_bb)
E_bound = data_bb[0]
n_bound = len(E_bound)
Wbb = (data_bb[1:])

data_bc = np.loadtxt(fname_bc)
Wbc     = (data_bc[:,1])  

data_be = np.loadtxt(fname_be)
E_cont  = data_be[0]
Wbe     = data_be[1:]

for i in range(len(ini_state)):
    n  = ini_state[i]
    st = colors[i] + lst
    plt.figure(51)
    plt.semilogy(E_bound[n:],Wbb[n][n:],st,label='n = '+str(n))
    
    plt.figure(52)
    plt.semilogy(E_cont,Wbe[n],st,label='n = '+str(n))
    
plt.figure(51)
plt.xlabel('E_bound [eV]')
plt.ylabel(r'$W_{n\rightarrow n+1}[1/s]$')
plt.legend()

plt.figure(52)
plt.xlabel('E_cont [eV]')
plt.ylabel(r'$W_{n\rightarrow \epsilon}[1/s]$')
plt.legend()