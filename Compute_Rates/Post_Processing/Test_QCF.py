#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 13:32:22 2021

@author: chaithanya
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from Params import *

if(len(sys.argv) == 1):
    T = 300
else:
    T = int(sys.argv[1])

plt.close('all')

dt    = 0.4
t_max = 40
label = str(dt) + "fs_" + str(t_max) + "ps"

# label='Ads-no'
Dir = "/media/chaithanya/Chaithanya_New/Surface_Chem/Desorption/Density_Matrix/Results/Oxygen/Edge_site/Model_3/"
Dir = Dir + label + "/"

# T = 1000
pref = "T_" + str(T) + "K"
beta = 1/(kb_eV*T)

def QCF_Std(beta,dE):
    QCF = 2/(1+np.exp(-beta*dE))
    return QCF
def QCF_Harm(beta,dE):
    QCF = beta*dE/(1-np.exp(-beta*dE))
    return QCF
def QCF_Scho(beta,dE):
    QCF = np.exp(beta*dE/2)
    return QCF

fname_bb = Dir + pref + ".Wbb"
fname_bc = Dir + pref + ".Wbc"
fname_be = Dir + pref + ".Wbe"

data_bb = np.loadtxt(fname_bb)
E_bound = data_bb[0]
n_bound = len(E_bound)
Wbb = data_bb[1:]

data_be = np.loadtxt(fname_be)
E_cont  = data_be[0]
n_cont  = len(E_cont)
Wbe     = data_be[1:]
Wbc = np.sum(Wbe,axis=1)

Wbb_Std  = np.zeros_like(Wbb)
Wbb_Harm = np.zeros_like(Wbb)
Wbb_Scho = np.zeros_like(Wbb)
Wbb_Egl  = np.zeros_like(Wbb)
Wbb_hyb  = np.zeros_like(Wbb)

Wbe_Std  = np.zeros_like(Wbe)
Wbe_Harm = np.zeros_like(Wbe)
Wbe_Scho = np.zeros_like(Wbe)
Wbe_Egl  = np.zeros_like(Wbe)
Wbe_hyb  = np.zeros_like(Wbe)

bs_no = np.arange(n_bound)

for m in range(n_bound):
    for n in range(m+1,n_bound):
        dE = -(E_bound[m]-E_bound[n])
        Wbb_Std[m][n] = Wbb[m][n] * QCF_Std(beta, dE)
        Wbb_Std[n][m] = Wbb_Std[m][n] * np.exp(beta*(E_bound[n]-E_bound[m]))
        
        Wbb_Harm[m][n] = Wbb[m][n] * QCF_Harm(beta, dE)
        Wbb_Harm[n][m] = Wbb_Harm[m][n] * np.exp(beta*(E_bound[n]-E_bound[m]))
        
        Wbb_Scho[m][n] = Wbb[m][n] * QCF_Scho(beta, dE)
        Wbb_Scho[n][m] = Wbb_Scho[m][n] * np.exp(beta*(E_bound[n]-E_bound[m]))
        
        Wbb_hyb[m][n] = Wbb[m][n] * ( QCF_Scho(beta, dE) * QCF_Harm(beta, dE) )**0.5
        Wbb_hyb[n][m] = Wbb_hyb[m][n] * np.exp(beta*(E_bound[n]-E_bound[m]))
        
    for e in range(n_cont):
        dE = -(E_bound[m]-E_cont[e])
        Wbe_Std[m][e]  = Wbe[m][e] * QCF_Std(beta, dE)
        Wbe_Harm[m][e] = Wbe[m][e] * QCF_Harm(beta, dE)
        Wbe_Scho[m][e] = Wbe[m][e] * QCF_Scho(beta, dE)
        Wbe_hyb[m][n]  = Wbe[m][n] * ( QCF_Scho(beta, dE) * QCF_Harm(beta, dE) )**0.5
        
        
Wbc_Std = np.sum(Wbe_Std,axis=1)
Wbc_Harm = np.sum(Wbe_Harm,axis=1)
Wbc_Scho = np.sum(Wbe_Scho,axis=1)
Wbc_hyb  = np.sum(Wbe_hyb,axis=1)
        
plt.figure(1)
b_jump = 1
plt.semilogy(bs_no[:-b_jump],np.diagonal(Wbb,b_jump),label='No QCF')
plt.semilogy(bs_no[:-b_jump],np.diagonal(Wbb_Std,b_jump),label='QCF: Standard')
plt.semilogy(bs_no[:-b_jump],np.diagonal(Wbb_Harm,b_jump),label='QCF: Harmonic')
# plt.semilogy(bs_no[:-b_jump],np.diagonal(Wbb_Scho,b_jump),label='QCF: Schofield')
# plt.semilogy(bs_no[:-b_jump],np.diagonal(Wbb_hyb,b_jump),label='QCF: harmonic/Schofield')
plt.legend()
plt.xlabel('n')
plt.ylabel(r'$W_{n\rightarrow m} [s^{-1}]$')

plt.figure(2)
plt.semilogy(bs_no,Wbc,label='No QCF')
plt.semilogy(bs_no,Wbc_Std,label='QCF: Standard')
plt.semilogy(bs_no,Wbc_Harm,label='QCF: Harmonic')
# plt.semilogy(bs_no,Wbc_Scho,label='QCF: Schofield')
# plt.semilogy(bs_no,Wbc_hyb,label='QCF: harmonic/Schofield')
plt.legend()
plt.xlabel('n')
plt.ylabel(r'$W_{n\rightarrow cont} [s^{-1}]$')

# plt.show()

if(1):
    fw_name_bb = Dir + pref + "-Harm" + ".Wbb"
    fw_name_bc = Dir + pref + "-Harm" + ".Wbc"
    
    file_bb = open(fw_name_bb,"w")
    file_bc = open(fw_name_bc,"w")
    
    for m in range(n_bound):
        file_bb.write("%4.6E "%(E_bound[m]))
    file_bb.write("\n")

    for m in range(n_bound):
        for n in range(n_bound):
            file_bb.write("%4.6E "%(Wbb_Harm[m][n]))
        file_bb.write("\n")

    # Bound - Continuum
    for m in range(n_bound):
        file_bc.write("%02d  %4.6E \n"%(m,Wbc_Harm[m]))

