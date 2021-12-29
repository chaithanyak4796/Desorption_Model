#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 18:16:46 2021

@author: ckondur
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
import sys
from mpmath import hyp1f1, hyperu
from scipy.special import factorial, gamma
from scipy.linalg import eigh
from scipy.optimize import curve_fit
from scipy import fftpack
from Params import  *
import logging

############### The system is solved in atomic units [Bo,Ha,.....] #####################

class Eigen:
    def __init__(self,Input):
        logging.info("Initializing the eigen class")
        
    def potential(self,x,D,alpha,r0):
        V = D*( np.exp(-2*alpha*(x-r0)) - 2*np.exp(-alpha*(x-r0))  ) 
        # V = 0.5*mass*(x - np.mean(x))**2
        return V

    def build_hamiltonian(self,pot):
        N = len(pot.z)
        H = np.zeros((N,N),dtype='float')
        dx = pot.z[1] - pot.z[0]
    
        for i in range(N):
            for j in range(N):
                pref = (hbar**2*(-1)**(i-j))/(2*pot.mass*dx**2)
                if(i == j):
                    H[i][j] = pref * (np.pi)**2/3 + pot.V0[i]
                else:
                    H[i][j] = pref * 2/(i-j)**2
        return H

    def get_eigenstate(self,pot,m):
        y = pot.evecs[:,m]
        return y
        
    def plot_eigenstates(self,pot,nstates=10):
        """Plot the first ``nstates`` eigenstates """
        x = pot.z
        V = pot.V0
        vecs = pot.evecs
        vals = pot.evals
        
        plt.figure(1)#,figsize=(10,8))
        conv_fac = Bo2Angs
        ax = plt.gca()
        out = ax.plot(x*conv_fac,V,'k')
        E_inf = V[-1]
        
        # vecs = -1*vecs
        for m in range(nstates):
            vec = vecs[:,m]
            c = np.max(vec)
            if(np.max(vec) < 1E-4):
                c = np.min(vec)
            vec = vec/c
            amp = vals[m+1] - vals[m]
            
            if(vals[m] <= E_inf):
                ax.plot(x*conv_fac, 0.3*amp*vec + vals[m],'-',linewidth=2, label=r"$n=%i$"%m)
            else:
                ax.plot(x*conv_fac, 0.3*amp*vec + vals[m],'-',linewidth=2,label=r"$n=%i$"%m)
            
        min_val = min(np.min(V),0)
        plt.xlabel(r"$z  [\AA]$")
        plt.ylabel(r"$\psi(z)$")
        plt.ylim(min_val, vals[m] + amp)   
        ax.axes.yaxis.set_ticks([])
        # plt.legend()
        return out

    def check_normalization(self,pot,m):
        """Checks the spatial normalization (For development purposes only) """
        state = pot.evecs[:,m]
        temp = state * state
        norm = sp.simps(temp,pot.z)
        
        print("<%d|%d> = %f"%(m,m,norm))

    def check_overlap(self,pot,m,n):
        eig_m = pot.evecs[:,m]
        eig_n = pot.evecs[:,n]
        over  = sp.simps(eig_m*eig_n,pot.z)
        print("<%d|%d> = %f"%(m,n,over))
    
    def normalize_eigenstates(self,pot):
        N = pot.evecs.shape[0]
        for i in range(N):
            temp = pot.evecs[:,i] * pot.evecs[:,i]
            norm = (sp.simps(temp,pot.z))**0.5
            pot.evecs[:,i] /= norm
            

    def check_morse_eigenstate(self,pot,m,corr=1):
        """ Checks the eigen states for the morse Oscillator """
        vals = pot.evals
        vecs = pot.evecs
        x = pot.z
        D = pot.D
        r0 = pot.r0
        alpha = pot.alpha
        mass = pot.mass
        E_inf = pot.E_inf
        
        E_num = vals[m]
        state_num = vecs[:,m] 
        # state_num /= (sp.simps(state_num*state_num,x))**0.5
        
        N = ((2*mass*D)**0.5/(hbar*alpha)) - 1/2 
        E_thr = -((hbar*alpha)**2/(2*mass)) * (N - m)**2 
        # print(k)  
        print("Energy (numerical)   [eV] = ",E_num/eV2Ha)
        
        # print(r0)
        z = (2*N+1)*np.exp(-alpha*(x-r0))
        self.z = z
    
        if (E_num < E_inf):
            print("Energy (theoretical) [eV] = ",E_thr/eV2Ha)
            state_thr = np.zeros_like(z)
            for i in range(len(z)):
                state_thr[i] = z[i]**(N-m) * np.exp(-z[i]/2) * hyp1f1(-m,1+2*N-2*m,z[i])
            norm = 2*(N-m)*gamma(2*N-m+1)/factorial(m)/(gamma(2*N-2*m+1))**2 * alpha
            norm = norm**(1/2)
            state_thr = state_thr * norm
        else: 
            print("E_inf                [eV] = ",E_inf/eV2Ha)  
            state_thr = np.zeros((len(z),),dtype='complex')
            k = (2*mass*E_num)**0.5/(hbar*alpha)
            print("k                     [-] = ",k)
            norm = abs(gamma(-N-1j*k))/np.pi * (k*np.sinh(2*np.pi*k))**0.5 
            for i in range(len(z)):
                state_thr[i] = z[i]**(-1j*k) * np.exp(-z[i]/2) * hyperu(-N-1j*k,1-2*1j*k,z[i])
            state_thr = np.real(state_thr*norm)
            
            # dkdeta = (vals[m+1]-vals[m]) /( (2*vals[m]/mass)**0.5 * hbar *alpha )
            # norm = (2*k**2*dkdeta/np.pi)**0.5 * abs(gamma(-N+1j*k)/gamma(1+2*1j*k))
            
            # print(dkdeta, norm)
            # temp = norm*state_thr
            # print("Test_norm = %f"%(sp.simps(abs(temp)**2,x)))

            
        self.state = state_thr
        state_thr = corr*state_thr
        norm_num = sp.simps(state_num*state_num,x)
        norm_thr = sp.simps(state_thr*state_thr,x)
        state_thr = state_thr/(sp.simps(state_thr*state_thr,x))**0.5
        state_num = state_num/(sp.simps(state_num*state_num,x))**0.5
        print("Norm_num = ",norm_num)
        print("Norm_thr = ",norm_thr)
        plt.close(3)
        plt.figure(3)
        plt.plot(x,state_num,'-',label='Numerical')
        plt.plot(x,state_thr,'--',label='Theoretical')
        plt.legend()

    def calc_morse_eigenvalues(self,pot):
        k = ((2*pot.mass*pot.D)**0.5/(hbar*pot.alpha))
        # print(" k =", k)
    
        E = np.zeros(int(k))
        for i in range(int(k)):
            E[i] = -((hbar*pot.alpha)**2/(2*pot.mass)) * (k - i - 0.5)**2 
        return E

    def save_H0(self,pot,Dir):
        z = pot.z

        fname_state = Dir + "Eigen.states"
        fname_en    = Dir + "Eigen.energy"

        fst = open(fname_state, "w")
        fen = open(fname_en, "w")

        for j in range(len(pot.z)):
            fst.write("%6.4E  "%(z[j]))
        fst.write("\n")
        for j in range(len(pot.z)):
            fst.write("%6.4E  "%(pot.V0[j]))
        fst.write("\n")
        
        for i in range(pot.n_bound + pot.n_cont):
            fen.write("%3d  %6.4E\n"%(i, pot.Energy[i]))
            state = self.get_eigenstate(pot,i)
            for j in range(len(pot.z)):
                fst.write("%6.4E  "%(state[j]))
            fst.write("\n")

            
        fst.close()
        fen.close()
        

