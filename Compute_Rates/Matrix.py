import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy import integrate, interpolate
from Params import *
from joblib import Parallel,delayed
from pwtools.signal import acorr
import sys

class Matrix:
    def __init__(self, pot, Input):
        self.n_jobs = Input.NProcs
        self.Wbb = np.zeros((pot.n_bound,pot.n_bound))
        self.Wbe = np.zeros((pot.n_bound,pot.n_cont))
        self.Wbc = np.zeros((pot.n_bound,))
        
        self.Temp      = Input.Temp
        self.beta      = 1/kb_Ha/Input.Temp 
        self.QCF_Model = Input.QCF_Model
        
        if(self.QCF_Model == 'Egl'):
            self.Egelstaff_initialize(pot.t)

    def compute_transition_rates(self, Input, pot, eigen):
        logging.info(" Computing the Bound - Bound transition rates using %d procs"%(self.n_jobs))
        self.Wbb = Parallel(n_jobs=self.n_jobs,backend='multiprocessing') (delayed(self.compute_bound_row)(pot,eigen,m) for m in range(pot.n_bound))
        self.Wbb = np.array(self.Wbb)
        logging.info(" Done Computing the Bound - Bound transition rates")
        
        logging.info(" Computing the Bound - Continuum transition rates using %d procs"%(self.n_jobs))
        self.Wbe = Parallel(n_jobs=self.n_jobs,backend='multiprocessing') (delayed(self.compute_cont_row)(pot,eigen,m) for m in range(pot.n_bound))
        self.Wbe = np.array(self.Wbe)
        self.Wbc = np.sum(self.Wbe,axis=1)
        logging.info(" Done Computing the Bound - Continuum transition rates")
        
        logging.info(" Applying Detailed balance for the bound-bound rates.")
        self.Apply_detailed_balance(pot)
        logging.info(" Done applying detailed balance.")
        
        logging.info(" Solving the Master Equation.")
        self.solve_ME(pot)
        
        logging.info(" Writing the transition rates output files.")
        self.Write_transition_rates(Input,pot,Input.Out_pref)
        

       
    def compute_ind_rate(self,pot,eigen,m,n,method=2):
        """ 
        Calculates the transition rate between state m and n 
        method = 1 : Calculate the autocorrelation of Vmn and integrate with cos(omega*t)
        method = 2 : Calculate the power spectrum directly.  
        """
        if(self.QCF_Model == 'Egl'):
            method = 1
            
        eig_n = eigen.get_eigenstate(pot,n)
        eig_m = eigen.get_eigenstate(pot,m)
        om_nm = (pot.Energy[n] - pot.Energy[m])/hbar

        Vnm = np.zeros(pot.Vf.shape[0])
        for i in range(len(Vnm)):
            f = eig_n * pot.Vf[i] * eig_m
            Vnm[i] = integrate.simps(f,pot.z)

        if(method == 1):
            Vnm_corr = acorr(Vnm,norm=False)/len(Vnm)
            wnm = 2 * integrate.trapz(Vnm_corr*np.cos(om_nm*pot.t),pot.t)
        elif(method == 2):
            psd = 0
            dt = pot.t[1] - pot.t[0]
            psd = np.sum(Vnm*np.exp(-1j*om_nm*pot.t))
            wnm = (abs(psd)*dt)**2/(pot.t[-1] - pot.t[0])
            
        if(self.QCF_Model != 'Egl'):
            wnm = wnm * self.compute_QCF(om_nm)
        else:
            wnm = wnm * self.compute_QCF_Egl(om_nm,wnm,Vnm_corr,pot.t)
            
        wnm = wnm/(hbar**2)
        
        return wnm

    def compute_bound_row(self,pot,eigen,m):
        Wn = np.zeros((pot.n_bound,))
        for n in range(m+1,pot.n_bound):
            Wn[n] = self.compute_ind_rate(pot,eigen,m,n,method=2)
        return Wn
         
    def compute_cont_row(self,pot,eigen,m):
        Wn = np.zeros((pot.n_cont,))
        for c in range(pot.n_bound, pot.n_bound+pot.n_cont):
            Wn[c-pot.n_bound] = self.compute_ind_rate(pot,eigen,m,c,method=2)
        return Wn
    
    def compute_QCF(self,om):
        if(self.QCF_Model == 'None'):
            return 1.0
        elif(self.QCF_Model == 'Std'):
            return 2.0/(1+np.exp(-self.beta*hbar*om))
        elif(self.QCF_Model == 'Harm'):
            return self.beta*hbar*om/(1-np.exp(-self.beta*hbar*om))
        elif(self.QCF_Model == 'Scho'):
            return np.exp(self.beta*hbar*om/2)   
        else:
            logging.error(" QCF_Model invalid.")
            sys.exit()
            
    def Egelstaff_initialize(self,t):
        alpha = self.beta*hbar/2
        self.Egl_tmax = ((t[-1])**2 - alpha**2)**0.5
        
        self.Egl_tau  = (t**2 + alpha**2)**0.5
        for i in range(1,len(t)-1):
            if(self.Egl_tau[-i] > self.Egl_tmax):
                self.Egl_tau[-i] = self.Egl_tmax
            else:
                break 
            
    def compute_QCF_Egl(self,om,w,corr,t):
        QCF = 2*np.exp(self.beta*hbar*om/2)/w
        f = interpolate.interp1d(t,corr)
        corr_shift = f(self.Egl_tau)
        for i in range(1,len(t)-1):
            if(self.Egl_tau[-i] >= self.Egl_tmax):
                corr_shift[-i] = 0
            else:
                break
        g = np.cos(om*t) * corr_shift
        QCF *= integrate.trapz(g,t)
        return QCF
    
    def Apply_detailed_balance(self,pot):
        for m in range(pot.n_bound):
            for n in range(m+1,pot.n_bound):
                om = (pot.E_bound[m] - pot.E_bound[n])
                self.Wbb[n][m] = self.Wbb[m][n] * np.exp(-self.beta*om)
                
    def solve_ME(self,pot):
        self.P0 = np.exp(-self.beta*pot.E_bound)
        self.P0 = self.P0/np.sum(self.P0)
        
        self.W = np.zeros((pot.n_bound,pot.n_bound))
        
        for n in range(pot.n_bound):
            for k in range(pot.n_bound):
                self.W[n][n] += self.Wbb[n][k]
            self.W[n][n] += self.Wbc[n]
            for m in range(pot.n_bound):
                self.W[n][m] -=  self.Wbb[m][n]
                
        self.tau_des = np.sum(np.linalg.inv(self.W)@self.P0)
        self.k_des   = 1/self.tau_des

        logging.info("\n   Overall desorption rate constant [1/s] = %6.8E"%(self.k_des*s2au))        
                
    def Write_transition_rates(self,Input,pot,pref):
        
        out_bb = Input.Out_Dir + pref + ".Wbb"
        out_bc = Input.Out_Dir + pref + ".Wbc"
        out_be = Input.Out_Dir + pref + ".Wbe"
        
        file_bb = open(out_bb,"w")
        file_bc = open(out_bc,"w")
        file_be = open(out_be,"w")
        
        logging.info(" Converting the transition rates to 1/s.")
        self.Wbb *= s2au
        self.Wbe *= s2au
        self.Wbc *= s2au
        
        # Bound - Bound 
        for m in range(pot.n_bound):
            file_bb.write("%4.6E "%(pot.E_bound[m]/eV2Ha))
        file_bb.write("\n")

        for m in range(pot.n_bound):
            for n in range(pot.n_bound):
                file_bb.write("%4.6E "%(self.Wbb[m][n]))
            file_bb.write("\n")

        # Bound - Continuum
        for e in range(pot.n_cont):
            file_be.write("%4.6E "%(pot.E_cont[e]/eV2Ha))
        file_be.write("\n")
        
        for m in range(pot.n_bound):
            file_bc.write("%02d  %4.6E \n"%(m,self.Wbc[m]))
            for e in range(pot.n_cont):
                file_be.write("%4.6E "%(self.Wbe[m][e]))
            file_be.write("\n")
        
        file_bb.close()
        file_bc.close()
        file_be.close()
        
    #_____________ Testing functions _________________#
    def Test_Corr(self,pot,eigen,m,n):
        eig_n = eigen.get_eigenstate(pot,n)
        eig_m = eigen.get_eigenstate(pot,m)
        om_nm = (pot.Energy[n] - pot.Energy[m])/hbar

        Vnm = np.zeros(pot.Vf.shape[0])
        for i in range(len(Vnm)):
            f = eig_n * pot.Vf[i] * eig_m
            Vnm[i] = integrate.simps(f,pot.z)

        
        Vnm_corr = acorr(Vnm,norm=False)/len(Vnm)
        plt.figure(10)
        plt.plot(pot.t/ps_au,Vnm_corr,'-',label="(%d,%d)"%(m,n))
        plt.legend()
        plt.xlabel('t [ps]')
        
        

