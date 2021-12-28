import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy import integrate, interpolate
from Params import *
from joblib import Parallel,delayed
from scipy.fftpack import fft, ifft
import sys

class Matrix:
    def __init__(self, pot, Input):
        self.n_jobs = Input.NProcs
        self.Wbb = np.zeros((pot.n_bound,pot.n_bound))
        self.Wbe = np.zeros((pot.n_bound,pot.n_cont))
        self.Wbc = np.zeros((pot.n_bound,))

        self.Wbb_corr = np.zeros_like(self.Wbb)
        self.Wbe_corr = np.zeros_like(self.Wbe)
        self.Wbc_corr = np.zeros_like(self.Wbc)
        
        self.Temp      = Input.Temp
        self.beta      = 1/kb_Ha/Input.Temp 
        self.QCF_Model = Input.QCF_Model
        
        self.use_Filon = True
        
        if(self.QCF_Model == 'Egl'):
            self.Egelstaff_initialize(pot.t)

    def compute_transition_rates(self, Input, pot, eigen):
        logging.info(" Computing the Bound - Bound transition rates using %d procs"%(self.n_jobs))
        temp = Parallel(n_jobs=self.n_jobs,backend='multiprocessing') (delayed(self.compute_bound_row)(pot,eigen,m) for m in range(pot.n_bound))
        self.temp = temp
        self.Wbb      = [item[0] for item in temp]
        self.Wbb_corr = [item[1] for item in temp]
        self.Wbb      = np.array(self.Wbb)
        self.Wbb_corr = np.array(self.Wbb_corr)
        logging.info(" Done Computing the Bound - Bound transition rates")
        
        logging.info(" Computing the Bound - Continuum transition rates using %d procs"%(self.n_jobs))
        temp = Parallel(n_jobs=self.n_jobs,backend='multiprocessing') (delayed(self.compute_cont_row)(pot,eigen,m) for m in range(pot.n_bound))
        self.Wbe      = [item[0] for item in temp]
        self.Wbe_corr = [item[1] for item in temp]
        self.Wbe      = np.array(self.Wbe)
        self.Wbe_corr = np.array(self.Wbe_corr)
        
        self.Wbc      = np.sum(self.Wbe,axis=1)
        self.Wbc_corr = np.sum(self.Wbe_corr,axis=1)
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
            
        if(self.use_Filon):
            if(len(pot.t) % 2 == 0):
                print(" ERROR : Odd number of time steps necessary for Filon quadrature implementation. ")
                sys.exit()
            dt       = pot.t[1] - pot.t[0]
            theta    = om_nm*dt            
            
        if(method == 1):
            Vnm_corr = acorr(Vnm,norm=False)/len(Vnm)
            if (self.use_Filon == False):
                wnm = 2 * integrate.trapz(Vnm_corr*np.cos(om_nm*pot.t),pot.t)
            else:
                if(om_nm == 0):
                    wnm = 2 * integrate.trapz(Vnm_corr*np.cos(om_nm*pot.t),pot.t)
                else:
                    alpha_filon, beta_filon, gamma_filon = self.get_filon_params(theta)

                    temp = Vnm_corr * np.cos(om_nm*pot.t)
                    C_e = sum(temp[::2]) - (temp[-1] + temp[0])/2
                    C_o = sum(temp[1::2])
    
                    wnm  = alpha_filon * ( Vnm_corr[-1]*np.sin(om_nm*pot.t[-1]) - Vnm_corr[0]*np.sin(om_nm*pot.t[0]) )
                    wnm += beta_filon  * C_e + gamma_filon * C_o
    
                    wnm = abs(wnm * 2 * dt)
                
        elif(method == 2):
            psd = 0
            if(self.use_Filon == False or om_nm == 0):
                dt = pot.t[1] - pot.t[0]
                psd = np.sum(Vnm*np.exp(-1j*om_nm*pot.t))
                wnm = (abs(psd)*dt)**2/(pot.t[-1] - pot.t[0])
            else:
                alpha_filon, beta_filon, gamma_filon = self.get_filon_params(theta)
                temp_s = Vnm * np.sin(om_nm*pot.t)
                temp_c = Vnm * np.cos(om_nm*pot.t)

                C_e = sum(temp_c[::2])  - (temp_c[-1] + temp_c[0])/2
                C_o = sum(temp_c[1::2])
                
                I_c  = alpha_filon* (Vnm[-1]*np.sin(om_nm*pot.t[-1]) - Vnm[0]*np.sin(om_nm*pot.t[0]))
                I_c += beta_filon*C_e + gamma_filon*C_o
                I_c *= dt

                S_e = sum(temp_s[::2])  - (temp_s[-1] + temp_s[0])/2
                S_o = sum(temp_s[1::2])
                
                I_s  = -1*alpha_filon* (Vnm[-1]*np.cos(om_nm*pot.t[-1]) - Vnm[0]*np.cos(om_nm*pot.t[0]))
                I_s += beta_filon*S_e + gamma_filon*S_o
                I_s *= dt
                
                wnm = (I_c**2 + I_s**2)/(pot.t[-1] - pot.t[0])                
            
            
        if(self.QCF_Model != 'Egl'):
            wnm_corr = wnm * self.compute_QCF(om_nm)
        else:
            wnm_corr = wnm * self.compute_QCF_Egl(om_nm,wnm,Vnm_corr,pot.t)
            
        wnm      = wnm/(hbar**2)
        wnm_corr = wnm_corr/(hbar**2)
        
        return wnm, wnm_corr
    
    def get_filon_params(self, theta):
         alpha_filon = 1/theta + np.sin(2*theta)/(2*theta**2) - 2*(np.sin(theta))**2/theta**3
         beta_filon  = 2*((1+(np.cos(theta))**2)/theta**2 - np.sin(2*theta)/theta**3)
         gamma_filon = (4/theta**2)*(np.sin(theta)/theta - np.cos(theta))
         
         return alpha_filon, beta_filon, gamma_filon

    def acorr(v, norm=False):
        nstep = v.shape[0]
        c = np.zeros((nstep,), dtype=float)
        _norm = 1 if norm else 0

        # Correlation via fft. After ifft, the imaginary part is (in theory) =                                                                                                                                                              
        # 0, in practise < 1e-16, so we are safe to return the real part only.                                                                                                                                                              
        vv = np.concatenate((v, np.zeros((nstep,),dtype=float)))
        c = ifft(np.abs(fft(vv))**2.0)[:nstep].real

        if norm:
            return c / c[0]
        else:
            return c

    def compute_bound_row(self,pot,eigen,m):
        Wn      = np.zeros((pot.n_bound,))
        Wn_corr = np.zeros_like(Wn)
        
        for n in range(m+1,pot.n_bound):
            Wn[n], Wn_corr[n] = self.compute_ind_rate(pot,eigen,m,n,method=2)
        return Wn, Wn_corr
         
    def compute_cont_row(self,pot,eigen,m):
        Wn      = np.zeros((pot.n_cont,))
        Wn_corr = np.zeros_like(Wn)
        
        for c in range(pot.n_bound, pot.n_bound+pot.n_cont):
            Wn[c-pot.n_bound], Wn_corr[c-pot.n_bound] = self.compute_ind_rate(pot,eigen,m,c,method=2)
        return Wn, Wn_corr
    
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
                self.Wbb[n][m]      = self.Wbb[m][n] * np.exp(-self.beta*om)
                self.Wbb_corr[n][m] = self.Wbb_corr[m][n] * np.exp(-self.beta*om)
                
    def solve_ME(self,pot):
        self.P0 = np.exp(-self.beta*pot.E_bound)
        self.P0 = self.P0/np.sum(self.P0)
        
        self.W      = np.zeros((pot.n_bound,pot.n_bound))
        self.W_corr = np.zeros_like(self.W)
        
        for n in range(pot.n_bound):
            for k in range(pot.n_bound):
                self.W[n][n]      += self.Wbb[n][k]
                self.W_corr[n][n] += self.Wbb_corr[n][k]
                
            self.W[n][n]      += self.Wbc[n]
            self.W_corr[n][n] += self.Wbc_corr[n]

            for m in range(pot.n_bound):
                self.W[n][m]      -=  self.Wbb[m][n]
                self.W_corr[n][m] -= self.Wbb_corr[m][n]

        self.tau_des = np.sum(np.linalg.inv(self.W)@self.P0)
        self.k_des   = 1/self.tau_des

        self.tau_des_corr = np.sum(np.linalg.inv(self.W_corr)@self.P0)
        self.k_des_corr   = 1/self.tau_des_corr

        logging.info("\n   Overall desorption rate constant [1/s]               = %6.8E"%(self.k_des*s2au))        
        logging.info("\n   Overall QCF corrected desorption rate constant [1/s] = %6.8E"%(self.k_des_corr*s2au))
        
    def Write_transition_rates(self,Input,pot,pref):
        
        out_bb = Input.Out_Dir + pref + ".Wbb"
        out_bc = Input.Out_Dir + pref + ".Wbc"
        out_be = Input.Out_Dir + pref + ".Wbe"
        
        file_bb = open(out_bb,"w")
        file_bc = open(out_bc,"w")
        file_be = open(out_be,"w")

        if(self.QCF_Model != "None"):
            out_bb_corr = Input.Out_Dir + pref + "-" + self.QCF_Model + ".Wbb"
            out_bc_corr = Input.Out_Dir + pref + "-" + self.QCF_Model + ".Wbc"
            out_be_corr = Input.Out_Dir + pref + "-" + self.QCF_Model + ".Wbe"
            
            file_bb_corr = open(out_bb_corr,"w")
            file_bc_corr = open(out_bc_corr,"w")
            file_be_corr = open(out_be_corr,"w")

                                                
        logging.info(" Converting the transition rates to 1/s.")
        self.Wbb *= s2au
        self.Wbe *= s2au
        self.Wbc *= s2au

        self.Wbb_corr *= s2au
        self.Wbe_corr *= s2au
        self.Wbc_corr *= s2au
        
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

        if(self.QCF_Model != "None"):
            # Bound - Bound 
            for m in range(pot.n_bound):
                file_bb_corr.write("%4.6E "%(pot.E_bound[m]/eV2Ha))
            file_bb_corr.write("\n")

            for m in range(pot.n_bound):
                for n in range(pot.n_bound):
                    file_bb_corr.write("%4.6E "%(self.Wbb_corr[m][n]))
                file_bb_corr.write("\n")

            # Bound - Continuum
            for e in range(pot.n_cont):
                file_be_corr.write("%4.6E "%(pot.E_cont[e]/eV2Ha))
            file_be_corr.write("\n")
        
            for m in range(pot.n_bound):
                file_bc_corr.write("%02d  %4.6E \n"%(m,self.Wbc_corr[m]))
                for e in range(pot.n_cont):
                    file_be_corr.write("%4.6E "%(self.Wbe_corr[m][e]))
                file_be_corr.write("\n")

            file_bb_corr.close()
            file_bc_corr.close()
            file_be_corr.close()
        
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

        
        Vnm_corr = acorr(Vnm,norm=False)/len(Vnm) * np.cos(om_nm*pot.t)
        plt.figure(10)
        plt.plot(pot.t/ps_au,Vnm_corr,'-o',label="(%d,%d)"%(m,n))
        plt.legend()
        plt.xlabel('t [ps]')
        
    def Test_Vmn(self,pot,eigen,m,n):
        eig_n = eigen.get_eigenstate(pot,n)
        eig_m = eigen.get_eigenstate(pot,m)
        om_nm = (pot.Energy[n] - pot.Energy[m])/hbar

        Vnm = np.zeros(pot.Vf.shape[0])
        for i in range(len(Vnm)):
            f = eig_n * pot.Vf[i] * eig_m
            Vnm[i] = integrate.simps(f,pot.z)
        
        plt.figure(11)
        plt.plot(pot.t/ps_au,Vnm,'.',label="(%d,%d)"%(m,n))
        plt.legend()
        plt.xlabel('t [ps]')

