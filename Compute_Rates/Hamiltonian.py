import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from Params import *
import sys
from scipy import integrate, interpolate
import logging


class Potentials:
    def __init__(self,Input):
        self.system = Input.system
        self.setup_Desorption(Input)

        self.E_inf   = self.V0[-1]
        self.E_bound = []
        self.E_cont  = []
        self.rho_E   = []
        self.E_max   = Input.E_max

    def setup_Desorption(self,Input):
        """                          
             System : Adatom-surface interactions in a thermally averaged potential with thermal fluctuations                                           
        """

        self.Inp_Dir = Input.Inp_Dir
        self.mass = Input.mass * mp  # Converting from amu to Hartree units

        if(Input.compute_pot):
            self.z, self.V0, self.Vf = self.Calc_avg_fluc(Input)
        else:
            self.z, self.V0, self.Vf = self.Read_avg_fluc(Input)

        self.z  = self.z / Bo2Angs
        self.V0 = self.V0 * eV2Ha
        self.Vf = self.Vf * eV2Ha

        self.dt = Input.dt
        self.dt = self.dt * (Input.n_skip+1)
        self.dt = self.dt * Input.interp_t

        self.z  = np.flip(self.z)
        self.V0 = np.flip(self.V0)
        self.Vf = np.flip(self.Vf,axis=-1)

        self.dt = self.dt * 1E-3  # Converting from fs to ps
        self.t = np.arange(0,self.Vf.shape[0]) * self.dt * ps_au

        if(len(self.t) % 2 == 0):
            logging.info(" Found even number of time steps. Discarding the final time step.")
            self.t  = self.t[:-1]
            self.Vf = self.Vf[:-1,:]

        logging.info(" Potentials : dz    [A]  = %.4f"%((self.z[1]-self.z[0])*Bo2Angs))
        logging.info(" Potentials : dt    [ps] = %.4f"%(self.dt))
        logging.info(" Potentials : t_max [ps] = %.4f"%(self.t[-1]/ps_au))

    def Calc_avg_fluc(self, Input):
        logging.info("\n Computing the Average and fluctuating potentials.")
        logging.info(" Name of the Merged Pot file = %s" %(Input.Merged_Pot_fname))

        max_steps = int(np.ceil(Input.t_max*1000/Input.dt)) + 1

        Data = np.loadtxt(Input.Merged_Pot_fname, usecols=(1,2))
        num_z = len( np.unique(Data[:,0]) )
        num_t = len(Data[:,0])//num_z

        logging.info("Total Number of time steps in File = %d"%(num_t))
        if(Input.pad_z > 0):
            z = Data[:num_z,0]
            if(Input.pad_z <= z[0]):
                logging.warning(" The padding limit is less than z_max. Skipping padding")
            logging.info("Padding the V and z arrays to z_max = %.4f" %(Input.pad_z))
            dz = z[0] - z[1]
            logging.info ("Original: z[0] = %f, z[-1] = %f, dz = %f"%(z[0],z[-1],dz))
            z = np.arange(Input.pad_z,z[-1],-1*dz)
            logging.info ("New     : z[0] = %f, z[-1] = %f, dz = %f"%(z[0],z[-1],dz))
            num_z_new = len(z)
            #print ("Length of new z array = ", num_z_new)
            V = np.zeros((num_t,num_z_new))
            temp = np.reshape(Data[:,1],(num_t,num_z))
            V[:,-num_z:] = temp
            for i in range(num_t):
                V[i,0:num_z_new-num_z+1] = V[i,num_z_new-num_z+1]
            num_z = num_z_new
        else:
            V = np.reshape(Data[:,1],(num_t,num_z))
            z = Data[:num_z,0]

        num_t = min(num_t,max_steps)
        V     = V[:num_t]

        logging.info("Total Number of time steps used = %d" %(num_t))

        if(Input.n_skip > 0):
            logging.info("Skipping every %d time steps"%(Input.n_skip))
            idx = np.arange(0,len(V[:,0]),Input.n_skip+1)
            V   = V[idx,:]
        if(Input.interp_t < 1.0):
            logging.info("Interpolating time data with a factor : %.4f" %(Input.interp_t))
            t    = np.arange(0,V.shape[0])
            tp   = np.arange(0,V.shape[0],Input.interp_t)[:-1*int(1.0/Input.interp_t)+1]
            temp = np.zeros((len(tp),V.shape[1]))
            for i in range(V.shape[1]):
                y  = V[:,i]
                f  = interpolate.interp1d(t,y)
                temp[:,i] = f(tp)
            V = temp
        if(Input.interp_z < 1.0):
            logging.info("Interpolating space data with a factor : %.4f" %(Input.interp_z))
            x   = np.arange(0,V.shape[1])
            xp  = np.arange(0,V.shape[1],Input.interp_z)[:-1*int(1.0/Input.interp_z)+1]
            temp = np.zeros((V.shape[0],len(xp)))
            for i in range(V.shape[0]):
                y = V[i]
                f = interpolate.interp1d(x,y)
                temp[i] = f(xp)
            V = temp
            f = interpolate.interp1d(x,z)
            z = f(xp)
            num_z = len(z)
            

        V_avg  = np.zeros((num_z))
        V_fluc = np.zeros_like(V)

        for i in range(V.shape[0]):
            V[i] = V[i] - V[i][0]
        V_avg = np.mean(V,axis=0)
        for i in range(V.shape[0]):
            V_fluc[i] = V[i] - V_avg

        logging.info("Done computing potentials\n")

        if(Input.save_pot):
            logging.info("Storing the computed potentials to disk")
            fname_avg  = Input.pot_prefix + ".avg"
            fname_fluc = Input.pot_prefix + ".fluc"
            fname_rms  = Input.pot_prefix + ".rms"
            logging.info("Avg  pot file name : %s"%(fname_avg))
            logging.info("Fluc pot file name : %s"%(fname_fluc))
            logging.info("RMS  pot file name : %s"%(fname_rms))

            f_avg  = open(fname_avg,"w")
            f_fluc = open(fname_fluc,"w")
            f_rms  = open(fname_rms,"w")

            V_rms = (np.mean(V_fluc**2,axis=0))**0.5
            
            for i in range(V_avg.shape[0]):
                f_avg.write("%4.6f %6.8f\n"%(z[i],V_avg[i]))
                f_rms.write("%4.6f %6.8f\n"%(z[i],V_rms[i]))
            for i in range(V_fluc.shape[0]):
                for j in range(V_fluc.shape[1]):
                    f_fluc.write("%6.8f "%(V_fluc[i][j]))
                f_fluc.write("\n")

            f_avg.close()
            f_fluc.close()
            f_rms.close()
            
        return z,V_avg,V_fluc

    def Read_avg_fluc(self, Input):
        logging.info("\n Reading the Average and fluctuating potentials")
        fname_avg  = Input.pot_prefix + ".avg"
        fname_fluc = Input.pot_prefix + ".fluc"

        data_avg = np.loadtxt(fname_avg)
        z  = data_avg[:,0] 
        V  = data_avg[:,1] 
        Vf = np.loadtxt(fname_fluc)

        if(Input.pad_z > 0):
            if(Input.pad_z <= z[0]):
                logging.warning(" The padding limit is less than z_max. Skipping padding")
            else:
                logging.info("Padding z,V,Vf arrays to z_max = %.2f"% (Input.pad_z) )
                dz = z[0] - z[1]
                logging.info("Original: z[0] = %f, z[-1] = %f, dz = %f"%(z[0],z[-1],dz))
                z = np.arange(Input.pad_z,z[-1],-1*dz)
                logging.info("New     : z[0] = %f, z[-1] = %f, dz = %f"%(z[0],z[-1],dz))
                num_z_new = len(z)

                V_new  = np.zeros((num_z_new,))
                Vf_new = np.zeros_like(V)

                V_new[:,-num_z]  = V
                Vf_new[:,-num_z] = Vf

                V_new[0:num_z_new-num_z+1]  = V_new[num_z_new-num_z+1]
                Vf_new[0:num_z_new-num_z+1] = Vf_new[num_z_new-num_z+1]

                V  = V_new
                Vf = Vf_new
        
        if(Input.n_skip > 0):
            logging.info("Skipping every %d time steps"%(Input.n_skip))
            idx = np.arange(0,len(Vf[:,0]),Input.n_skip+1)
            Vf   = Vf[idx,:]
        if(Input.interp_t < 1):
             logging.info("Interpolating time data with a factor of %.4f" %(Input.interp_t) )
             t    = np.arange(0,Vf.shape[0])
             tp   = np.arange(0,Vf.shape[0],Input.interp_t)[:-1*int(1.0/Input.interp_t)+1]
             temp = np.zeros((len(tp),Vf.shape[1]))
             for i in range(Vf.shape[1]):
                 y  = Vf[:,i]
                 f  = interpolate.interp1d(t,y)
                 temp[:,i] = f(tp)
             Vf = temp
        if(Input.interp_z < 1.0):
            logging.info("Interpolating space data with a factor : %.4f" %(Input.interp_z))
            x   = np.arange(0,Vf.shape[1])
            xp  = np.arange(0,Vf.shape[1],Input.interp_z)[:-1*int(1.0/Input.interp_z)+1]
            temp = np.zeros((Vf.shape[0],len(xp)))
            for i in range(Vf.shape[0]):
                y = Vf[i]
                f = interpolate.interp1d(x,y)
                temp[i] = f(xp)
            Vf = temp
            f  = interpolate.interp1d(x,V)
            V  = f(xp)
            f = interpolate.interp1d(x,z)
            z = f(xp)
            num_z = len(z)

        return z,V,Vf

    def solve_H0(self,eigen):
        
        H_mat = eigen.build_hamiltonian(self)
        self.evals, self.evecs = eigh(H_mat)
        logging.info(f"\n First few energy eigen values [eV]:{self.evals[:10]/eV2Ha} ")
        
        ##### Plot the eigenstates ######
        #eigen.plot_eigenstates(self, 40)
        eigen.normalize_eigenstates(self)
        #plt.show()
        
        ##### Get the bound and continuum states and density of continuum states
        for i in range(len(self.evals)):
            if(self.evals[i] < self.E_inf):
                self.E_bound.append(self.evals[i])
            elif(self.evals[i] >= self.E_inf and self.evals[i] <= self.E_max):
                self.E_cont.append(self.evals[i])
                if(len(self.rho_E) == 0):
                    self.rho_E.append(1/(self.evals[i+1]-self.evals[i]))    # Forward differencing at lowest energy
                else:
                    self.rho_E.append(2/(self.evals[i+1]-self.evals[i-1]))  # Central differencing @ higher energies

        self.E_bound = np.array(self.E_bound)
        self.E_cont  = np.array(self.E_cont)
        self.rho_E   = np.array(self.rho_E)
    
        self.Energy = np.concatenate((self.E_bound,self.E_cont))
        self.n_bound = len(self.E_bound)
        self.n_cont  = len(self.E_cont)

        logging.info("\n Maximum value of continuum energy [eV] = %.4f"%(self.E_max/eV2Ha))
        logging.info(" Total # of bound     levels = %d" %(self.n_bound))
        logging.info(" Total # of continuum levels = %d" %(self.n_cont))
        logging.info(" Total # of levels = %d" %(self.n_bound + self.n_cont))
        logging.info(" ")
    
