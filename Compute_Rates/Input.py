import numpy as np
from Params import *
from pathlib import Path
from configparser import ConfigParser, ExtendedInterpolation
import os

class Input:
    def __init__(self,fname):
        self.Inp_fname = fname
        parser = ConfigParser(os.environ,interpolation=ExtendedInterpolation())
        parser.read(fname)

        ###### General settings #############
        self.Temp        = 300     # Default value. Can be overridden in [POTENTIALS]
        self.system      = str(parser['GENERAL']['System'])

        self.Inp_Dir   = str(parser['GENERAL']['Inp_Dir'])
        self.Out_Dir   = str(parser['GENERAL']['Out_Dir'])
        self.Log_fname = str(parser['GENERAL']['Log_fname'])
        self.Out_pref  = str(parser['GENERAL']['Out_pref'])

        ######## Desorption settings ##############
        
        self.Temp     = float(parser['POTENTIALS']['Temp'])
        self.mass     = float(parser['POTENTIALS']['Mass'])
        self.dt       = float(parser['POTENTIALS']['dt'])  
        self.t_max    = float(parser['POTENTIALS']['t_max'])
        
        self.compute_pot = parser['POTENTIALS'].getboolean('Compute_pot')  
        self.save_pot    = parser['POTENTIALS'].getboolean('Save_pot')
        
        self.pot_prefix       =  str(parser['POTENTIALS']['Pot_prefix']) 
        self.Merged_Pot_fname =  str(parser['POTENTIALS']['Merged_Pot'])

        self.n_skip       = int(parser['POTENTIALS']['n_skip'])
        self.interp_t     = float(parser['POTENTIALS']['interp_t'])
        self.pad_z        = float(parser['POTENTIALS']['pad_z'])
        self.interp_z     = float(parser['POTENTIALS']['interp_z'])
        self.interp_z_int = float(parser['POTENTIALS']['interp_z_int'])
        
        ######## MISC Options #############
    
        self.E_max       = float(parser['MISC']['E_max']) * eV2Ha
        self.NProcs      = int(parser['MISC']['NProcs'])
        self.QCF_Model   = str(parser['MISC']['QCF_Model'])

        ### Finally create  the temp directory if it does not exist    
        Path(self.Out_Dir).mkdir(parents=True, exist_ok=True)  # mkdir -p self.Out_Dir
