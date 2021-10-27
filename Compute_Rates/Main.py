#!/usr/bin/env python3

# Python modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import logging

# Add my modules to path
modulepath = "./"
sys.path.append(modulepath)

# My modules
from Params import *
import Input
import eigen_states_DVR
import Hamiltonian
import Matrix

plt.close('all')

if(len(sys.argv) == 1):
    Input_fname = './Input_files/Input.cfg'
elif(len(sys.argv) == 2):
    Input_fname = str(sys.argv[1])
else:
    print("ERROR : Too many inputs provided. Provide only the Input file name.")
    sys.exit()

start = time.time()

######## Initialize the Input class ######
Input = Input.Input(Input_fname)
print(Input.Inp_Dir)
#sys.exit()

##### Create the logger ###########
Log_fname=Input.Log_fname
#logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logging.basicConfig(filename=Log_fname, filemode='w', level=logging.INFO, format='%(levelname)s: %(message)s')
logging.info('Begining the run')

######## Create the Eigen class ###########
eigen = eigen_states_DVR.Eigen(Input)

####### Solve the zero-order problem ########
pot = Hamiltonian.Potentials(Input)
pot.solve_H0(eigen)
print(" Done solving the TISE.")

###### Initialize the matrices ############
soln = Matrix.Matrix(pot,Input)
print(" Done initializng the solution matrix.")

##### Compute the transition rates ########
soln.compute_transition_rates(Input,pot,eigen)

##########################################  
end = time.time()
logging.info("\n Total run time [s] = %.2f"%(end - start) )
print("\n Total run time [s] = %.2f"%(end - start))
