import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from Params import *

Dir = "/media/chaithanya/Chaithanya_New/Surface_Chem/Desorption/Density_Matrix/Results/Oxygen/Bridge_site/Test_dz/Method_1/Test/"

fname_en = Dir + "Eigen.energy"
fname_st = Dir + "Eigen.states"

data_en = np.loadtxt(fname_en)
data_st = np.loadtxt(fname_st)

En = data_en[:,1]

z      = data_st[0]
V0     = data_st[1]

states = data_st[2:]

nstates = min(50, len(states)-1)

if(1):
    plt.figure(3)
    # plt.plot(z*Bo2Angs, V0, 'k')
    plt.plot(z*Bo2Angs, 1*states[-1],'.')
    plt.xlim([0,4])

if(0):
    ax  = plt.gca()
    out = ax.plot(z*Bo2Angs,V0,'k')

    for m in range(nstates):
        vec = states[m]
        c   = max(vec)

        if(c < 1E-4): c = min(vec)
        vec = vec/c
        amp = En[m+1] - En[m]

        if(En[m] <= V0[-1]):
            ax.plot(z*Bo2Angs, 0.3*amp*vec + En[m], '-', linewidth=2, label=r"$n=%i$"%m)
        else:
            ax.plot(z*Bo2Angs, 0.3*amp*vec + En[m], '-', linewidth=2, label=r"$n=%i$"%m)

    min_V = min(min(V0),0)
    plt.ylim(min_V, En[m] + amp)
    ax.axes.yaxis.set_ticks([])
    plt.xlabel(r"$z  [\AA]$")
    plt.ylabel(r"$\psi(z)$")
    plt.show()

