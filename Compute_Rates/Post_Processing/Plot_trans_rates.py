import numpy as np
import matplotlib.pyplot as plt
# plt.close('all')

Dir = "/nobackupp2/ckondur/Desorption/Desorption_Model/Results/Oxygen/Bridge_site/Ads-no/"

Temp = np.array([300,400,600,800,1000])
Des_rate = np.zeros_like(Temp,dtype=float)

plot_cont_en = False
n_en = [0,5,10,20,25]

kb_eV = 8.617333262E-5  # eV/K
au2s  = 2.4188843265857E-17
s2au  = 1/au2s 

for i in range(len(Temp)):
    T = Temp[i]
    pref = "T_" + str(T) + "K"
    fname_bb = Dir + pref + ".Wbb"
    fname_bc = Dir + pref + ".Wbc"

    data_bb = np.loadtxt(fname_bb)
    E_bound = data_bb[0]
    n_bound = len(E_bound)
    Wbb = (data_bb[1:])

    data_bc = np.loadtxt(fname_bc)
    Wbc     = (data_bc[:,1]) 

    bs_no = np.arange(n_bound)

    if(0):
        plt.figure(1)
        b_jump = 1
        plt.semilogy(bs_no[:-b_jump],np.diagonal(Wbb,b_jump),label=pref)
        #plt.semilogy(bs_no[:-2],np.diagonal(Wbb,2),label='n->n+2')
        plt.legend()
        plt.xlabel('n')
        plt.ylabel(r'$W_{n\rightarrow m} [s^{-1}]$')

        plt.figure(2)
        plt.semilogy(bs_no,Wbc,label=pref)
        plt.legend()
        plt.xlabel('n')
        plt.ylabel(r'$W_{n\rightarrow cont} [s^{-1}]$')
    
    if(plot_cont_en):
        fname_be = Dir + pref + ".Wbe"
        data_be  = np.loadtxt(fname_be)
        E_cont = data_be[0]
        
        for j in range(len(n_en)):
            n = n_en[j]
            plt.semilogy(E_cont,data_be[n+1],label=str(n))
        plt.xlabel("E_cont [eV]")
        plt.legend()
        
    
    beta = 1/kb_eV/T
    P0 = np.exp(-beta*E_bound)
    P0 = P0/sum(P0)
    
    W = np.zeros((n_bound,n_bound))
    for n in range(n_bound):
        for k in range(n_bound):
            W[n][n] += Wbb[n][k]
        W[n][n] += Wbc[n]
        for m in range(n_bound):
            W[n][m] -= Wbb[m][n]
    tau   = np.sum(np.linalg.inv(W)@P0)
    k_des = 1/tau 
    print("k_des [1/s] = %.4E"%(k_des))
    
    Des_rate[i] = k_des

plt.figure(5)
plt.semilogy(1000/Temp,Des_rate,'-o',label='Semi-classical Model-One_Sided')
plt.xlabel('1000/T')
plt.ylabel('Rate [1/s]')


if(0):
    Model = ['Gauss_25_Renorm', 'Lorentzian_100_Bare', 'Lorentzian_100_Renorm']
    for i in range(len(Model)): 
        fname = '../../Data/CO_Cu/' + Model[i] + '.csv'
        data = np.loadtxt(fname,delimiter=',')
        plt.semilogy(data[:,0],np.exp(data[:,1]),'o',label=Model[i])
    plt.legend()
    # plt.ylim([1E-4,1E11])

plt.show()

if(len(Temp)>1):
    log_rate = np.log(Des_rate)
    inv_temp = 1/(kb_eV*Temp)
    plt.figure(6)
    plt.plot(inv_temp,log_rate)
    
    kin = np.polyfit(inv_temp,log_rate,1)
    print("E_A [eV] = %.4f"%(-1*kin[0]))
    print("A  [1/s] = %6.4E"%(np.exp(kin[1])))
