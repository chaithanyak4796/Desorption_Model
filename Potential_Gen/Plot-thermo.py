import numpy as np
import matplotlib.pyplot as plt

Dt = [0.5]   #fs

plt.figure(figsize=(16,9))
for dt in Dt:
    #fname = "thermo_" + str(dt) +"fs.dat"
    fname = "thermo.dat"
    Data = np.loadtxt(fname,skiprows=2)

    Temp = Data[:,1]
    En   = Data[:,2]

    t = Data[:,0] * dt * 0.001  # ps
    
    l = int(len(Temp)/2)
    print("Avg_Temp = ",np.mean(Temp[l:]))
    print("Std_Temp = ",np.std(Temp[l:]))

    #print("Avg_TE   = ",np.mean(En[l:]))
    #print("Std_TE   = ",np.std(En[l:]))

    zEn = np.polyfit(t[l:],En[l:],1)
    zT  = np.polyfit(t[l:],Temp[l:],1)
    print("slope in TE   = ",zEn[0])
    print("slope in Temp = ",zT[0])

    plt.subplot(121)
    plt.plot(t,Temp,label=str(dt))
    plt.legend()

    plt.subplot(122)
    plt.plot(t,En,label=str(dt))
    plt.legend()

plt.suptitle("reaxFF_CHO")
plt.show()
