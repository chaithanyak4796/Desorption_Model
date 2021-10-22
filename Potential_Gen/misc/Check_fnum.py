import numpy as np

fname = "pot_fnum.dump"

Data = np.loadtxt(fname,dtype='int')
#Data = int(Data)
print(Data[0])
N = np.arange(2001,4000)
miss = [] 

for i in N:
    if i in Data:
        continue
    else:
        miss.append(i)

miss = np.array(miss)
print("%d elements missing"%(len(miss)))
print(miss) 
