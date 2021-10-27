import numpy as np
import sys

def Read_Lat_Const(fname):
    pot_name = []
    Lat_const = []
    num_pot = 0
    fr = open(fname,'r')
    line = fr.readline().strip("\n")
    while(line != 'End'):
        if(line [0] == '*'):
            pot = line.strip(" *\n")
            # print (pot)
            pot_name.append(pot)
            num_pot += 1
            Lat_const.append([''])
            line = fr.readline().strip("\n") 
            while(True):
                line = fr.readline().strip("\n") 
                if(len(line) > 0): 
                    # print(line)
                    temp = [float(j) for j in line.split()]
                    Lat_const[num_pot-1].append(temp)
                else:
                    break
        
        line = fr.readline().strip("\n")

    pot_name = np.array(pot_name)
    for i in range(num_pot) :
        Lat_const[i].pop(0)
        Lat_const[i] = np.array(Lat_const[i])
    
    return pot_name,Lat_const

class Particle:
    def __init__(self, idx=1, a_type=1, a_name='C'):
        self.pos  = np.zeros(3)
        self.type = a_type
        self.name = a_name
        self.idx  = idx

def Get_site_idx(site, num_rep, particle_list, info_fname, model):
    if(site == 'Bridge' or site == 'bridge'):
        site = "bridge"
        site_idx = np.zeros(2,dtype=int)
        i = num_rep[0]//2 - 1
        j  = num_rep[1]//2 - 1
        k = num_rep[2] 
        base_idx = (i*num_rep[1]*num_rep[2] + j*num_rep[2] + k) * 4
        site_idx[0] = base_idx - 2 + 1
        site_idx[1] = base_idx - 1 + 1
        for i in range(len(site_idx)):
            print(site_idx[i])
        print(particle_list[site_idx[0] - 1].pos)
        print(particle_list[site_idx[1] - 1].pos)

    elif(site == "Top" or site == 'top'):
        site = "top"
        site_idx = np.zeros(1,dtype=int)
        i = num_rep[0]//2 - 1
        j = num_rep[1]//2 - 1
        k = num_rep[2]
        base_idx = (i*num_rep[1]*num_rep[2] + j*num_rep[2] + k) * 4
        site_idx[0] = base_idx - 2 + 1
        print(site_idx)
        print(particle_list[site_idx[0] - 1].pos)
    else:
        print(" Error : invalid site type.")
        sys.exit(0)

    fw = open(info_fname,"w")
    fw.write("%d\n"%(model))            # Interaction model
    fw.write("%s\n"%(site))             # Site type
    for i in range(len(site_idx)):
        fw.write("%d\n"%(site_idx[i]))  # Lattice indices
    for i in range(len(site_idx)):      # Lattice coordinates
        fw.write("%.6f  %.6f  %.6f \n"%(particle_list[site_idx[i] - 1].pos[0],particle_list[site_idx[i] - 1].pos[1],particle_list[site_idx[i] - 1].pos[2]))
    
    fw.close()

    return site, site_idx
