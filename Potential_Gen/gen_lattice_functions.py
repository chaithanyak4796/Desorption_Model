import numpy as np
import sys

kb = 1.3806E-23  # J/K
NA = 6.022E23

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
    def __init__(self, mass=12.0107, idx=1, a_type=1, a_name='C'):
        self.pos  = np.zeros(3)
        self.vel  = np.zeros(3)
        self.type = a_type
        self.name = a_name
        self.idx  = idx
        self.mass = mass

def Get_site_idx(site, num_rep, particle_list, info_fname, model):
    if(site == 'Bridge' or site == 'bridge'):
        site = "bridge"
        site_idx = np.zeros(2,dtype=int)
        i = num_rep[0]//2 - 1
        j = num_rep[1]//2 - 1
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
        
    elif(site == "Edge" or site == "edge"):
        site = "edge"
        site_idx = np.zeros(1,dtype=int)
        return site, site_idx


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

def Get_edge_site_idx(particle_list, info_fname, model, top_layer, lat_a, pit_D, pit_H, center, edge_type):
     # Find the edge sites
    surr = []
    edge = []
    zig_zag   = []
    arm_chair = []
    
    # Sweep the top layer to identify atoms close to the edge sites
    for j in top_layer:
        pos = particle_list[j].pos
        extra_a = lat_a
        pit = np.array([pit_D/2+lat_a,pit_D/2+lat_a,pit_H+2])
        pos = (pos-center)/pit
        dis = np.linalg.norm(pos)
        if(dis <= 1.0):
            particle_list[j].name = 'C'
            surr.append(j)
    print("Number of potential edge sites after first sweep = ", len(surr))
    
    # Idenitfy the edge sites by caomputing the coordination numbers of the identified atoms
    num_edge = 0
    for i in surr:
        coord_num = 0;
        pos_i = particle_list[i].pos
        for j in top_layer:
            pos_j = particle_list[j].pos
            dis = np.linalg.norm(pos_i-pos_j)
            if(dis <= 1.05*lat_a/3**0.5):   # 5% tolerance
                coord_num += 1
        if(coord_num == 3):   # 3 because we're not excluding i=j
            # print(i)
            particle_list[i].name = 'N'
            edge.append(i)
            num_edge += 1
    print("Number of edge sites = ", num_edge)
            
    # Classify the edge sites as zig-zag or arm-chair by computing the partial coordination numbers 
    for i in edge:
        coord_num = 0
        pos_i = particle_list[i].pos
        for j in edge:
            pos_j = particle_list[j].pos
            dis = np.linalg.norm(pos_i-pos_j)
            if(dis <= 1.05*lat_a/3**0.5):   # 5% tolerance
                coord_num += 1
        if(coord_num == 2):
            particle_list[i].name = 'N'
            arm_chair.append(i)
        else:
            particle_list[i].name = 'O'
            zig_zag.append(i)
    print(" Number of zig-zag   configurations = ",len(zig_zag))
    print(" Number of arm-chair configurations = ",len(arm_chair))
    
    site_idx = np.zeros(1,dtype=int)
    if(edge_type == 'zig-zag'):
        site_idx[0] = zig_zag[0] + 1   # LAMMPS idx starts from 1
    else:
        site_idx[0] = arm_chair[0] + 1
        
    print(site_idx)
    print(particle_list[site_idx[0] - 1].pos)
    
    fw = open(info_fname,"w")
    fw.write("%d\n"%(model))            # Interaction model
    fw.write("%s\n"%('edge'))             # Site type
    for i in range(len(site_idx)):
        fw.write("%d\n"%(site_idx[i]))  # Lattice indices
    for i in range(len(site_idx)):      # Lattice coordinates
        fw.write("%.6f  %.6f  %.6f \n"%(particle_list[site_idx[i] - 1].pos[0],particle_list[site_idx[i] - 1].pos[1],particle_list[site_idx[i] - 1].pos[2]))
    
    fw.close()

    return site_idx, zig_zag, arm_chair

def get_velocities(num_rep,mass,Temp):
    """ Returns the velocities array in m/s"""
    num_per_layer = 2 * num_rep[0] * num_rep[1] 
    num_layers = 2*num_rep[2]
    N = num_per_layer * num_layers

    m = mass/1000/NA
    sig = (kb*Temp/m)**0.5

    vel = np.zeros((num_layers,num_per_layer,3))
    vel = np.random.normal(0.0,sig,vel.shape).reshape(vel.shape)

    for i in range(num_layers):  # Zeroing out the com velocity of each layer
        v_com   = np.sum(vel[i],axis=0)/num_per_layer
        vel[i] -= v_com

        #print("v_com  =", np.sum(vel[i],axis=0)/num_per_layer)

    T = m*np.sum(vel**2)/(3*N*kb)
    vel = vel*(Temp/T)**0.5
    #print(T)

    T = m*np.sum(vel**2)/(3*N*kb)
    #print(T)

    return vel
