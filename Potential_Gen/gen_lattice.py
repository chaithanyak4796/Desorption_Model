import numpy as np
import matplotlib.pyplot as plt
from gen_lattice_functions import *
import sys

if(len(sys.argv) == 5):
    Pot  = str(sys.argv[1])
    site = str(sys.argv[2])
    Temp = float(sys.argv[3])
    interaction_model = int(sys.argv[4])
    if (site == 'Edge' or site == 'edge'):
        print("Error: Etch pit diameter and height not given.")
        sys.exit(0)
if(len(sys.argv) == 8):
    Pot  = str(sys.argv[1])
    site = str(sys.argv[2])
    edge_type = str(sys.argv[3])
    Temp = float(sys.argv[4])
    interaction_model = int(sys.argv[5])
    pit_D = float(sys.argv[6])
    pit_H = float(sys.argv[7])
    
elif(len(sys.argv) == 1):
    Pot  = 'COMB3'
    site = "Edge"
    edge_type = "arm-chair"
    Temp = 1000
    interaction_model = 3
    pit_D = 16
    pit_H = 6
    
else:
    print("Error: Incorrect number of arguments.")
    sys.exit(0)
    

num_rep     = [16,16,3]
mark_layers = False   # True : Each graphene layer is assigned a different atom_type
write_xyz   = True
adatom_name = 'O'
mass_adatom = 15.999

vac_height      = 20
Lat_const_fname = './Lat_const.dat'
Pos_prefix      = "graphite_slab"
info_fname      = './info.lattice'
mass_C          = 12.0107

# Get the Lattice constants 
units = 'metal'
if(Pot == 'reaxFF' or Pot == 'reaxFF-CHO' or Pot == 'reaxFF-CHON'):
    units = 'real'

pot_name, Lat_const = Read_Lat_Const(Lat_const_fname)

idx_pot = np.where(pot_name == Pot)[0]
if(len(idx_pot) == 0):
    print(" Error : Potential not found.")
    sys.exit()
idx_pot = idx_pot[0]
idx_Temp = np.where(Lat_const[idx_pot][:,0] == Temp)[0]
if(len(idx_Temp) == 0):
    print("\n Error : Temp not found in the potential.")
    print(" Assigning T = 300 K\n")
   
    idx_Temp = np.where(Lat_const[idx_pot][:,0] == 300)[0]
    #sys.exit() 
idx_Temp = idx_Temp[0]

lat_param = Lat_const[idx_pot][idx_Temp]
print("Lattice params found = ",lat_param)

lat_a = lat_param[1]*3**0.5
lat_c = lat_param[2]*2

# LAMMPS box data 
unit_cell_len = np.array([lat_a, lat_a*3**0.5, 2*lat_c])

lattice_vec       = np.zeros((3,3))
lattice_vec[0][0] = lat_a
lattice_vec[1][0] = 0.5 * lat_a
lattice_vec[1][1] = 0.5 * 3**0.5 * lat_a
lattice_vec[2][2] = lat_c

xyz = np.zeros((3,))
xyz[0] = 0.5*lat_a * num_rep[1] - 0.00001

box = np.zeros((6,))
box[1] = lat_a * num_rep[0]
box[3] = 0.5 * 3**0.5 * lat_a * num_rep[1]
box[4] = -1 * vac_height
box[5] =  1 * vac_height + lat_c * num_rep[2] 

# Generate the atoms pos and vels 
num_atoms = 4 * np.product(num_rep)
particle_list = []
pos = np.zeros((3,))
base_lattice = np.zeros((3,))

vels = get_velocities(num_rep,mass_C,Temp) # m/s
T    = (mass_C/1000/NA)*np.sum(vels**2)/(2*4*np.product(num_rep)*kb)
print ("Temp [K] = ",Temp) 
if (units == 'metal'):
    vels *= 1E-2
else:
    vels *= 1E-5

idx = 1
idx_layer = np.zeros(2*num_rep[2],dtype=int)

for i in range(num_rep[0]):
    for j in range(num_rep[1]):
        for k in range(num_rep[2]):
            base_lattice[0] = i*lattice_vec[0][0] + j*lattice_vec[1][0] + k*lattice_vec[2][0]
            base_lattice[1] = i*lattice_vec[0][1] + j*lattice_vec[1][1] + k*lattice_vec[2][1]
            base_lattice[2] = i*lattice_vec[0][2] + j*lattice_vec[1][2] + k*lattice_vec[2][2]
            
            pos = np.copy(base_lattice)
            vel = vels[2*k][idx_layer[2*k]]
            idx_layer[2*k] += 1
            #atom_type = 1
            atom_type = 2*k + 1
            particle_list.append(Particle(mass_C, idx, atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            particle_list[-1].vel = np.copy(vel)

            pos[0] = base_lattice[0] + 0.5*lat_a
            pos[1] = base_lattice[1] + (0.5/3**0.5) * lat_a
            pos[2] = base_lattice[2] + 0
            vel = vels[2*k][idx_layer[2*k]]
            idx_layer[2*k] += 1
            particle_list.append(Particle(mass_C, idx+1, atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            particle_list[-1].vel = np.copy(vel)

            pos[0] = base_lattice[0] + 0.5*lat_a
            pos[1] = base_lattice[1] + (0.5/3.0**0.5) * lat_a
            pos[2] = base_lattice[2] + 0.5*lat_c
            vel = vels[2*k+1][idx_layer[2*k+1]]
            idx_layer[2*k+1] += 1
            atom_type = 2*k + 2
            particle_list.append(Particle(mass_C, idx+2, atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            particle_list[-1].vel = np.copy(vel)
            
            pos[0] = base_lattice[0] + 1.0*lat_a
            pos[1] = base_lattice[1] + (1.0/3.0**0.5) * lat_a
            pos[2] = base_lattice[2] + 0.5*lat_c
            vel = vels[2*k+1][idx_layer[2*k+1]]
            idx_layer[2*k+1] += 1
            particle_list.append(Particle(mass_C, idx+3,atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            particle_list[-1].vel = np.copy(vel)

            z_max = np.copy(pos[2])
            idx += 4

num_atoms = len(particle_list)

# Adjust the z coordinates
box[4] -= z_max
box[5] -= z_max
for i in range(num_atoms):
    particle_list[i].pos[2] -= z_max

if (mark_layers == False):
    for i in range(num_atoms):
        particle_list[i].type = 1
        if (particle_list[i].pos[2] == 0): 
            particle_list[i].type = 2
            atom_type = 2

top_layer = []
# Create an etch pit
if (site == 'Edge' or site == 'edge'):
    print("Original number of atoms = ",num_atoms)
    center = np.array([(box[1]-box[0]+xyz[0])/2,(box[3]-box[2])/2,0])
    j = 0
    for i in range(num_atoms):
        pos = particle_list[j].pos
        
        pit = np.array([pit_D/2,pit_D/2,pit_H])
        pos = (pos-center)/pit
        dis = np.linalg.norm(pos)
        if(dis <= 1.0):
            particle_list.pop(j)
        else:
            if (pos[2] == 0):
                top_layer.append(j)
                particle_list[j].name='C'
            j = j+1

    num_atoms = len(particle_list)
    print("New number of atoms = ",num_atoms)
    for i in range(num_atoms):
        particle_list[i].idx = i+1
    

# Adding the adatom
if(interaction_model > 1):
    # Get the atom-indices for the adatom site
    site, site_idx = Get_site_idx(site, num_rep, particle_list, info_fname, interaction_model)
    pos = np.zeros(3)
    if(site == 'bridge'):
        pos = (particle_list[site_idx[0]-1].pos + particle_list[site_idx[1]-1].pos)/2
        pos[2]    += 1.0
        atom_type += 1
        particle_list.append(Particle(mass_adatom, particle_list[-1].idx+1, atom_type, adatom_name))
        particle_list[-1].pos = np.copy(pos)
        particle_list[-1].vel = np.zeros(3)
        #particle_list[-1].vel = np.copy(particle_list[-2].vel)
    elif(site == 'top'):
        pos = np.copy(particle_list[site_idx[0]-1].pos)
        pos[2]    += 1.3
        atom_type += 1
        particle_list.append(Particle(mass_adatom, particle_list[-1].idx+1, atom_type, adatom_name))
        particle_list[-1].pos = np.copy(pos)
        particle_list[-1].vel = np.zeros(3)
        
    elif(site == 'edge'):
        site_idx, zig_zag, arm_chair = Get_edge_site_idx(particle_list, info_fname, interaction_model,
                                                         top_layer, lat_a, pit_D, pit_H, center, edge_type)
        pos = np.copy(particle_list[site_idx[0]-1].pos)
        pos[2]    += 1.5
        atom_type += 1
        particle_list.append(Particle(mass_adatom, particle_list[-1].idx+1, atom_type, adatom_name))
        particle_list[-1].pos = np.copy(pos)
        particle_list[-1].vel = np.zeros(3)


## Zeroing out the com velocities
com_vel    = np.zeros(3)
total_mass = 0

for i in range(len(particle_list)):
    com_vel    += particle_list[i].vel * particle_list[i].mass
    total_mass += particle_list[i].mass

com_vel = com_vel/total_mass
print("com_vel before zeroing out : ", com_vel)
for i in range(len(particle_list)):
    particle_list[i].vel -= com_vel

# Write pos files
fname_lmp = Pos_prefix + ".lmp"
num_atoms = len(particle_list)
if(write_xyz):
    fname_xyz = Pos_prefix + ".xyz"
    
    fxyz = open(fname_xyz,"w")
    fxyz.write("%d \n"%(len(particle_list)))
    fxyz.write("Atoms\n")
    for i in range(num_atoms):
        fxyz.write("%s %6.8f %6.8f %6.8f\n"%(particle_list[i].name, particle_list[i].pos[0], particle_list[i].pos[1], particle_list[i].pos[2]))
    fxyz.close()

flmp = open(fname_lmp,"w")
flmp.write(" Hexagonal lattice structure of dimension %d x %d x %d with a = %.4f and c = %.4f\n"%(num_rep[0],num_rep[1],num_rep[2],lat_a,lat_c))
flmp.write("%d atoms\n"%(num_atoms))
flmp.write("%d atom types\n\n"%(atom_type))
flmp.write("%6.8E %6.8E xlo xhi\n"%(box[0],box[1]))
flmp.write("%6.8E %6.8E ylo yhi\n"%(box[2],box[3]))
flmp.write("%6.8E %6.8E zlo zhi\n\n"%(box[4],box[5]))
flmp.write("%6.8E  %6.8E  %6.8E  xy xz yz\n\n"%(xyz[0],xyz[1],xyz[2]))
flmp.write("Atoms\n\n")
for i in range(num_atoms):
    flmp.write("%5d  %2d  %6.4E  %6.8E  %6.8E  %6.8E\n"%(particle_list[i].idx,particle_list[i].type,0.0,particle_list[i].pos[0],particle_list[i].pos[1],particle_list[i].pos[2]))

flmp.write("\nVelocities\n\n")
for i in range(num_atoms):
    flmp.write("%5d  %6.8E  %6.8E  %6.8E\n"%(particle_list[i].idx,particle_list[i].vel[0],particle_list[i].vel[1],particle_list[i].vel[2]))

flmp.close()
