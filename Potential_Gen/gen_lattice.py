import numpy as np
import matplotlib.pyplot as plt
from gen_lattice_functions import *
import sys

if(len(sys.argv) == 5):
    Pot  = str(sys.argv[1])
    site = str(sys.argv[2])
    Temp = float(sys.argv[3])
    interaction_model = int(sys.argv[4])
elif(len(sys.argv) == 1):
    Pot  = 'COMB3'
    site = "Bridge"
    Temp = 1000
    interaction_model = 2
else:
    print("Error: INcorrect number of arguments.")
    sys.exit(0);

num_rep     = [10,10,2]
mark_layers = True   # True : Each graphene layer is assigned a different atom_type
write_xyz   = True
adatom_name = 'O'

vac_height  = 20
Lat_const_fname = './Lat_const.dat'
Pos_prefix      = "graphite_slab"
info_fname      = './info.lattice'

pot_name, Lat_const = Read_Lat_Const(Lat_const_fname)

idx_pot = np.where(pot_name == Pot)[0]
if(len(idx_pot) == 0):
    print(" Error : Potential not found.")
    sys.exit()
idx_pot = idx_pot[0]
idx_Temp = np.where(Lat_const[idx_pot][:,0] == Temp)[0]
if(len(idx_Temp) == 0):
    print(" Error : Temp not found in the potential.")
    sys.exit() 
idx_Temp = idx_Temp[0]

lat_param = Lat_const[idx_pot][idx_Temp]
print("Lattice params found = ",lat_param)

lat_a = lat_param[1]*3**0.5
lat_c = lat_param[2]*2

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

num_atoms = 4 * np.product(num_rep)
particle_list = []
pos = np.zeros((3,))
base_lattice = np.zeros((3,))

idx = 1
for i in range(num_rep[0]):
    for j in range(num_rep[1]):
        for k in range(num_rep[2]):
            base_lattice[0] = i*lattice_vec[0][0] + j*lattice_vec[1][0] + k*lattice_vec[2][0]
            base_lattice[1] = i*lattice_vec[0][1] + j*lattice_vec[1][1] + k*lattice_vec[2][1]
            base_lattice[2] = i*lattice_vec[0][2] + j*lattice_vec[1][2] + k*lattice_vec[2][2]
            
            pos = np.copy(base_lattice)
            atom_type = 1
            if(mark_layers): atom_type = 2*k + 1
            particle_list.append(Particle(idx, atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            
            pos[0] = base_lattice[0] + 0.5*lat_a
            pos[1] = base_lattice[1] + (0.5/3**0.5) * lat_a
            pos[2] = base_lattice[2] + 0
            particle_list.append(Particle(idx+1, atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            
            pos[0] = base_lattice[0] + 0.5*lat_a
            pos[1] = base_lattice[1] + (0.5/3.0**0.5) * lat_a
            pos[2] = base_lattice[2] + 0.5*lat_c
            if(mark_layers): atom_type = 2*k + 2
            particle_list.append(Particle(idx+2, atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            
            pos[0] = base_lattice[0] + 1.0*lat_a
            pos[1] = base_lattice[1] + (1.0/3.0**0.5) * lat_a
            pos[2] = base_lattice[2] + 0.5*lat_c
            particle_list.append(Particle(idx+3,atom_type,'C'))
            particle_list[-1].pos = np.copy(pos)
            z_max = np.copy(pos[2])
            idx += 4

num_atoms = len(particle_list)
# Adjust the z coordinates
box[4] -= z_max
box[5] -= z_max
for i in range(num_atoms):
    particle_list[i].pos[2] -= z_max

# Get the atom-indices for the adatom site
site, site_idx = Get_site_idx(site, num_rep, particle_list, info_fname, interaction_model)

# Adding the adatom
if(interaction_model > 1):
    pos = np.zeros(3)
    if(site == 'bridge'):
        pos = (particle_list[site_idx[0]-1].pos + particle_list[site_idx[1]-1].pos)/2
        pos[2]    += 1.0
        atom_type += 1
        particle_list.append(Particle(idx,atom_type,adatom_name))
        particle_list[-1].pos = np.copy(pos)
    elif(site == 'top'):
        pos = (particle_list[site_idx[0]-1].pos)
        pos[2]    += 1.0
        atom_type += 1
        particle_list.append(Particle(idx,atom_type,adatom_name))
        particle_list[-1].pos = np.copy(pos)

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
flmp.close()

