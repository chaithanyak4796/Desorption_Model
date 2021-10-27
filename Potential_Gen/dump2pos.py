import numpy as np
import sys

if (len(sys.argv) != 4):
    sys.exit("!!! Input file, output file and z_shift needs to be given as input !!!!!!!")
else:
    In_file  = sys.argv[1]
    Out_file = sys.argv[2]
    z_shift = float(sys.argv[3])
    
info_fname = 'info.lattice'
fr = open(info_fname,"r")
line = fr.readline().strip("\n")
interaction_model = int(line)
print("Interaction model = ",interaction_model)
line = fr.readline().strip("\n")
bond_style = line
if(bond_style == "Bridge" or "bridge"):
    bond_style = "bridge"
    C_idx = np.array([0,0])
    line = fr.readline().strip("\n")
    C_idx[0] = float(line)
    line = fr.readline().strip("\n")
    C_idx[1] = float(line)
elif(bond_style == "Top" or "top"):
    bond_style = "top"
    C_idx = np.array([0])
    line = fr.readline().strip("\n")
    C_idx[0] = float(line)
fr.close()

#bond_style = "bridge"  # bridge, top
#C_idx  = np.array([398,559])
z_max  = 5.00   # Angstrom

print("\n\n !!!!!!!!! WARNING : Assuming the modify ... sort id was turned on in LAMMPS script !!!!!!!!!\n\n")

print("dump_file : ",In_file)
print("pos_file  : ",Out_file)

fr = open(In_file,"r")
fw = open(Out_file,"w")

line = fr.readline()   # ITEM: TIMESTEP
line = fr.readline().split()   # timestep
t = int(line[0])
print("Timestep = ",t)

line = fr.readline()  # ITEM
line = fr.readline().split() # Number of atoms
natoms = int(line[0])

print("Original # of atoms = ",natoms)

xyz = np.zeros(3)
line = fr.readline() # ITEM
line = fr.readline().split()  # box_x
xlo = float(line[0])
xhi = float(line[1])
xyz[0] = float(line[2])

line = fr.readline().split()  # box_y
ylo = float(line[0])
yhi = float(line[1])
xyz[1] = float(line[2])

line = fr.readline().split()  # box_z
zlo = float(line[0])
zhi = float(line[1])
xyz[2] = float(line[2])

print("Box dimensions x : %6.8E %6.8E"%(xlo,xhi))
print("Box dimensions y : %6.8E %6.8E"%(ylo,yhi))
print("Box dimensions z : %6.8E %6.8E"%(zlo,zhi))

fr.close()

atom_data = np.loadtxt(In_file,skiprows=9)

atom_types = np.unique(atom_data[:,1])
n_type = len(atom_types)

print("Original # of atom types = ",n_type)

################################################################
atom_data[:,1] = 1
atom_data[:,5] -= z_shift
#################################################################

if(bond_style == "top"):
    if(len(C_idx) != 1):
        sys.exit("!!! Incorrect # of C atoms for top site config. !!!!!!")
    else:
        C_idx = C_idx - 1 # Converting to python idx style
        atom_data[C_idx[0],1] = 2
        idx = C_idx[0]
        O_pos = np.array([ atom_data[idx][3], atom_data[idx][4], z_max ])
elif(bond_style == "bridge"):
    if(len(C_idx) != 2):
        sys.exit("!!! Incorrect # of C atoms for bridge site config. !!!!!!")
    else:
        C_idx = C_idx - 1 # Converting to python idx style
        idx1 = C_idx[0]
        idx2 = C_idx[1]
        atom_data[idx1,1] = 2
        atom_data[idx2,1] = 2
        
        if(interaction_model == 1 or interaction_model == 2):
            O_pos = np.array([ (atom_data[idx1][3]+atom_data[idx2][3])/2 , (atom_data[idx1][4]+atom_data[idx2][4])/2 , z_max ])
        else:
            O_pos = np.array([ atom_data[natoms-1][3] , atom_data[natoms-1][4] , z_max ])
else:
    sys.exit("!!! bond_style unknown. Accepted values are 'top' and 'bridge' !!!!!!")

print("\nPosition of O atom : %6.8E %6.8E %6.8E"%(O_pos[0],O_pos[1],O_pos[2]))

if(interaction_model == 1):
    O_data = np.array([natoms+1,3,0,O_pos[0],O_pos[1],O_pos[2]])
    atom_data = np.vstack([atom_data,O_data])
else:
    O_data = np.array([natoms,3,0,O_pos[0],O_pos[1],O_pos[2]])
    atom_data[-1] = O_data

natoms = len(atom_data[:,0])
n_types = len(np.unique(atom_data[:,1]))

print("Final # of atoms = ",natoms)
print("Final # of atom types = ",n_types)

print("\nWriting the updated position file")

fw.write("# Displaced lattice with adatom. Timestep = %d\n" %(t) )
fw.write("   %d atoms \n"%(natoms) )
fw.write("   %d atom types\n" %(n_types) )
fw.write("\n ")

fw.write("%6.8E  %6.8E  xlo xhi\n" %(xlo,xhi) )
fw.write("%6.8E  %6.8E  ylo yhi\n" %(ylo,yhi) )
fw.write("%6.8E  %6.8E  zlo zhi\n" %(zlo,zhi) )
fw.write("\n")
fw.write("%6.8E  %6.8E  %6.8E xy xz yz\n\n"%(xyz[0], xyz[1], xyz[2]))
fw.write("Atoms\n\n")

for i in range(natoms):
    fw.write( "%04d %d 0.000000  %4.6E  %4.6E  %4.6E\n"%(atom_data[i][0],atom_data[i][1],atom_data[i][3],atom_data[i][4],atom_data[i][5]) )

    
print("Done writing \n")

fw.close()
