#!/bin/bash 
source /usr/share/modules/init/bash

start=`date +%s`

icase=$1
z_shift=$2
n_atoms=$3

nprocs=4

source ~/.bash_profile
module load comp-intel/2018.3.222
module load mpi-sgi/mpt
module load python3/Intel_Python_3.6_2018.3.222

#echo "icase     = "$icase
#echo "z_shift   = "$z_shift
#echo "num_atoms = "$n_atoms

pos_fname="./pos_files/coords$icase.pos"
lat_fname="./pos_files/lattice$icase.pos"
pot_fname="./pot_files_ind/PE_$icase.dat"
lammps_in="./lammps_potential.in"

#need to create PBS_NODEFILE to run mpiexec
echo `hostname` >> nodefile_$icase
echo `hostname` >> nodefile_$icase
echo `hostname` >> nodefile_$icase
echo `hostname` >> nodefile_$icase
echo `hostname` >> nodefile_$icase
export PBS_NODEFILE=`pwd`/nodefile_$icase

# need to disable pinning to avoid having multiple processes run
# on the same set of CPUs
export MPI_DSM_DISTRIBUTE=0

l1=$(echo "$icase*($n_atoms+9) +1" | bc)
l2=$(echo "($icase+1)*($n_atoms+9)" | bc)

sed -n "${l1},${l2}p" coords.pos > $pos_fname

temp_out_file="./output_${icase}.out"
python dump2pos.py $pos_fname $lat_fname $z_shift > $temp_out_file

export lat_fname
export pot_fname

mpiexec -np $nprocs $lammps_mpi -in $lammps_in >> $temp_out_file

#rm $temp_out_file $pos_fname $lat_fname

end=`date +%s`
runtime=$((end-start))

echo "icase = ${icase},  time = ${runtime}"

# remove nodefile at end of run
rm nodefile_$icase
