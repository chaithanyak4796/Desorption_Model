#!/bin/bash 
source /usr/share/modules/init/bash

start_time=`date +%s`

icase=$1
z_shift=$2
n_atoms=$3
beg_idx=$4
split_fname=$5
pot_dir=$6

nprocs=4

pos_fname="./pos_files/coords$icase.pos"
lat_fname="./pos_files/lattice$icase.pos"
pot_fname=$pot_dir"/PE_$icase.dat"
lammps_in="./lammps_potential.in"

#num_lines=478  # Modify this if z array params are changed
num_lines=$(( $(grep "run" ${lammps_in} | awk '{print $2}') + 3 ))

Error_file="Error.$beg_idx"

source ~/.bash_profile
module load comp-intel/2018.3.222
module load mpi-sgi/mpt
module load python3/Intel_Python_3.6_2018.3.222

#need to create PBS_NODEFILE to run mpiexec
for i in $(eval echo "{1..$nprocs}")
do 
    echo `hostname` >> nodefile_$icase
done
export PBS_NODEFILE=`pwd`/nodefile_$icase

# need to disable pinning to avoid having multiple processes run
# on the same set of CPUs
export MPI_DSM_DISTRIBUTE=0

l1=$(echo "($icase-$beg_idx)*($n_atoms+9) +1" | bc)
l2=$(echo "($icase-$beg_idx+1)*($n_atoms+9)" | bc)

sed -n "${l1},${l2}p" $split_fname > $pos_fname

temp_out_file="./output_${icase}.out"
python dump2pos.py $pos_fname $lat_fname $z_shift > $temp_out_file

export lat_fname
export pot_fname

mpiexec -np $nprocs $lammps_mpi -in $lammps_in >> $temp_out_file

######### Checking if the LAMMPS script ran successfully ############
perfect=1
if test -f "$pot_fname"; then      # First check if file exists
    message=$(tail -n2 $temp_out_file | head -n1)
    if [[ "$message" != "Normal termination"  ]]; then  # Check for termination
 	echo "--------- Termination error in $pot_fname "
	perfect=0
    fi
    num_lines_act=$(wc -l < $pot_fname)
    if [ $num_lines != $num_lines_act ]; then  # Then check if the number of lines in the potential file matches the expected value
    echo "------ Incorrect number of lines in $pot_fname"
    echo "expected number = $num_lines"
    echo "actual number   = $num_lines_act"
    perfect=0
    fi
    
else   # If file does not exist, just print that info to the error file
    echo "$pot_fname does not exist"
    echo $beg_idx $icase >> $Error_file
fi

if (( $perfect == 0)); then
    echo "$beg_idx $icase " >> $Error_file
fi

if (( $perfect == 1 )); then
    rm  $pos_fname $lat_fname $temp_out_file
fi

end_time=`date +%s`
runtime=$((end_time-start_time))

echo "icase = ${icase},  time = ${runtime}"

# remove nodefile at end of run
rm nodefile_$icase
