#PBS -S /bin/bash
#PBS -N Gen_Lattice_1000
#PBS -o Gen_Lattice_1000.txt

#PBS -l select=4:ncpus=16:mpiprocs=16:model=san
#PBS -l walltime=04:00:00

#PBS -q normal
#PBS -W group_list=e1543

export Temp=1000
export Potential="COMB3"

export site="edge"  # Bridge, top, edge
export edge_type="zig-zag"  # zig-zag or arm-chair : only used when site==edge
export interaction_model=3     # 1 : Lattice does not contain adatom                                                                                  
                               # 2 : Lattice contains adatom. Bridge site is determined by mean of lattice positions       
                               # 3 : Lattice contains adatom. Bridge site is taken from Stage 1
cd $PBS_O_WORKDIR

module load comp-intel/2018.3.222
module load mpi-sgi/mpt

source ~/.bash_profile

if test -f "info.lattice" ; then
    echo "Removing info.lattice"
    rm info.lattice
fi


struc_file="gen_lattice.py"

if [ ${interaction_model} == 1 ]; then
    lammps_file="lammps_generate_lattice_1.in"
else
    lammps_file="lammps_generate_lattice_2.in"
fi

echo "Lammps file : "$lammps_file

if [ ${site} == "bridge" ] || [ ${site} == "top" ] ; then
    python $struc_file $Potential $site $Temp $interaction_model

elif [ ${site} == "edge" ]; then
    python $struc_file $Potential $site $edge_type $Temp $interaction_model 16.0 6.0
fi

export struc_file
export lammps_file

if test -f "info.dat" ; then
    echo "Removing info.dat"
    rm info.dat
fi


mpiexec -np 64 $lammps_mpi -in $lammps_file

echo $site >> info.dat
echo $Temp >> info.dat
echo $lammps_file >> info.dat
