#PBS -S /bin/bash
#PBS -N Gen_Lattice_1000
#PBS -o Gen_Lattice_1000.txt

#PBS -l select=4:ncpus=16:mpiprocs=16:model=san
#PBS -l walltime=02:00:00

#PBS -q normal
#PBS -W group_list=e1543

export site="bridge"  # Bridge, top, edge
export Temp=1000
export Potential="COMB3"
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

if [ ${site} == "bridge" ] || [ ${site} == "top" ] ; then
    lammps_file="lammps_generate_lattice_bridge.in"
    struc_file="gen_lattice.py"
    python $struc_file $Potential $site $Temp $interaction_model

elif [ ${site} == "edge" ]; then
    struc_file="genC_lammps_etch_pit.cc"
    lammps_file="lammps_generate_lattice_etch.in"
    g++ -std=c++11 ${struc_file}
    ./a.out 16 6
fi

export struc_file
export lammps_file

mkdir -p  pos_files pot_files

if test -f "info.dat" ; then
    echo "Removing info.dat"
    rm info.dat
fi


mpiexec -np 64 $lammps_mpi -in $lammps_file

echo $site >> info.dat
echo $Temp >> info.dat
