#PBS -S /bin/bash

#PBS -l select=1:ncpus=4:model=san
#PBS -l walltime=00:30:00

#PBS -q devel
#PBS -W group_list=e1543

cd $PBS_O_WORKDIR

#module load comp-intel/2018.3.222
#module load mpi-sgi/mpt
#module load python3/Intel_Python_3.6_2018.3.222

start=`date +%s`


Dump_dir="./dump_coords"
Pos_dir="./pos_files"
Pot_dir="./pot_files_ind"

export Pos_dir
export Pot_dir

##_______________________________________________________________________________________________________________________##

echo " "
echo " Computing the average z_com using z_surf.dat"
z_shift=$( awk  ' NR>=3 { total+=$4; count++ } END {print total/count}' z_surf.dat)
echo "z_shift = $z_shift "



##_______________________________________________________________________________________________________________________##

echo " "
echo " Running LAMMPS to generate potential energy curves for each displaced lattice"
#echo " !!!! This is currently running in serial. Needs to be parrallelized !!!!!"

#N_proc=16


begin=$(sed -n 1p info.dat)
end=$(sed -n 2p info.dat)
incr=$(sed -n 3p info.dat)
n_at=$(sed -n 4p info.dat)

#beg=0
#end=$(echo "$end-$begin" | bc)
#num=$(echo "($end-$beg)/$incr" | bc)

arr=( 1366, 3039, 11030, 17033, 29005, 38945, 39005, 45002)
arr=( 1369, 3037, 11031, 17034, 29004, 38944, 39006, 45001)
#arr=(1, 2, 3)
for i in "${arr[@]}"
do
    printf -v j "%06d" $i
    ./RunLAMMPS_Parallel_ind.sh $j $z_shift $n_at
done 
#seq -f "%06g"  ${begin} ${end} | parallel -j 4 --sshloginfile "$PBS_NODEFILE" "cd $PWD; ./RunLAMMPS_Parallel.sh {} $z_shift $n_at"

#seq -f "%06g" ${beg} ${num}  | parallel -j 10 "cd $PWD; ./RunLAMMPS_Parallel.sh {} $z_shift $n_at"

echo " "
echo " Done collecting the potential files"

##_______________________________________________________________________________________________________________________##  

end=`date +%s`
runtime=$((end-start))
echo "Total RunTime = "${runtime}"s"

