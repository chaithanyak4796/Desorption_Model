#PBS -S /bin/bash
#PBS -N Merge_all
#PBS -o Merge_all.txt

#PBS -l select=1:ncpus=1:model=san
#PBS -l walltime=01:00:00

#PBS -q normal
#PBS -W group_list=e1543

cd $PBS_O_WORKDIR

source ~/.bash_profile

start=`date +%s`

Pos_dir="./pos_files"
Pot_dir="./pot_files"

#----------------- Reading the necessary details for output file names ------------------#

site=$(sed -n 5p info.dat)
if [ ${site} == "bridge" ] || [ ${site} == "top" ] ; then
    struc_file="genC_lammps_surface_charge.cc"
    lammps_file="lammps_generate_lattice_bridge.in"
elif [ ${site} == "edge" ]; then
    struc_file="genC_lammps_etch_pit.cc"
    lammps_file="lammps_generate_lattice_etch.in"
fi

dt=$(grep "dt equal" ${lammps_file} | grep -Eo '[0-9]([.][0-9]+)?')
unit=$(grep "units" ${lammps_file} | awk '{print$ 2}')
if [[ unit -eq "metal" ]]; then
    dt=$(bc -l <<< "${dt}*1000" )
    dt=$(printf "%.1f\n" $dt)
else
    dt=$(printf "%.1f\n" $dt)
fi

Temp=$(sed -n 6p info.dat)

echo "units    = "$unit
echo "dt [fs]  = "$dt
echo "Temp [K] = "$Temp

#-________________________ Merging files _______________________________#

readarray beg_idx < info.jobs
echo "Original job_idx:" ${beg_idx[@]}

cat << EOF > sort.py
#!/usr/bin/env python3
import numpy as np
data = np.loadtxt("info.jobs")
data = np.sort(data)
np.savetxt("info.jobs.sorted",data,fmt="%d")
EOF

python sort.py

rm sort.py

readarray sorted < info.jobs.sorted
echo "Sorted job_idx:  " ${sorted[@]}

Final_merged="./Merged_Pot_${dt}fs_${Temp}K.dat"

if test -f "$Final_merged"; then
    echo "File exists. Deleting it first"
    rm $Final_merged
fi

for num in "${sorted[@]}"
do
    i=$(echo $num | tr -dc '0-9')
    file="./pot_files/Merged_Pot_"$i".dat"
    echo "File : $file"
    cat $file >> $Final_merged
    
done

echo "Done merging all files"

#________________________________________________________________________#

echo " "
echo "Computing the average and fluctuating potentials "
python Compute_pot.py $dt $Temp $Final_merged
echo " "

#________________________________________________________________________#

echo "  "
echo "Deleting temporary directories"
rm -r $Pos_dir $Pot_dir
echo "Done"
echo " "

end=`date +%s`
runtime=$((end-start))
echo "Total RunTime = "${runtime}"s"
