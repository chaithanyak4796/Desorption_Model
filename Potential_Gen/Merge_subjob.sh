#!/bin/bash 
source /usr/share/modules/init/bash

module load comp-intel/2018.3.222
#module load mpi-sgi/mpt
module load python3/Intel_Python_3.6_2018.3.222
source ~/.bash_profile

start=`date +%s`

begin=$1
end=$2
pot_dir=$3

num=$(echo "($end-$begin+1)" | bc)

#### ------------ Verifying the number of potenial files ------------#
pot_num=$(find $pot_dir/ -name "*.dat" | wc -l)
if [ $num -ne $pot_num ]
then
    echo " ****** Number of potential files do not match *****"
    echo " Expected number = "$num
    echo " Actual   number = "$pot_num
    echo " ***************************************************"
    exit 
fi
##_______________________________________________________________________________________________________________________##  

#### ------------- Verifying that all the subjobs ran properly ----------#
if test -f "Error.files"; then
    perfect=1
    echo "Testing if this job array has any error files."
    declare $(awk -v awkbeg="$begin" '{if ($1 == awkbeg) {print "perfect=0"} }' "Error.files")

    if (( $perfect == 0 )); then
        echo  "This job array contains  atleast one error file."
        echo  "Exiting without merging"
	exit 
    fi
    if (( $perfect == 1 )); then
        echo " This job array ran successfully."
    fi
fi
#_______________________________________________________________________________________________________________________##

echo " Merging the potential files"
beg=$begin
merged_pot="./pot_files/Merged_Pot_"$begin".dat"
file=$pot_dir"/PE_"$(printf %06g $beg)".dat"
echo "First file : " $file
tail -n+3 $file > $merged_pot

for (( i=$beg+1; i<=$end; i++ ))
do
    file=$pot_dir"/PE_"$(printf %06g $i)".dat"
    echo $file
    tail -n+3 $file >> $merged_pot
    rm $file
done

echo "Final file : "$file

echo " "
echo "Done merging potential files for this sub-job"


end=`date +%s`
runtime=$((end-start))
echo "RunTime for merging= "${runtime}"s"

