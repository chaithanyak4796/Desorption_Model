#!/bin/bash 
source /usr/share/modules/init/bash

module load comp-intel/2018.3.222
module load python3/Intel_Python_3.6_2018.3.222
source ~/.bash_profile

start_time=`date +%s`

begin=$1         # Begining index of the subjob array
end=$2           # Final    index of the subjob array
pot_dir=$3       # Directory where the potential files are stored
chunk_idx=$4     # Index of the chunk in the subjob array [0-nchunks-1]
nchunks=$5       # Number of chunks
c_beg=$6         # Begining index of the chunk in the subjob array
c_end=$7         # Final    index of the chunk in the subjob array

num=$(echo "($c_end-$c_beg+1)" | bc)
#### ------------ Verifying the number of potenial files ------------#
pot_num=$(find $pot_dir/ -name "PE*.dat" | wc -l)
echo $num $pot_num
if [ $num -ne $pot_num ]; then
    echo " ****** Number of potential files do not match *****"
    echo " Expected number = "$num
    echo " Actual   number = "$pot_num
    echo " ***************************************************"
    exit 
fi
##_______________________________________________________________________________________________________________________##  

Error_file="Error.$begin"
#### ------------- Verifying that all the subjobs ran properly ----------#
if test -f "$Error_file"; then
    perfect=1
    echo "Testing if this job array has any error files."
    declare $(awk -v awkbeg="$begin" '{if ($1 == awkbeg) {print "perfect=0"} }' "${Error_file}")

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
m_beg=$c_beg
m_end=$c_end
merged_pot="./pot_files/Merged_Pot_"$begin".dat."$chunk_idx

file=$pot_dir"/PE_"$(printf %06g $m_beg)".dat"
echo "First file : " $file
tail -n+3 $file > $merged_pot
rm $file

for (( i=$m_beg+1; i<=$m_end; i++ ))
do
    file=$pot_dir"/PE_"$(printf %06g $i)".dat"
    #echo $file
    tail -n+3 $file >> $merged_pot
    rm $file
done

echo "Final file : "$file

echo " "
echo "Done merging potential files for this sub-job and chunk_idx $chunk_idx"

if [ $chunk_idx -eq $((nchunks-1)) ]; then
    echo "Merging all the chunks in this subjob array"
    merged_pot="./pot_files/Merged_Pot_"$begin".dat"
    
    file="./pot_files/Merged_Pot_"$begin".dat.0"
    echo "First file : " $file
    cat $file > $merged_pot
    rm $file
    
    for (( i=1; i<$nchunks; i++ ))
    do
	file="./pot_files/Merged_Pot_"$begin".dat."$i
	echo $file
	cat $file >> $merged_pot
	rm $file
    done
    echo "Done merging all the chunks in this subjob array"
fi


end_time=`date +%s`
runtime=$((end_time-start_time))
echo "RunTime for merging= "${runtime}"s"

