#

Cur_Dir=$(pwd)
Dir="pot_files/temp-2001/"
cd $Dir

find ./ -name "*.dat" > $Cur_Dir/pot_fname.dump

begin=$(sed -n 1p $Cur_Dir/info.dat)
end=$(sed -n 2p $Cur_Dir/info.dat)
incr=$(sed -n 3p $Cur_Dir/info.dat)
n_at=$(sed -n 4p $Cur_Dir/info.dat)

beg=0
end=$(echo "$end-$begin" | bc)
num=$(echo "($end-$beg)/$incr" | bc)

echo ${beg} ${incr} ${end} $num

pot_num=$(find ./ -name "*.dat" | wc -l)

if [ $num -ne $pot_num ]
then
    echo " ****** Number of potential files do not match *****"
    echo " Expected number = "$num
    echo " Actual   number = "$pot_num

    #return 1
fi

while read -r line
do
    fnum=$(echo "$line" | tr -dc '0-9')
    #fnum=$(echo "$line" | grep -Eo '[0-9]+$')
    echo $fnum >> $Cur_Dir/pot_fnum.dump

done < $Cur_Dir/pot_fname.dump

cd $Cur_Dir
