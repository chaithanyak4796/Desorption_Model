#!/bin/bash

exp="Normal termination"
file="output_000233.out"
act=$(tail -n2 $file | head -n1)

if [[ "$act" != "Normal termination" ]]; then
    echo "Strings are not equal"
else
    echo "Strings are equal"
fi

num_lines=478
pot_fname="../pot_files/temp-1/PE_000001.dat"
perfect=1
num_lines_act=$(wc -l < $pot_fname)
if [ $num_lines_act != $num_lines ]; then  # Then check if the number of lines in the potential file matches the expected value
    echo "------ Incorrect number of lines in $pot_fname"
    echo "expected number = $num_lines"
    echo "actual number   = $num_lines_act"
    perfect=0
fi

echo "perfect = $perfect"
