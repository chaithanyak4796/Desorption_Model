#!/bin/bash

if test -f "info.jobs"; then
    echo "Deleting info.jobs"
    rm info.jobs
fi

if test -f "info.jobs.sorted"; then
    echo "Deleting info.jobs.sorted"
    rm info.jobs.sorted
fi

if test -f "Error.files"; then
    echo "Deleting Error.files"
    rm Error.files
fi

job0=$(qsub Lattice_Gen.sh)
echo $job0

job1=$(qsub -W depend=afterok:$job0 Extract_Pot.sh)
#job1=$(qsub Extract_Pot.sh)
echo $job1 

job2=$(qsub -W depend=afterok:$job1 Merge_all.sh)
echo $job2
