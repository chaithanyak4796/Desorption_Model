#!/bin/bash

begin=15001

if test -f "Error.files"; then
    perfect=1
    echo "begin = " $begin
    echo "Testing if this job array has any error files."
#    awk -v awkvar="$begin" '{print $1 " " awkvar}' "Error.files"
    declare $(awk -v awkbeg="$begin" '{if ($1 == awkbeg) {print "perfect=0"} }' "Error.files")
   
    if (( $perfect == 0 )); then
	echo  "This job array contains  atleast one error file."
	echo  "Exiting without merging"
       
    fi
    if (( $perfect == 1 )); then
	echo " This job array ran successfully."
    fi

    echo "********* Hello ******************"
fi
