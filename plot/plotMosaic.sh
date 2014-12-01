#!/bin/bash
HDF5FILE=data/SFR.h5
sort=("MorphType" "McorSD" "u-r" "Mr" "Mcor")
for s in $sort
do  
    i=0
    while [ $i != 40 ]
    do
        echo "SORT: ${s} iT: ${i}"
        ./plotMosaic.py ${HDF5FILE} ${i} ${s}
       let "i = i + 1"
   done
done
