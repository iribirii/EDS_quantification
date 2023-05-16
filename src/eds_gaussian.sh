#!/bin/bash

##############################################################################
#Script Name: eds_gaussian.sh
#Description: This script calculates the Electron Density Shift (EDS) by 
#             substracting the electron density of each individual fragment 
#             to the electron density of the whole system
#
#             EDS = rho_sys - rho_frag1 - rho_frag2 - [...]
#
#             Files must be named <filename>_sys.chk, <filename>_frag1.chk, ...
#Version: 1.0
#Authors: Inigo Iribarren, Goar Sanchez, Ibon Alkorta 
#Date: 2020-11-06
#Usage: eds_gaussian.sh <filename> <n_points>
##############################################################################

GAUSS=/path/to/gaussian
export GAUSS

# Convert CHK files tot FCHK
name=$1
n_points=$2
all_files=`ls -1 $1*.chk | awk -F.chk '{print $1}'`

for i in $all_files
do
    echo 'Converting ' $i'.chk into '$i'.fchk'
    $GAUSS/formchk $i.chk $i.fchk
done

echo '------------------------------------------------------------------------'
echo 'All files transformed to FCHK'
echo '------------------------------------------------------------------------'

# Calculating EDS 
# FCHK -> cub

files_fchk=`ls -1 ${name}_sys.fchk | awk -F_sys.fchk '{print 1}'`

for i in $files_fchk
do
    echo 'Working with '$i
    # Check if MP2 density or SCF
    mp2=$(grep MP2 ${i}_sys.fchk | wc -l)
    if [ $mp2 != 0 ]
    then
        density='MP2'
    else
        density='SCF'
    fi

    # Convert SYS FCHKfile into cub and taking the options for the fragX files
    echo 'Converting '${i}'_sys.fchk into '${i}'_sys.cub'
    $GAUSS/cubegen 0 density=$density ${i}_sys.fchk ${i}_sys.cub $n_points h

    head -n4 ${i}_sys.cub | tail -n2 | awk 'print "-"$1"   "$2"   "$3"   "$4' > options.dat
    head -n6 ${i}_sys.cub | tail -n2 | awk 'print $1"   "$2"   "$3"   "$4' >> options.dat

    # Convert each FRAG FCHK to cub
    for file_fchk in `ls ${i}_frag*fchk`
    do
        name=$(echo $file_fchk | cut -d"." -f1)
        echo 'Converting '${name}'.fchk to '${name}'.cub'
        $GAUSS/cubegen 0 density=$density ${name}.fchk ${name}.cub -1 h < options.dat
    done

    # Get number of fragments 
    n_frag=$(ls ${i}_frag*cub | wc -l)

    # SYS - frag1
    echo 'SU'              > temp_1.temp
    echo ${i}_sys.cub     >> temp_1.temp
    echo 'yes'            >> temp_1.temp
    echo ${i}_frag1.cub   >> temp_1.temp
    echo 'yes'            >> temp_1.temp
    echo ${i}_inter_1.cub >> temp_1.temp
    echo 'yes'            >> temp_1.temp
    $GAUSS/cubman < temp_1.temp

    # Intermediate EDSs - fragX
    for frag in `seq 1 $n_frag`
    do
        let pre_frag=${frag}
        let actual_frag=${frag}+1
        echo 'SU'                               > temp_${actual_frag}.temp
        echo ${i}_inter_${pre_frag}.cub        >> temp_${actual_frag}.temp
        echo 'yes'                             >> temp_${actual_frag}.temp
        echo ${i}_frag${actual_frag}.cub       >> temp_${actual_frag}.temp
        echo 'yes'                             >> temp_${actual_frag}.temp
        if [[ $actual_frag == $n_frag ]]
        then
            echo ${i}_final.cub                >> temp_${actual_frag}.temp   
        else
            echo ${i}_inter_${actual_frag}.cub >> temp_${actual_frag}.temp   
        fi
        echo 'yes'                             >> temp_${actual_frag}.temp   

        $GAUSS/cubman < temp_${pre_frag}.temp
    done

    rm -f *temp
done





