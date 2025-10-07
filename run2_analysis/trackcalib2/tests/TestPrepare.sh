#!/bin/bash


#Make sure you are in the right directory
cd /afs/cern.ch/work/r/rekopecn/TrackCalib2/trackcalib2


 
# Declare a string array with type
declare -a MethodArray=("Long")
declare -a YearArray=("2012_50ns_strip" "2015_25ns" "2017_25ns") 
declare -a ModeArray=("MC" "Data")


# Loop over year, mode and method
for year in ${YearArray[@]}; do
    for method in ${MethodArray[@]}; do
    	for mode in "${ModeArray[@]}"; do
  	    python3.9 trackcalib.py prepare -mode "$mode" -year "$year" -method "$method" -max_entries 1000 --force -test -high_multi
    	done   	
    done
done 

