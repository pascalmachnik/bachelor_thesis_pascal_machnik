#!/usr/bin/bash
set -euo pipefail

# Declare a string array with type
declare -a MethodArray=("Velo" "T" "Long")
declare -a YearArray=("2011_50ns_strip" "2012_50ns_strip" "2016_25ns_strip" "2017_25ns_strip" "2015_50ns" "2015_25ns" "2016_25ns" "2017_25ns" "2018_25ns") 
declare -a ModeArray=("Data" "MC")

# Loop over year, mode and method
for year in ${YearArray[@]}; do
    for method in ${MethodArray[@]}; do
    	for mode in "${ModeArray[@]}"; do
  	    conda run -n trackcalib2 python trackcalib.py fit -mode "$mode" -year "$year" -method "$method" -sim_fit 
    	done   	
    done
done 

