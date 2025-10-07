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
			conda run -n trackcalib2 python trackcalib.py prepare -mode "$mode" -year "$year" -method "$method" -v -max_entries 1000 2>&1 | tee -a results/logs/logs_dask/prepare_"$mode"_"$year"_"$method".log
		done   
	done
done
