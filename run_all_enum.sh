#!/bin/bash

mkdir -p out
touch "out\enum_results.csv"

conda.bat activate env-redist

echo "first_player,num_rows,num_cols,num_districts,grid_instance,definer_util" >> "out/enum_results.csv"

for first_player in R D
do
	for grid_instance in grid_instances/4x6*
	do
		# Four districts
		python definecombine/enum_dcp.py 6 4 4 $grid_instance $first_player >> "out/enum_results.csv"
		# Six districts
		python definecombine/enum_dcp.py 6 4 6 $grid_instance $first_player >> "out/enum_results.csv"
	done
	
	for grid_instance in grid_instances/6x6*
	do
		# Three districts
		python definecombine/enum_dcp.py 6 6 3 $grid_instance $first_player >> "out/enum_results.csv"
		# Six districts
		python definecombine/enum_dcp.py 6 6 6 $grid_instance $first_player >> "out/enum_results.csv"
	done
done

echo "Type anything and press Enter to exit."
read tempvar
