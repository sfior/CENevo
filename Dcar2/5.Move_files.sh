#!/bin/bash

# Important note: It's important to run pooled_sumstats_calculator_OBS.py in separate directories (for each Chr) because the script outputs a temporary file "summary_stats_temp.txt" with a common name, and these can thus be mixed up!
# Thus, we move files into separate folders (based on chromosome)

for i in `seq 1 15`; do
	mkdir Chr_${i}
	mv *.Chr${i}.* Chr_${i}
	mv Chr${i}_filelist.txt Chr_${i}
done


