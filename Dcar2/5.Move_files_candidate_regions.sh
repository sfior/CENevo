#!/bin/bash

# Important note: It's important to run pooled_sumstats_calculator_OBS.py in separate directories (for each Chr) because the script outputs a temporary file "summary_stats_temp.txt" with a common name, and these can thus be mixed up!
# Thus, we move files into separate folders (based on chromosome)

for i in TCF1 FT; do
	mkdir ${i}
	mv *.${i}.* ${i}
	mv ${i}_filelist.txt ${i}
done

