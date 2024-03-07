#!/bin/bash

# This script removes observations (windows) with 0 segsites. We remove these observation because they are inherently uninformative and will lead to misinformative inferences in downstream ABCestimation.

for i in $(seq 1 16); do 

	# Change into Chr directory
	cd Chr_${i}

	# Define input variables
	obs_file=$(ls *windowsize5000_M2_Q0_1kbOverlapWindows.obs)
	scaflist_file=$(ls *windowsize5000_M2_Q0_1kbOverlapWindows.scaffoldlist)
	prefix=$(basename "${obs_file}" .obs)

	# We note down all lines with S = 0.0
	# We want to remove these because these observations (which are uninformative) will lead to misinformative inference in ABCestimation.
	# See: https://stackoverflow.com/questions/17842903/delete-specific-rows-based-on-specific-word-in-column
	awk '{print $1}' ${obs_file} | grep -n "\b0.0\b" | awk -F ":" '{print $1}' > zerosegsites_rows_to_remove_wHeader.list

	# The obs file has +1 rows (header line) compared to the scaffoldlist file. 
	cat zerosegsites_rows_to_remove_wHeader.list | awk '{print $1 - 1 }' >  zerosegsites_rows_to_remove_noHeader.list

	# Condition removal on presence of 0 segsite windows.
	num_zerosegsitewindows=$(cat zerosegsites_rows_to_remove_wHeader.list | wc -l)
	if [ ${num_zerosegsitewindows} -gt 0 ]; then
		# We remove rows in the obs and scaffoldlist files based on the line numbers indicated in the removal file.
		awk 'NR==FNR{a[$1];next}!(FNR in a)' zerosegsites_rows_to_remove_wHeader.list ${obs_file} > ${prefix}.zeroSegSitesRemoved.obs
		awk 'NR==FNR{a[$1];next}!(FNR in a)' zerosegsites_rows_to_remove_noHeader.list ${scaflist_file} > ${prefix}.zeroSegSitesRemoved.scaffoldlist
		# Or to do directly without separate file (but only for the obs file)
		#awk '($1 != 0 ) ' ${obs_file} > ${prefix}.zeroSegSitesRemoved.obs
	elif [ ${num_zerosegsitewindows} -eq 0 ]; then
		cp ${obs_file} ${prefix}.zeroSegSitesRemoved.obs
		cp ${scaflist_file} ${prefix}.zeroSegSitesRemoved.scaffoldlist
	fi

	# Change back into parent directory
	cd ..

done

