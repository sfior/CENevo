#!/bin/bash

### It may be the case (albeit very rarely) that some windows are returned with NA in the significance columns. This will lead to failure when plotting.
### A solution to this is to simply remove these NA windows, using this script.

# Define files and working directory
working_dir="GenomeScan_resultsTables"
file="genomeScan_table.Chr16.txt"
stem=$(basename "${file}" .txt)

# Identify lines with NAs in the significance columns. NAs will result in error during plotting.
cat ${working_dir}/${file} | awk '{ print $6 }' | grep -n "NA" >> ${working_dir}/Chr16_lines_to_remove.Nans.temp
cat ${working_dir}/${file} | awk '{ print $7 }' | grep -n "NA" >> ${working_dir}/Chr16_lines_to_remove.Nans.temp
cat ${working_dir}/${file} | awk '{ print $8 }' | grep -n "NA" >> ${working_dir}/Chr16_lines_to_remove.Nans.temp

# Reformat and make unique list
cat ${working_dir}/Chr16_lines_to_remove.Nans.temp | awk -F':' '{ print $1 }' > ${working_dir}/Chr16_lines_to_remove.Nans.temp2
cat ${working_dir}/Chr16_lines_to_remove.Nans.temp2 | sort -V | uniq > ${working_dir}/Chr16_lines_to_remove.Nans

# Print number of lines removed
num_lines_removed=$(cat ${working_dir}/Chr16_lines_to_remove.Nans | wc -l)
echo ${num_lines_removed} "window(s) were removed because of NAs"

# Make a copy of original file
cp ${working_dir}/${file} ${working_dir}/${stem}.NAsRemoved.txt

# Remove lines with NAs. Remember to sort (remove) in reverse order of index!
for line in $(cat ${working_dir}/Chr16_lines_to_remove.Nans | sort -rV); do
sed -e "${line}d" ${working_dir}/${stem}.NAsRemoved.txt > ${working_dir}/${stem}.NAsRemoved.temp.txt
mv ${working_dir}/${stem}.NAsRemoved.temp.txt ${working_dir}/${stem}.NAsRemoved.txt
done

# Remove temporary intermediate files
rm ${working_dir}/Chr16_lines_to_remove.Nans.temp ${working_dir}/Chr16_lines_to_remove.Nans.temp2


