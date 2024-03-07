#!/bin/bash

### Define used paths
results_raw="/cluster/scratch/lhirzi/ABCtoolbox_results"
results_summary="/cluster/scratch/lhirzi/ABCtoolbox_results/Simulation_results_summary"
folder_prefix="run_5000_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_"

### Concatenate results from mutiple, parallel simulations into 1 file
head -n1 ${results_raw}/${folder_prefix}1/*sampling1.txt > ${results_summary}/concatenated_results_temp.txt 
for folder in $(ls -d ${results_raw}/${folder_prefix}*); 
do tail -n+2 ${folder}/*sampling1.txt >> ${results_summary}/concatenated_results_temp.txt;
done

### Modify concatenated file. Revise numbering in first column so each row index is unique.
wordcount_raw=$(cat ${results_summary}/concatenated_results_temp.txt | wc -l)
wordcount=$((${wordcount_raw} - 1))
echo Sim > ${results_summary}/first_column
seq 1 ${wordcount} >> ${results_summary}/first_column
# Print all columns except the first one. Another e.g.: awk '{$6=$8=""; print $0}' file to print all but the 6th and 8th columns. We also set the output field separator (OFS) to a tab.
awk '{$1=""; print $0}' ${results_summary}/concatenated_results_temp.txt > ${results_summary}/concatenated_results_temp2.txt
paste -d' ' ${results_summary}/first_column ${results_summary}/concatenated_results_temp2.txt > ${results_summary}/concatenated_results.txt

### Remove temporary files
rm ${results_summary}/first_column ${results_summary}/concatenated_results_temp2.txt ${results_summary}/concatenated_results_temp.txt
