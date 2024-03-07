#!/bin/bash

working_dir="GenomeScan_resultsTables"

### Concatenate results from mutiple, parallel simulations into 1 file
head -n1 ${working_dir}/genomeScan_table.Chr1.txt > ${working_dir}/genomeScan_table.allChr1_15.txt
for file in $(ls ${working_dir}/genomeScan_table.Chr* | sort -V | head -15); 
do tail -n+2 ${file} >> ${working_dir}/genomeScan_table.allChr1_15.txt;
done
