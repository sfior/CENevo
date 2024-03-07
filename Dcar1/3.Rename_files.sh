#!/bin/bash

# Make folders. Move files into folders. Rename files
for i in `seq 1 20`; do
	mkdir split_${i}
	mv *.Chr16.${i}.* split_${i}
	mv split_${i}/CGTACG.Chr16.${i}.mpileup split_${i}/HIGH_CGTACG.Chr16.${i}.mpileup
	mv split_${i}/GTGGCC.Chr16.${i}.mpileup split_${i}/HIGH_GTGGCC.Chr16.${i}.mpileup
	mv split_${i}/GTTTCG.Chr16.${i}.mpileup split_${i}/HIGH_GTTTCG.Chr16.${i}.mpileup
	mv split_${i}/ATTCCT.Chr16.${i}.mpileup split_${i}/LOW_ATTCCT.Chr16.${i}.mpileup
	mv split_${i}/GATCAG.Chr16.${i}.mpileup split_${i}/LOW_GATCAG.Chr16.${i}.mpileup
	mv split_${i}/TTAGGC.Chr16.${i}.mpileup split_${i}/LOW_TTAGGC.Chr16.${i}.mpileup
done
