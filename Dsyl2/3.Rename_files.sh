#!/bin/bash

# Make folders. Move files into folders. Rename files
for i in `seq 11 15`; do
	mkdir split_${i}
	mv *.Chr16.${i}.* split_${i}
	mv split_${i}/ACAGTG.Chr16.${i}.mpileup split_${i}/HIGH_ACAGTG.Chr16.${i}.mpileup
	mv split_${i}/CAGATC.Chr16.${i}.mpileup split_${i}/HIGH_CAGATC.Chr16.${i}.mpileup
	mv split_${i}/CTTGTA.Chr16.${i}.mpileup split_${i}/HIGH_CTTGTA.Chr16.${i}.mpileup
	mv split_${i}/AGTCAA.Chr16.${i}.mpileup split_${i}/LOW_AGTCAA.Chr16.${i}.mpileup
	mv split_${i}/GCCAAT.Chr16.${i}.mpileup split_${i}/LOW_GCCAAT.Chr16.${i}.mpileup
	mv split_${i}/GTGAAA.Chr16.${i}.mpileup split_${i}/LOW_GTGAAA.Chr16.${i}.mpileup	
done
