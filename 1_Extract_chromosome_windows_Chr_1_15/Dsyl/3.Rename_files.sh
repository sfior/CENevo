#!/bin/bash

# Rename files

#for i in `seq 1 10`; do
for i in 2 3 4 5 6 7 8 9 10 11 12 14 15; do 
	mv ACAGTG.Chr${i}.mpileup HIGH_ACAGTG.Chr${i}.mpileup
	mv CAGATC.Chr${i}.mpileup HIGH_CAGATC.Chr${i}.mpileup
	mv CTTGTA.Chr${i}.mpileup HIGH_CTTGTA.Chr${i}.mpileup
	mv AGTCAA.Chr${i}.mpileup LOW_AGTCAA.Chr${i}.mpileup
	mv GCCAAT.Chr${i}.mpileup LOW_GCCAAT.Chr${i}.mpileup
	mv GTGAAA.Chr${i}.mpileup LOW_GTGAAA.Chr${i}.mpileup
done
