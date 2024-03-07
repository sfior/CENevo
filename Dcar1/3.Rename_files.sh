#!/bin/bash

# Rename files
for i in `seq 1 15`; do
	mv CGTACG.Chr${i}.mpileup HIGH_CGTACG.Chr${i}.mpileup
	mv GTGGCC.Chr${i}.mpileup HIGH_GTGGCC.Chr${i}.mpileup
	mv GTTTCG.Chr${i}.mpileup HIGH_GTTTCG.Chr${i}.mpileup
	mv ATTCCT.Chr${i}.mpileup LOW_ATTCCT.Chr${i}.mpileup
	mv GATCAG.Chr${i}.mpileup LOW_GATCAG.Chr${i}.mpileup
	mv TTAGGC.Chr${i}.mpileup LOW_TTAGGC.Chr${i}.mpileup
done

