#!/bin/bash

# Make filelist for each chromosome. Note that the order is important! The order is actually alphabetical so we can simply sort by name.

outdir="/cluster/project/gdc/people/lhirzi/GenomeScan/Dcar/Chr_16"

for i in `seq 1 20`; do
	cd ${outdir}/split_${i}
	ls *Chr16.${i}.mpileup | sort > Chr16.${i}_filelist.txt
done


