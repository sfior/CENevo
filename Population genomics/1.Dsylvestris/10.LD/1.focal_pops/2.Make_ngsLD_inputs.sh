#!/bin/bash

# This script generates the input (genotype and position) files for ngsLD from ANGSD beagle GLF outputs


#### CEN
# 95gwas data
prefix="GL_95gwas_CEN_region_w100kbFlanks_pVal12maf20depth2"

# Make genotype file
# As input, ngsLD accepts both genotypes, genotype likelihoods (GP) or genotype posterior probabilities (GP). Genotypes must be input as gziped TSV with one row per site and one column per individual n_sites.n_ind and genotypes coded as [-1, 0, 1, 2]. 
# As for GL and GP, ngsLD accepts both gzipd TSV and binary formats, but with 3 columns per individual 3.n_sites.n_ind and, in the case of binary, the GL/GP coded as doubles.
zcat ${prefix}.beagle.gz | cut -f 4- | tail -n+2 | gzip -c > ${prefix}.GLF.gz

# Make position file
# FILE: input file with site coordinates (one per line), where the 1st column stands for the chromosome/contig and the 2nd for the position (bp); remaining columns will be ignored but included in output; --posH assumes there is a header.
zcat ${prefix}.beagle.gz | cut -f 1 | awk -F "_" -v OFS="_" '{ print $1, $2}' | tail -n+2 > scaffolds.temp #CEN
zcat ${prefix}.beagle.gz | cut -f 1 | awk -F "_" -v OFS="\t" '{ print $3}' | tail -n+2 > positions.temp #CEN
paste -d"\t" scaffolds.temp positions.temp > ${prefix}.pos
rm scaffolds.temp positions.temp


