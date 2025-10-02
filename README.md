# CENevo

<img src="https://github.com/sfior/CENevo/blob/main/Figures/D.sylvestris_small.jpg" width="800">

This repository contains code and datasets associated with the paper "Ancient alleles drive contemporary climate adaptation in an alpine plant" ([Simone, Luqman et al. 2025](https://www.science.org/doi/10.1126/science.adp5717)). They are sorted into the following directories:

## Demography

ABC coalmod speciation data and scripts: Approximate Bayesian Computation (ABC) to estimate the speciation time for _D. sylvestris_ and the time to coalescence at the DsCEN/2 locus from a geographic sample of haplotypes.

Onset of selection: ABC workflow to estimate the time of onset of selection for candidate regions.

## Distribution modelling

Scripts to perform allele distribution models, with past and future temporal projections.

## Ecology

Data and scripts assessing phenotypic divergence in common gardens, and for evaluating DsCEN/2 allelic effects in transplant experiments.

## Genome assembly tools

Scaffold and gene names used in the study, and workflow to extract annotation for candidate loci. The genome assemblies, related annotation, and linkage maps are available [here](https://datadryad.org/dataset/doi:10.5061/dryad.x0k6djhng).

## Genome scans

LSD scans for detecting loci under selection and characterising genetic trade-offs, as parameterised under an empirically-informed continent-island model. The paper introducing the method (Luqman et al. 2021) can be found [here](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13415) and the GitHub repository is available [here](https://github.com/hirzi/LSD). 

## Genomic offset

Gene-environment association models based on gradient forests for predicting shifts in genome-wide adaptive composition (i.e. the genomic offset and glacial genomic offset+). The latter metric was derived in previous work ([Luqman et al. 2023](https://www.nature.com/articles/s41467-023-36631-9)). 

## Haplotype data

Data and scripts for analyses of CEN/2 and TCF1 haplotypes.

## Molecular biology

Alignment of FT gene family. Transcripts, data and scripts used in the _Arabidopsis_ transformation experiment.

## Population genomics

Population genomic analyses including [Poolstat](https://bitbucket.org/wegmannlab/poolstat/wiki/Home), computation of SFS, F~ST~ genome scans using HMM, computation of thetas, F~ST~ and θ plots for candidate regions, candidate annotation tables, GWAS on flowering time, LD for focal populations and across the species' range, environmental associations and IBD.

## Radiation introgression

ABBA-BABA data and scripts: tests for introgression across the _Dianthus_ radiation with D-statistics.

