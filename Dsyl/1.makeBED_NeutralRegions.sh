######## THIS SCRIPT OUTPUTS NEUTRAL REGIONS FROM A GFF ANNOTATION FILE ########

# 1. Load modules
module load python/3.6.1 bedtools/2.27.1 samtools/1.8
export PATH=$PATH:/cluster/project/gdc/people/lhirzi/BEDOPS/bin

# 2. Make and define output directory
mkdir 5000_1000
output="./5000_1000"

### 3. Remove items from gff annotation file if necessary. In the gff file all rows may be considered as annotated regions, apart from in this case rows containing “contig” (in the 3rd column) as these were created by Stefan’s pipeline as it worked on contigs instead of scaffolds. These actually report any possible annotated contig, so if you remove those possibly nothing is left. Hence let's remove this for now
cat Ds.final.28.10.14.contigs.v7.maker.contig.to.scaffold.gff3 | sed '/contig/d' > Ds_trimmed.gff
# It seems that the above gff file is still too liberal in defining annotated regions (i.e. once we do the steps below, nothing is left behind). As an alternative, lets grep only lines containing 'maker'. Maker is a second-step validation program which validates the first 'first' discovery of annotated regions (in second column).
#grep maker Ds.final.28.10.14.contigs.v7.maker.contig.to.scaffold.gff3 > Ds_trimmed.gff
# In case you want something more liberal than maker (numerous false negatives) but more conservative than all matches (few false positives), try an intermediate combination. Try playing around with the filtering, e.g. knowing that Augustus is very conservative in calling an annotated region while SNAP and Genemarker are very trigger happy. Here, let's try e.g. with maker annotated regions w/ 10kb flanks + all matches w/0kb flanks.

### 4. Convert GFF to bed format. We'll need a bed format to filter/extract from the bam file.
gff2bed < Ds_trimmed.gff > Ds_trimmed.bed

### 5. Maker validated regions have a high certainty of being correct (few false positives). Other annotated regions (non-maker approved) have a low certainty of being correct. Given this difference in certainty/accuracy of annotation, it makes sense to deal with these two groups of regions differently, e.g. apply long flanking regions for maker-validated regions and short or no flanking regions for non maker-validated regions. To do this, we thus split the bed file into maker and non-maker (others) validated bed files.
cat Ds_trimmed.bed | sed '/maker/d' > Ds_trimmed_others.bed
grep maker Ds_trimmed.bed > Ds_trimmed_maker.bed

### 6. Merge (union and flatten) sites. Many entries are overlapping or even duplicate, so we want to simplify the bed file by merging (taking the union) of entries
bedops --merge Ds_trimmed_others.bed  > Ds_trimmed_others_merged_temp.bed
bedops --merge Ds_trimmed_maker.bed  > Ds_trimmed_maker_merged_temp.bed

### 7. There appears to be inconsistency between the gff position indexing and bed. The gff2bed script converts 1-based, closed [start, end] General Feature Format v3 (GFF3) to sorted, 0-based, half-open [start-1, end) extended BED-formatted data.
# The problem here is that, for some reason, the gff file were working with is of the format closed [start+1, end+1], i.e. chromStart = 2 and chromEnd = 101 (instead of the expected chromStart = 1 and chromEnd = 100; recall gff2bed would then 'correct' this to chromStart = 1 and chromEnd = 101 instead of chromStart = 0 and chromEnd = 100, which is the correct BED format). Hence we need to subtract 1 from both chromStart and Chromend in the resulting bed file. 
python fix_gff2bed.py Ds_trimmed_others_merged_temp.bed -o Ds_trimmed_others_merged.bed
python fix_gff2bed.py Ds_trimmed_maker_merged_temp.bed -o Ds_trimmed_maker_merged.bed

### 8. Add flanking regions if desired. If we're interested in region far away from annotated regions/genes, we need to apply a flanking region to the annotated regions so that we can remove them later with the annotated regions. Here we run a python script whereby the the length of the flanking region (in bp) is defined by the -f tag.
# 5kb is a suitable length for the flanking regions (for maker-validated regions). 1kb is a suitable length for others.
python addFlankingRegions2BED.py Ds_trimmed_others_merged.bed -f 1000 -o ${output}/Ds_trimmed_others_merged_flanked.bed
python addFlankingRegions2BED.py Ds_trimmed_maker_merged.bed -f 5000 -o ${output}/Ds_trimmed_maker_merged_flanked.bed

### 9. Concatenate the two bed files
cat ${output}/Ds_trimmed_maker_merged_flanked.bed ${output}/Ds_trimmed_others_merged_flanked.bed > ${output}/Ds_trimmed_flanked.bed
#cat ${output}/Ds_trimmed_maker_merged_flanked.bed Ds_trimmed_others_merged.bed > ${output}/Ds_trimmed_flanked.bed

### 10. Sort bed file (https://www.biostars.org/p/64687/). Remember that sorting needs to precede merging.
cat ${output}/Ds_trimmed_flanked.bed | sort -k1,1V -k2,2n -k3,3n > ${output}/Ds_trimmed_flanked_sorted.bed

### 11. Merge again, as there is new potential overlap after the addition of flanking regions.
bedops --merge ${output}/Ds_trimmed_flanked_sorted.bed  > ${output}/Ds_final_temp.bed
cat ${output}/Ds_final_temp.bed | sort -k1,1V -k2,2n -k3,3n > ${output}/Ds_final.bed

### 12. Use bedtool complement to find the complement of the bed file. In case you want to use bedtools intersect (without -v) or samtools view -L; options which consider overlap of regions listed in a bed.
bedtools complement -i ${output}/Ds_final.bed -g Dsyl_assembly_homozygous.genome > ${output}/Ds_complement.bed


##### The following steps (13-15) are are performed by running the following script, subsequent to running this script ###
# bsub < FilterBAMnMpileup.lsf

#### 13a. Remove regions from a bam file using a bed file. Use bedtools intersect with the -v tag (only report those entries in A that have no overlap in B). ## -v -g with Ds_final and -f 1 with Ds_complement are equivalent.
## The BAM files are located at: /cluster/project/gdc/shared/p219/poolstat_AF_4_Hirzi/Dsylvestris/bamfiles
# Run 
#bedtools intersect -a /cluster/project/gdc/shared/p219/poolstat_AF_4_Hirzi/Dsylvestris/bamfiles/ACAGTG_realn_recal_sorted.bam -b ${output}/Ds_complement.bed -f 1.0 -sorted > ${output}/ACAGTG_regions_removed_bedtools.bam
#bedtools intersect -a /cluster/project/gdc/shared/p219/poolstat_AF_4_Hirzi/Dsylvestris/bamfiles/ACAGTG_realn_recal_sorted.bam -b ${output}/Ds_final.bed -v -sorted -g Dsyl_assembly_homozygous.genome > ${output}/ACAGTG_regions_removed_bedtoolsVg.bam\

### 13b. Note that another alternative, samtools view with the -L argument, which is equivalent to bedtools intersect with Ds_complement, doesn't give what you want, because it doesn't allow you to set the fraction of overlap as is possible with the bedtools -f and -F tags. i.e. it gives all partially overlapping alignments/reads.
#samtools view -L Ds_complement.bed /cluster/project/gdc/shared/p219/poolstat_AF_4_Hirzi/Dsylvestris/bamfiles/ACAGTG_realn_recal_sorted.bam -o ./ACAGTG_regions_removed.bam

### 13c. For fast sorting and runtimes, we need the genome index, which can be produced from the fasta index file (see https://www.biostars.org/p/70795/). E.g. to make a genome file (for bedtools) using reference genome
# 1) Use samtools to generate fasta index (.fai file)
#samtools faidx test_genome.fa
# 2) This index file won't work as genome file due to file format issue (mainly more than required number of columns), so we modify the index file
#awk -v OFS='\t' {'print $1,$2'} /cluster/project/gdc2/special/shared/p461/Assemblies/Dsylvestris/assembly_homozygous.fa.fai > ./Dsyl_assembly_homozygous.genome
# This prints the 1st and 2nd column of a fai index file and separates the column by tab (OFS flag). Use this file as genome file in bedtools.
# If space desired between columns do this
#awk {'print $1,"",$2'} test_genome.fa.fai > test.genome
# if 'chr' needs to be added infront of the chromosome/scaffold names do this
#awk {'print "chr"$1,"",$2'} test_genome.fa.fai > test.genome

### 14. To produce mpileup file from BAM
#samtools mpileup ACAGTG_regions_removed_bedtoolsVg.bam > ACAGTG_regions_removed_bedtoolsVg.mpileup

### 15. To filter for anchored scaffolds (21369 scaffolds total; 1359 are anchored). To get a unique list of all scaffolds: cut -f1 -d$'\t' Dsyl_assembly_homozygous.genome | sort | uniq 
#for i in $(cut -f1 -d$'\t' Dsyl_all_scaffolds_v2.txt | tail -n+2); do grep ${i} ${output}/ACAGTG_regions_removed_bedtoolsVg.mpileup >> ${output}/ACAGTG_regions_removed_bedtoolsVg_anchored.mpileup; done



###################### To test this various aspects of this script, you can work on a simulated BAM file (https://bitbucket.org/phaentu/atlas/wiki/Auxiliary%20Tools:%20simulate). writeTrueGenotypes outputs the the true genotypes and bed files ######################
#atlas task=simulate out=example depth=3 chrLength=1000 writeTrueGenotypes
## General commands
## To convert BAM to SAM and view
#samtools view -H
## To view header of BAM
#samtools view -h

### Debugging.
## If you have: "Error: Sorted input specified, but the file Ds_final.bed has the following record with a different sort order than the genomeFile Dsyl_assembly_homozygous.genome
#scaffold1000_size94247	0	94247", you need to make sure the sort order of the bed file and genome file are the same.