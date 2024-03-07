##### This script generates a genome index file from the fasta index file (of the reference assembly) ####

### This is useful for fast sorting and runtimes. For instructions and reference of how to generate a genome index file from a fasta index file, see: https://www.biostars.org/p/70795/)

ref_directory="/cluster/project/gdc2/special/shared/p461/Assemblies/Dsylvestris"

# 1) Use samtools to generate fasta index (.fai file)
#samtools faidx ${ref_directory}/assembly_homozygous.fa
# 2) This index file won't work as genome file due to file format issue (mainly more than required number of columns), so we modify the index file
awk -v OFS='\t' {'print $1,$2'} ${ref_directory}/assembly_homozygous.fa.fai > ./Dsyl_assembly_homozygous.genome
# This prints the 1st and 2nd column of a fai index file and separates the column by tab (OFS flag). Use this file as genome file in bedtools.
# If space desired between columns do this
#awk {'print $1,"",$2'} ${ref_directory}/assembly_homozygous.fa.fai > ./Dsyl_assembly_homozygous.genome
# if 'chr' needs to be added infront of the chromosome/scaffold names do this
#awk {'print "chr"$1,"",$2'} ${ref_directory}/assembly_homozygous.fa.fai > ./Dsyl_assembly_homozygous.genome
