#This scripts extracts gene information for candidate regions from genome annotation files

#CANDIDATE_REGIONS.bed is a list of genomic regions with format:
#chromosome	start	end

#Annotation files are available at Dryad repository https://doi.org/10.5061/dryad.x0k6djhng:
#FILE_annot.gff3
#FILE_maker_proteins.fa
#FILE_abinitio_proteins.fa
#FILE_maker_blast2go.txt
#FILE_abinitio_blast2go.txt
# make intersection with bedtools
   bedtools intersect -loj -a CANDIDATE_REGIONS.bed   -b FILE_annot.gff3  > intersection.txt 

# then e.g. do:
   grep -P "\tmaker\tgene" intersection.txt > intersection.maker.gene.txt
   grep -P "\tgene\t" intersection.txt > intersection.gene.txt 
   grep -P "\tmaker\t" intersection.txt > intersection.all.maker.predicted.txt
   grep -v "Parent=" intersection.txt |grep -v repeatmasker |cut -f1,2,3,7,8 > intersection.norepeat.short.txt
   grep -v "Parent=" intersection.txt |grep -v repeatmasker > intersection.norepeat.long.txt
paste intersection.norepeat.short.txt intersection.norepeat.long.txt > intersection.norepeat.comb.txt
   sort -k1,5 intersection.norepeat.comb.txt > intersection.norepeat.comb.sorted.txt


# create list of abinitio non-overlapping genes
   grep  -P "\tmatch\t" intersection.norepeat.comb.sorted.txt > intersection.abinit.non-overlap.gene.txt

# add gene name 
paste <(rev intersection.maker.gene.txt |cut -f1 -d"=" |rev) intersection.maker.gene.txt  > intersection.maker.gene.name.txt
paste <(rev intersection.abinit.non-overlap.gene.txt |cut -f1 -d"=" |rev ) intersection.abinit.non-overlap.gene.txt > intersection.abinit.non-overlap.gene.name.txt

# link in the protein data from maker
   ln -s ../../../maker-results/dc/FILE_maker_proteins.fa
   ln -s ../../../maker-results/dc/FILE_abinitio_proteins.fa
# convert "abinit-gene" into "processed-gene" in gene names
   perl -p -i -e 's/abinit-gene/processed-gene/g' intersection.abinit.non-overlap.gene.name.txt    

# final list of maker approved genes
   cp intersection.maker.gene.name.txt final.maker.gene.name.txt

# create the final list of abinitio, non-overlapping genes. Not all genes in above intersection.abinit.non.overlap.gene.name.txt are actually non-overlapping !!!
   for g in `cut -f1 intersection.abinit.non-overlap.gene.name.txt `; do echo $g ; grep "$g" FILE_abinitio_proteins.fa ; done |grep ">"  |cut -f1 -d" " |tr -d ">" > final.non-overlapping.txt
 for f in `cat final.non-overlapping.txt `; do grep $f intersection.abinit.non-overlap.gene.name.txt ; done > final.abinit.non-overlap.gene.name.txt


# pull out proteins. Use extractFastaSeqs.pl from here: https://github.com/douglasgscofield/bioinfo
# maker approved:
   for t in `cut -f1 intersection.maker.gene.name.txt` ; do grep $t FILE_maker_proteins.fa ; done |cut -f1 -d " " > tmp 
	cat tmp | sed 's/>//g' > tmp2
	/My_programs/douglasgscofield/bioinfo/scripts/extractFastaSeqs.pl -i FILE_maker_proteins.fa -n tmp2 -o OUTFILE.maker.proteins.fasta

# non-overlapping:
   for t in `cut -f1 final.abinit.non-overlap.gene.name.txt `; do grep $t FILE_abinitio_proteins.fa  ; done |cut -f1 -d" " > tmp
	cat tmp | sed 's/>//g' > tmp2
	/My_programs/douglasgscofield/bioinfo/scripts/extractFastaSeqs.pl -i FILE_abinitio_proteins.fa  -n tmp2 -o OUTFILE.maker.non_overlapping_ab_initio.proteins.fasta 
# pull out the blast2go annotation
	for g in `cut -f1 intersection.maker.gene.name.txt ` ; do echo -ne "$g\t" ; grep $g   FILE_maker_blast2go.txt  ; echo ; done |grep .  > OUTFILE.maker.gene.blast2go.annot.txt 

    for g in `cut -f1 final.abinit.non-overlap.gene.name.txt ` ; do echo -ne "$g\t" ; grep $g   FILE_abinitio_blast2go.txt  ; echo ; done |grep . > OUTFILE.abinitio.non-overlap.gene.blast2go.annot.txt 
