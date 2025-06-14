# Takes sequences and makes input files

file1='Exons_CEN.fasta'

for n in {1..500}
do
	> inputFile.Dsuite${n}.fasta
cut -f1 SETS${n}.txt | grep -F -f - "$file1" -A1 | grep -v '^--$' >> inputFile.Dsuite${n}.fasta
done