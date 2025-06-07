
for n in {1..500}
do
python3 fasta_to_vcf.py inputFile.Dsuite${n}.fasta inputFile.Dsuite${n}.vcf
done