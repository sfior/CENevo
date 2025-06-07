# Take Dsyl seqs including 3 major haplotypes + 3 radom seqs


file1='0_sample_list_1165.txt'
file2='Exons_CEN.fasta'

for n in {1..500}
do
cat $file1 | head -n3 > SETS${n}.txt
cat $file1 | tail -n 162 | grep Ds | gshuf | head -n3 >> SETS${n}.txt

#Take other groups
groups='Outgroup Group_2 Group_3 Group_5 Group_6'
for g in $groups
do
	#echo $g
	cat $file1 | grep $g | gshuf | head -n6 >> SETS${n}.txt
done
done