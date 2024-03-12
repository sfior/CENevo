# This script takes Fst values from the output format of ANGSD for glabal Fst and makes a simmetric matrix.
# Input files must contain ' *_global.fst'

rm global.fst.matrix.tmp
touch global.fst.matrix.tmp

echo -e '\t'pop1'\t'pop2'\t'pop3'\t'pop4'\t'pop5'\t'pop6 >> global.fst.matrix.tmp

for m in {1..6}
	do
	rm tmp
	touch tmp
	echo 'pop'$m >> tmp
	
	for n in {1..6}
		do
		echo $m $n
		if [ $m = $n ]; then
			echo '0' >> tmp
		else
			fstfile=$(ls -1 *_global.fst | grep $m | grep $n)
			echo $fstfile
			fstvalue=$(cat $fstfile | cut -f2 -d' ')
			echo $fstvalue >> tmp
			
		fi
		done
	mat_line=$(cat tmp)
	echo $mat_line >> global.fst.matrix.tmp
done			


cat global.fst.matrix.tmp | tr ' ' '\t' > global.fst.matrix.txt
rm global.fst.matrix.tmp
rm tmp

cat global.fst.matrix.txt