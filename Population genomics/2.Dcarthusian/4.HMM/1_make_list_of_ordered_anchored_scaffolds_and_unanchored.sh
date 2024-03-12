## This script makes list of scaffolds present in the Fst estimates (i.e., drops scaffold with no coverage)

echo '############ creating list of all scaffolds in each chromosome ###############'
for i in {1..15}
do
echo '**** chromosome ' ${i} ' ****'
cat Dcart_all_anchored_scaffolds_v2.txt | grep 'Dcart_'${i}'_ml' | sort -k3 -g | cut -f1 > 'chr_'${i}'.tmp'
done

echo '############ now creating list of scaffolds with windows in each chromosome ###############'
for i in {1..15}
do
echo '**** chromosome ' ${i} ' ****'
rm 'chr_'${i}
touch 'chr_'${i}
SCAFFS=$(cat 'chr_'${i}'.tmp')
for scaff in $SCAFFS
do
x=$(cat list_of_scaffolds_all.txt | grep ${scaff} | wc -l)
	if [ $x = 1 ]; then
	echo ${scaff} >> 'chr_'${i}
	else
	echo ${scaff} is missing in windows file - not writing it in final chromosme list for hmm
	fi
	done
done
rm chr*.tmp


## List of unanchored scaffolds in chr_16


echo '############ now creating list of unanchored scaffolds with windows  ###############'

cp list_of_scaffolds_all.txt tmp1
n=0
for i in {1..15}
do 
echo '**** chromosome ' ${i} ' ****'
SCAFFS=$(cat 'chr_'${i} | cut -f1)
	for scaff in ${SCAFFS}
	do
	n=$(( $n + 1 ))
	o=$(( $n + 1 ))
	#echo $n
	#echo ${scaff}
	cat tmp${n} | grep -v ${scaff} > tmp${o}	
done
done

mv tmp${o} chr_16
rm tmp*

