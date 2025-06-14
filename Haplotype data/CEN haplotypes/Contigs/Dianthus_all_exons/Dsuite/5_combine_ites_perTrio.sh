# This script collects trios from the iterations 

n=500   #set number of iterations

s='BBAA'

outfile="collected_trios_${s}.txt"
> $outfile  

# Writes header
cat SETS1_${s}.txt | head -n1 >> $outfile

trios=$(cat "unique_trio_names_${s}.txt")
echo ">>>>>> Unique trios:" $trios
for t in $trios
do
	echo ">>>>>>>>> Taking this trio:" $t "<<<<<<<<<<<<"
	P1=$(echo $t | cut -f1 -d',')
	P2=$(echo $t | cut -f2 -d',')
	P3=$(echo $t | cut -f3 -d',')
	#echo $P1 
	#echo $P2
	#echo $P3
	eval cat SETS{1.."$n"}_${s}.txt | grep "$P1" | grep "$P2" | grep "$P3" >> $outfile
	echo >> $outfile
	echo
done
