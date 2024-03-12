# This scipt makes temporary folders and links bam files used downstream


DIR='xxx'  # specify path of directory containign bam files


NJOBS=199

for n in $(seq 1 $NJOBS)
do
T='./tmp.'${n}'/'
mkdir ${T}
for i in $(ls -1 ${DIR}*sorted.bam* |  awk -F '/' '{print $NF}' )
do
#echo ${i}
ln -s ${DIR}${i} ${T}${i}
done
done
