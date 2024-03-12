# Run steps in sequence from command line


for i in $(seq 1 16)
do
dir='hmm.chr'${i}
mkdir ${dir}
mv *'_'${i}'.hmm.txt' ${dir}
mv *'_'${i}'.Summarypar.txt' ${dir}
done


cwd=$(pwd)
for i in $(seq 1 16)
do
cd 'hmm.chr'${i}
ln -s ../fst_comps.txt
pairs=$(cat fst_comps.txt )
pairs_hmm=$(echo $pairs | sed "s/ /_chr_${i}.hmm.txt /g" | sed "s/pop6$/pop6_chr_${i}.hmm.txt/")
paste $pairs_hmm > combined_hmm_all_pairs.txt
cat combined_hmm_all_pairs.txt | cut -f1,2,3,4,6,10,12,16,18,22,24,28,30,34,36,40,42,46,48,52,54,58,60,64,66,70,72,76,78,82,84,88,90 > combined_hmm_all_pairs_reduced.txt 
cd $cwd
done

# command 4
module load gdc r/3.2.2
cwd=$(pwd)
for i in $(seq 1 16)
do
cd 'hmm.chr'${i}
Rscript  ../4_overlap_high_differentiation_states.R combined_hmm_all_pairs_reduced.txt "pop1.pop2","pop3.pop5","pop4.pop6","pop1.pop4","pop1.pop5","pop2.pop3","pop2.pop6","pop3.pop4","pop5.pop6" Overlapping_Outliers_ite_500
cd $cwd
done


cwd=$(pwd)
for i in $(seq 1 16)
do
cd 'hmm.chr'${i}
cat Overlapping_Outliers_ite_500_pop1.pop2_pop3.pop5_pop4.pop6_pop1.pop4_pop1.pop5_pop2.pop3_pop2.pop6_pop3.pop4_pop5.pop6.txt | cut -f1 | sort | uniq | grep scaff > list_scaffolds_with_hmm_outliers.txt
cd $cwd
done

# command 5
module load gdc python/2.7.6
cwd=$(pwd)
for i in $(seq 1 16)
do
cd 'hmm.chr'${i}
python ../5_Extract_overlapping_windows_v2.py list_scaffolds_with_hmm_outliers.txt Overlapping_Outliers_ite_500_pop1.pop2_pop3.pop5_pop4.pop6_pop1.pop4_pop1.pop5_pop2.pop3_pop2.pop6_pop3.pop4_pop5.pop6.txt 500
cd $cwd
done


# command 6
module load gdc python/2.7.6
cwd=$(pwd)
for i in $(seq 1 16)
do
cd 'hmm.chr'${i}
python ../6_define_continuous__sweeps_because_of_lack_of_coverage.py    list_scaffolds_with_hmm_outliers.txt   ../combined_pairs_filter_250bp_sorted.fst   500   candidate.hmm.continous.sweeps.txt
cd $cwd
done



outfile='candidate.hmm.continous.sweeps.per.chrm.txt'
rm ${outfile}
touch ${outfile}
for i in $(seq 1 16)
do
echo 'chromosome '${i} >> ${outfile}
cat 'hmm.chr'${i}'/candidate.hmm.continous.sweeps.txt' >> ${outfile}
cd $cwd
done


dir='hmm.combined.chrs'
mkdir ${dir}
out=${dir}'/hmm_combined_all_chr_reduced.txt'
rm ${out}
touch ${out}
cat hmm.chr1/combined_hmm_all_pairs_reduced.txt | head -n1 >> ${out}
for i in $(seq 1 16)
do
cat 'hmm.chr'${i}'/combined_hmm_all_pairs_reduced.txt' | grep -v 'chr' >> ${out}
done
cat candidate.hmm.continous.sweeps.per.chrm.txt | grep scaffold | cut -f1 > ${dir}/list_scaffolds_in_hmm_peaks.txt









