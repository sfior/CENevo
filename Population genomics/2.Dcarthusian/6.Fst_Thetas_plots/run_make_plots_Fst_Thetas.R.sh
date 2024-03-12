#module load r/3.2.2
# for testing

Rscript combined_significance_of_loess_v3_15pw_larger_plots.R \
list_selected_scaffolds.txt \
combined_pairs_filter_250bp_sorted .fst \
.001 .1 1000 \
"V1","V2","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18" \
"Pair1_2","Pair1_3","Pair1_4","Pair1_5","Pair1_6","Pair2_3","Pair2_4","Pair2_5","Pair2_6","Pair3_4","Pair3_5","Pair3_6","Pair4_5","Pair4_6","Pair5_6" \
"Pair1_2","Pair1_4","Pair1_5","Pair2_3","Pair2_6","Pair3_4","Pair3_5","Pair4_6","Pair5_6","Pair2_4","Pair2_5","Pair4_5","Pair1_3","Pair1_6","Pair3_6" \
.05 .95 \
pdfyes 500 \
_w500_s500_filter_250bp_sorted .tP list_selected_scaffolds.txt \
"Pair1_2","Pair1_4","Pair1_5","Pair2_3","Pair2_6","Pair3_4","Pair3_5","Pair4_6","Pair5_6","Pair2_4","Pair2_5","Pair4_5","Pair1_3","Pair1_6","Pair3_6" '' \
1000 .05 .95 .001 .1 'tPfiles/' \
_w500_s500_filter_250bp_sorted .Tajima list_selected_scaffolds.txt \
"Pair1_2","Pair1_4","Pair1_5","Pair2_3","Pair2_6","Pair3_4","Pair3_5","Pair4_6","Pair5_6","Pair2_4","Pair2_5","Pair4_5","Pair1_3","Pair1_6","Pair3_6" '' \
1000 .05 .95 .001 .1 'Tajimafiles/' \
'Dcar' \
intersect.txt 500 500 50 0.5 .001 .1 \
"pop1.pop2","pop1.pop4","pop1.pop5","pop2.pop3","pop2.pop6","pop3.pop4","pop3.pop5","pop4.pop6","pop5.pop6","pop2.pop4","pop2.pop5","pop4.pop5","pop1.pop3","pop1.pop6","pop3.pop6" \
'snps.weighted.mean' \
"HL","HL","HL","HL","HL","HL","HL","HL","HL","HH","HH","HH","LL","LL","LL" \
combined_scaffolds_evidence.pdf
