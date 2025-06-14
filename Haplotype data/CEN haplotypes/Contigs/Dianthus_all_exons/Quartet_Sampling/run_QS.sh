# run QS
python3 quartet_sampling.py --tree prunedTree.nwk --align inputFile.QS.phy --reps 200 --threads 2 \
--engine raxml --engine-exec raxmlHPC-AVX --result-prefix RESULT_200reps

# Plot tree
Rscript plot_QC_ggtree.R -c RESULT_200reps.labeled.tre.qc -d RESULT_200reps.labeled.tre.qd -i RESULT_200reps.labeled.tre.qi -o ./
