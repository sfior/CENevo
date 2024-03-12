bsub -n1 -W4:00 -R "rusage[mem=20000]" "Rscript make_1DSFS_poolstat.R poolstatoutput.1-20_sites_lowerQ20eq0.txt"
