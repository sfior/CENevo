POPS=$(cat Dsyl_barcodes.txt | cut -f2 -d' ')

for pop in $POPS
do
cat ${pop}_w500_s500.thetasWindow.pestPG |  awk '$14 > 250' > ${pop}'_w500_s500_filter_250bp.thetasWindow.pestPG'
done

