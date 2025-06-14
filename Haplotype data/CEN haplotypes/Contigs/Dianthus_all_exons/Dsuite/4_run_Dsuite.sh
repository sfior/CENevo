# Run Dtrios for iterations

for n in {1..500}
do
Dsuite/Build/Dsuite Dtrios inputFile.Dsuite${n}.vcf SETS${n}.txt
done
