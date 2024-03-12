# This script set sites with likelihood of beeing polymorphix < 20 as invariant

zcat poolstatoutput.1-20_sites.txt.gz | awk 'BEGIN {FS=OFS="\t"}
{if($5< 20)
{print $1,$2,$3,$4,$5,"0","0","0","0","0","0"}
else {print $0}
}' > poolstatoutput.1-20_sites_lowerQ20eq0.txt

