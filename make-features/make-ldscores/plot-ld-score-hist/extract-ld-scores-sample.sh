ldfiles=/home/em13383/gwas-enrichment/data/ldscores/*.l2.ldscore.gz
out=sample-ld-scores.txt.gz
n=100000

zcat $ldfiles | awk 'NR != 1 {print $4}' | shuf | head -n $n | gzip -c > $out
