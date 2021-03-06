#!/bin/sh 

cd $PBS_O_WORKDIR

script=/panfs/panasas01/sscm/em13383/programs/ldsc/ldsc.py
geno_pattern=/panfs/panasas01/sscm/em13383/phd/alspac/data/ldscores/make-ldscores/temp_geno/data_chrCHROM

# Make merge list
for chrom in $(seq -w 1 23); do
	inbfile=${geno_pattern/CHROM/$chrom}
	cmd="python $script --bfile $inbfile --l2 --ld-wind-cm 1 --out $chrom"
	echo $cmd
done
