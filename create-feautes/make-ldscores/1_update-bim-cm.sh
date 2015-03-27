#!/bin/sh 
#PBS -l nodes=1:ppn=16
#PBS -l walltime=00:00:10:00
#PBS -N n-update-cm
#PBS -q testq

cd $PBS_O_WORKDIR

geno_pattern=/panfs/panasas01/sscm/em13383/alspac_data/27Feb2015/data/genotypes/bestguess/subset_maf0.01_info0.8/data_chrCHROM
map_pattern=/panfs/panasas01/sscm/em13383/phd/alspac/data/ldscores/make-ldscores/1000GP_genetic-maps/genetic_map_chr@_combined_b37.txt
exclude_samples=/panfs/panasas01/sscm/em13383/phd/alspac/data/filter-files/mothers-and-custom-related_150323.txt
outdir=temp_geno

mkdir $outdir

# Make merge list
gmnum=0
for chrom in $(seq -w 1 22); do
	gmnum=$[$gmnum +1]
	inbfile=${geno_pattern/CHROM/$chrom}
	#ingm=${map_pattern/CHROM2/$gmnum}
	outbfile=$outdir/$(basename $inbfile)
	cmd="plink --bfile $inbfile --cm-map $map_pattern --remove $exclude_samples --make-bed --out $outbfile"
	echo $cmd
done | time parallel -j 1
