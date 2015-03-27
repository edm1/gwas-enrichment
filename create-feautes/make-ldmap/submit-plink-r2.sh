#!/bin/sh 
#PBS -l nodes=1:ppn=16
#PBS -l walltime=00:06:00:00
#PBS -N n-plink-r2
#PBS -m abe
#PBS -M edward.mountjoy@bristol.ac.uk
##PBS -q testq

cd $PBS_O_WORKDIR

geno_pattern=/panfs/panasas01/sscm/em13383/alspac_data/27Feb2015/data/genotypes/bestguess/data_chrCHROM
exclude_samples=/panfs/panasas01/sscm/em13383/phd/alspac/data/filter-files/mothers-and-custom-related_150323.txt

# Make merge list
ml=merge_list.txt
for chrom in $(seq -w 1 23); do
	inbfile=${geno_pattern/CHROM/$chrom}
	echo $inbfile
done > $ml

# Merge and extract
outbfile=unrelated-children_maf0.01_info0.8
time plink --merge-list "$ml" --out "$outbfile" --remove "$exclude_samples" --r2 --ld-window-kb 1000 --ld-window-r2 0.5
