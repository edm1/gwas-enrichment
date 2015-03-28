#!/usr/bin/env bash
#
# Convert ensembl gene gtf to bed format

insumstats=../../data/visacutiy_best_acuity_ex4sd_unrel1_warped_f7_chr20-21.tsv.gz
outbed=alspac_maf0.01_info0.8_snps.bed

zcat $insumstats | python 1_snp-to-bed.py > $outbed
