#!/usr/bin/env bash
#
# Convert ensembl gene gtf to bed format

ingtf=Homo_sapiens.GRCh37.75.gtf.gz
outbed=Homo_sapiens.GRCh37.75_exons.bed

zcat $ingtf | awk -v OFS="\t" '$3 == "exon" {print $1, $4, $5, ".", ".", $7}' > $outbed
