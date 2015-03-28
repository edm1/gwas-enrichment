#!/usr/bin/env bash
#
# Convert ensembl gene gtf to bed format and find overlapping SNPs
# Usage:
#   bash 2_gtf-to-bed.sh <gtf> <snp bed> <out pref>

# In and out names

ingtf=$1    # GTF file with genomic features
insnpbed=$2 # SNP positions in bed format
outpref=$3  # Output prefix

# ingtf=Homo_sapiens.GRCh37.75.gtf.gz      # GTF file with genomic features
# insnpbed=alspac_maf0.01_info0.8_snps.bed # SNP positions in bed format
# outpref=Homo_sapiens.GRCh37.75           # Output prefix

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract genomic features in bed format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Source type - to see possibilities use:
#   zcat <gtf.gz file> | awk '!/^#/ {print $2}' | sort | uniq
sourcetype=protein_coding

# Extract transcripts in bed format
featuretype=transcript
outtrans=$outpref"_"$sourcetype"_"$featuretype".bed"
zcat $ingtf | awk -v OFS="\t" -v stype=${sourcetype} -v ftype=${featuretype} '$2 == stype && $3 == ftype {print $1, $4, $5, ".", ".", $7}' > $outtrans

# Extract exons in bed format
featuretype=exon
outexon=$outpref"_"$sourcetype"_"$featuretype".bed"
zcat $ingtf | awk -v OFS="\t" -v stype=${sourcetype} -v ftype=${featuretype} '$2 == stype && $3 == ftype {print $1, $4, $5, ".", ".", $7}' > $outexon

# Subtract exons from transcripts to get introns in bed format
outintron=$outpref"_"$sourcetype"_introns.bed"
bedtools subtract -a $outtrans -b $outexon > $outintron

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find SNPs that overlap with genomic features
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find snps that are in exons
outname=snpinexon.txt
bedtools intersect -wa -a $insnpbed -b $outexon | awk '{print $4}' | sort | uniq > $outname

# Find snps that are in introns
outname=snpinintron.txt
bedtools intersect -wa -a $insnpbed -b $outintron | awk '{print $4}' | sort | uniq > $outname
