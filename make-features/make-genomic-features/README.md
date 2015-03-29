# Make addition genomic features

#### Dependancies
- bedtools

#### Convert SNP positions to bed format

The script `1_snp-to-bed.py` can be used to convert a GWAS summary stat file into a list of SNP position in bed format. Open the script and edit fields that specify columns to use.

The script reads from stdin and writes to stdout.

```
cat <in file> | python 1_snp-to-bed.py > <out prefix.bed>
```
