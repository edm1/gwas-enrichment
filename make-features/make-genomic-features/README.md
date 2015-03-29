## Make additional genomic features

#### Dependancies
- bedtools

#### Convert SNP positions to bed format

The script `1_snp-to-bed.py` can be used to convert a GWAS summary stat file into a list of SNP position in bed format. Open the script and edit fields that specify columns to use.

The script reads from stdin and writes to stdout.

```
cat <in file> | python 1_snp-to-bed.py > <outprefix_snps.bed>
```

#### Extract features from GTF file and find overlapping SNPs

An annotation of gene positions can be downloaded in GTF format from the [Ensembl FTP](http://www.ensembl.org/info/data/ftp/). If using 1000 genomes SNPs then use GRCh37 release 75.

```
bash 2_gtf-to-bed.sh <gtf gz> <snps.bed> <outprefix>
```

#### Output

- snpinexon.txt - list of SNP ids that are located in exons
- snpinintron.txt - list of SNP ids that are located in introns

You can then use `--features snpinexon.txt snpinintron.txt` to include the features in the enrichment analysis. The `2_gtf-to-bed.sh` could be adapted to extract different feature sets.
