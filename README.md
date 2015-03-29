# GWAS Enrichment
Cross-trait enrichment using GWAS summary stats

#### Dependancies
- Python (>=3.2)
- pandas
- numpy
- HTSeq
- statsmodels

Require python packages can be install using pip

```
pip install numpy pandas HTSeq statsmodels
```

#### Usage

```
usage: gwas-enrichment.py [-h] [--features [<txt file> [<txt file> ...]]]
                          [--out <out prefix>] [--log]
                          [--test {lower,upper,two}] [--topthresh <float>]
                          [--nullsize <int>] [--testposrange <int kb>]
                          [--ldmap <plink ld file>] [--maxr2 <float>]
                          [--mafrange <float>] [--posrange <int kb>]
                          [--ldscbins <int>] [--testcol <str>]
                          [--snpcol <str>] [--bpcol <str>] [--chrcol <str>]
                          [--pcol <str>] [--mafcol <str>] [--sep <str>]
                          [--missing <str> [<str> ...]]
                          <reference sum stats> <test sum stats> <ld score
                          file> [<ld score file> ...]

GWAS stats enrichment

positional arguments:
  <reference sum stats>
                        GWAS summary stats for selecting top SNPS and creating
                        null.
  <test sum stats>      GWAS summary stats to test for enrichment.
  <ld score file>       List of LD score files to be merged. LD scores can be
                        generated using ldscr from
                        https://github.com/bulik/ldsc/#where-can-i-get-ld-
                        scores

optional arguments:
  -h, --help            show this help message and exit
  --features [<txt file> [<txt file> ...]]
                        Genomic features to match profile on. Each file should
                        be a list of SNP names which share a feature.
  --out <out prefix>    Output prefix. (default: out)
  --log                 Store log rather than print to stdout.
  --test {lower,upper,two}
                        Select one-tailed [upper/lower] or two-tailed test.
                        (default: lower)
  --topthresh <float>   SNPs below this p-value threshold will be used for
                        test. (default: 1e-6)
  --nullsize <int>      Number of SNPs to sample for null distribution.
                        (default: 100,000)
  --testposrange <int kb>
                        Test SNPs are excluded if there is another test SNP
                        with higher p-value within this genomic range. Useful
                        to reduce multiple testing burden. (default: None)
  --ldmap <plink ld file>
                        LD map output from plink --r2. Used to exclude test
                        SNPs in high LD. Useful to reduce multiple testing
                        burden.
  --maxr2 <float>       Test SNPs with r2 higher than this will be excluded.
                        Only used in conjunction with --ldmap. (default: 0.7)
  --mafrange <float>    Max MAF deviation from reference SNP for inclusion in
                        null. (default: 0.02)
  --posrange <int kb>   Null sample SNPs must be further than this distance
                        from reference SNP. (default: 1000 kb)
  --ldscbins <int>      Number of bins to make for LD scores. (default: 10)
  --testcol <str>       Column with test-statistic. (default:
                        frequentist_add_pvalue)
  --snpcol <str>        Column with SNP id. (default: rsid)
  --bpcol <str>         Column with SNP position. (default: position)
  --chrcol <str>        Column with chromosome number. (default:
                        alternate_ids)
  --pcol <str>          Column with p-value significance. (default:
                        frequentist_add_pvalue)
  --mafcol <str>        Column with MAF. (default: all_maf)
  --sep <str>           Column separater. (default: tab)
  --missing <str> [<str> ...]
                        List of values to use as missing. (default: . "" NA)
```