#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import sys
import argparse
import pandas as pd
import numpy as np
from HTSeq import GenomicInterval
import gzip
from operator import itemgetter

def main():

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load first phenotype and sample null SNPs with similar profile to test
    # SNPs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Turn off annoying pandas warning
    pd.options.mode.chained_assignment = None 

    # Parse args
    global args
    args = parse_arguments()

    # Set up log
    if args.log == True:
        outname = outname = args.out + "_log.txt"
        log_h = open(outname, "w")
        sys.stdout = log_h

    # Load reference df
    print("Loading reference sum stats...")
    refdf = load_sumstats(args.refstats)
    print(" done")

    # Merge LD scores to reference df
    print("Loading LD scores...")
    load_merge_ldscores(refdf, args.ldscstats)
    print(" done")

    # Load additional genomic features
    # TODO

    # Split refdf into test snps vs null snps
    print("Splitting test and null snps...")
    istestsnp = refdf["p"] < args.topthresh
    refdf_test = refdf.loc[istestsnp, :]
    refdf_null = refdf.loc[~istestsnp, :]
    print(" {0} test snps with p-value < {1} found".format(refdf_test.shape[0],
                                                           args.topthresh))

    # Exclude test SNPs based on their proximity to each other
    # refdf_test = exclude_test_snps(refdf_test)

    # Give warning if nullsize is greater than pool of null snps
    if args.topthresh > refdf_null.shape[0]:
        message = ("Warning: --nullsize is greater than the total number of"
                   "null snps.")
        print(message)
    
    # Calc how many null SNPs to sample per test snp
    nullpertest = int(args.nullsize / refdf_test.shape[0]) + 1
    print(" {0} null snps will be sampled per test snp".format(nullpertest))
    
    # For each test snp sample null snps with a similar profile
    sampled_df = sample_null_dist(refdf_test, refdf_null, nullpertest)

    # Output test snps and sample dataframes
    print("Saving copy of sampled SNPs and test SNPs...")
    outname = args.out + "_test-snps.tsv"
    refdf_test.to_csv(outname, sep="\t", header=True, index=False)
    outname = args.out + "_null-dist-snps.tsv"
    sampled_df.to_csv(outname, sep="\t", header=True, index=False)

    # Make a plot of the null
    # from ggplot import *
    # outname = args.out + "_null-hist.png"
    # p = ggplot(sampled_df, aes(x='stat')) + geom_histogram()
    # ggsave(p, outname)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load second phenotype and compare test SNPs to null distribution
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Load data
    print("Loading test dataset...")
    testdf = load_sumstats(args.teststats)

    # Get test snps
    testdf_test = testdf.loc[refdf_test.index, :]

    # Get np array of null dist stats
    null_dist = np.array(sampled_df["stat"])

    # For each test SNP compare stat is null distribution
    results = compare_tests_to_null(testdf_test, null_dist)

    # Write results to file
    print("Writing results...")
    header = ["snp", "chr", "pos", "maf", "stat", "enrich_pval", "null_size"]
    outname = args.out + "_enrichment.tsv"
    with open(outname, "w") as out_h:
        # Write header
        out_h.write("\t".join(header) + "\n")
        # Write each line
        for line in sorted(results, key=itemgetter(5)):
            line = [str(x) for x in line]
            out_h.write("\t".join(line) + "\n")

    # Close log
    if args.log == True:
        log_h.close()

    return 0

# def exclude_test_snps(refdf_test):
#     """ Removes SNPs that are close to each other.
#     """
#     # Sort by P-value so that we are always keeping the 
#     refdf_test = refdf_test.sort("p")
#     print(refdf_test.head())
#     sys.exit()

def compare_tests_to_null(testdf_test, null_dist):
    """ For each SNP in the test dataset, compare it to the null and produce
        results table
    """
    print("Comparing test SNPs to null distribution...")
    results = []
    null_size = len(null_dist)
    for i in range(testdf_test.shape[0]):

        # Get test snp info
        test_row = testdf_test.iloc[i, ]
        test_name = test_row["snp"]
        test_maf = test_row["maf"]
        test_chrom = test_row["chr"]
        test_bp = test_row["bp"]
        test_stat = test_row["stat"]

        # Count how many times test stat is more extreme than null dist
        if args.test == "upper":
            counts = test_stat > null_dist
        elif args.test == "lower":
            counts = test_stat < null_dist
        elif args.test == "two":
            counts = abs(test_stat) < np.absolute(null_dist)
        # Covert to pval
        pval = float(sum(counts)) / null_size
        if args.test == "two":
            pval = pval / 2
        pval = 1 - pval

        # Make output
        outline = [test_name,
                   test_chrom,
                   test_bp,
                   test_maf,
                   test_stat,
                   pval,
                   null_size]
        results.append(outline)

    return results

def sample_null_dist(refdf_test, refdf_null, nullpertest):
    """ For each SNP in test set, finds SNPs with similar profile in null set
        and samples them to produce a null distribution.
    """
    print("Sampling distribution of null for each test snp...")
    firstiter = True
    for i in range(refdf_test.shape[0]):
        
        # Get profile of test snp
        test_row = refdf_test.iloc[i, ]
        test_name = test_row["snp"]
        test_maf = test_row["maf"]
        test_chrom = test_row["chr"]
        test_bp = test_row["bp"]
        test_p = -np.log10(test_row["p"])
        test_l2 = test_row["l2"]

        # Increase test interval by posrange
        test_interval = test_row["interval"]
        test_interval_exp = GenomicInterval(str(test_chrom),
                                            test_bp - args.posrange,
                                            test_bp + args.posrange,
                                            ".")
        # Print info
        # print("-Test snp: {0}".format(test_name))
        # print("  -log10(p-val) = {0}".format(test_p))
        # print("  maf = {0}".format(test_maf))
        # print("  pos = {0}:{1}".format(test_chrom, test_bp))

        # Copy df of candidate SNPs
        candidate_df = refdf_null.copy()
        print(candidate_df.shape)

        # Keep SNPs where maf in range
        mafinrange = candidate_df.loc[:, "maf"].apply(maf_in_range,
                                                      test_maf=test_maf,
                                                      mafrange=args.mafrange)
        candidate_df = candidate_df.loc[mafinrange, :]

        # Keep SNPs where l2 score is in range
        l2inrange = candidate_df.loc[:, "l2"].apply(maf_in_range,
                                                    test_maf=test_l2,
                                                    mafrange=args.ldscrange)
        candidate_df = candidate_df.loc[l2inrange, :]

        # Keep SNPs where interval outside test_interval_exp
        iscontained = candidate_df.loc[:, "interval"].apply(
            lambda x: x.is_contained_in(test_interval_exp))
        candidate_df = candidate_df.loc[~iscontained, :]
        print(candidate_df.shape)

        # Warn if sample n is larger than candidate_df
        if nullpertest > candidate_df.shape[0]:
            message = (" Warning: nullpertest is greater than the number of "
                       "SNPs matching test SNP profile.")
            print(message)
        
        # Warn and skip if no snps in candidate
        if candidate_df.shape[0] == 0:
            message =  " Warning: No null candidates matching {0} profile".format(test_name)
            print(message)
            continue

        # Take random sample of candidate_df snps
        rand_idx = np.random.choice(candidate_df.shape[0], nullpertest,
                                    replace=True)
        candidate_sample_df = candidate_df.iloc[rand_idx, :]
        candidate_sample_df["test_snp"] = test_name

        # Concatenate to overall sample
        if firstiter:
            sampled_df = candidate_sample_df
            firstiter = False
        else:
            sampled_df = sampled_df.append(candidate_sample_df)

    return sampled_df

def maf_in_range(cand_maf, test_maf, mafrange):
    """ Given pandas row containing maf, checks if maf is within the interval
        defined in --mafrange.
    """
    if cand_maf > (test_maf - mafrange) and cand_maf < (test_maf + mafrange):
        return True
    else:
        return False

def load_merge_ldscores(refdf, ldscfiles):
    """ Parses each of the LD score files and gets LD scores for SNPs.
    """
    missingvals = set(args.missing)
    snpstoget = set(refdf.index)

    # Iterate over ldsc files
    rows = []
    for filen in ldscfiles:
        with get_handle(filen, "r") as in_h:
            # Skip header
            in_h.readline()
            # Get parts
            for line in in_h:
                chrom, snp, bp, l2 = line.rstrip().split("\t")
                # If SNP is missing, make new snp name
                if snp in missingvals:
                    snp = ":".join(chrom, bp)
                # Add to rows if in refdf
                if snp in snpstoget:
                    rows.append([snp, float(l2)])
    
    # Make a dict
    l2_dict = dict(rows)
    keys = set(l2_dict.keys())
    missing = snpstoget.difference(keys)

    # Warn about missing LD scores
    if len(keys) < len(snpstoget):
        outname = args.out + "_missing-LD-score.tsv"
        message = [" Warning: {0} snps in {1}.".format(len(snpstoget),
                                                      args.refstats),
                   " Warning: {0} snps have LD scores.".format(len(keys)),
                   " Warning: {0} snps missing LD scores.".format(len(missing)),
                   " Warning: Snps with missing LD scores will be dropped. ",
                   " Warning: Missing LD score SNPs will be written to file."]
        print("\n".join(message))
        with open(outname, "w") as out_h:
            for snp in missing:
                out_h.write(snp + "\n")

    # Merge l2 scores to refdf
    refdf["l2"] = refdf.loc[:, "snp"].apply(get_l2, l2_dict=l2_dict)
    
    # Drop rows with missing LD scores
    isnan = pd.isnull(refdf["l2"])
    refdf = refdf.loc[~isnan, :]

    return refdf

def get_l2(snp, l2_dict):
    """ Returns l2_dict[snp], else nan.
    """
    try:
        return l2_dict[snp]
    except KeyError:
        return np.nan

def get_handle(filen, rw="r"):
    """ Returns gzip handle if ext is gz.
    """
    ext = filen.rsplit(".", 1)[-1]
    if ext == "gz":
        return gzip.open(filen, rw)
    else:
        return open(filen, rw)

def load_sumstats(filen):
    """ Load GWAS summary stats. Returns pandas df of required columns.
    """
    # Check if compressed
    ext = filen.rsplit(".", 1)[-1]
    if ext == "gz":
        compr = "gzip"
    elif ext == "bz2":
        compr = "bz2"
    else:
        compr = None
    
    # Load dataframe
    df = pd.read_csv(filen, sep=args.sep, compression=compr,
                     na_values=args.missing)
    
    # Check that required columns are present
    cols = [args.chrcol, args.bpcol, args.snpcol, args.testcol, args.mafcol,
            args.pcol]
    for col in cols:
        if not col in df.columns:
            message = "{0} not found in {1}".format(col, filen)
            sys.exit(message)
    
    # Take only required columns
    df = df.loc[:, cols]
    df.columns = ["chr", "bp", "snp", "stat", "maf", "p"]
    
    # If snp col is NaN then create an ID from chr:pos
    isnan = pd.isnull(df["snp"])
    newname = df.loc[isnan, ["chr", "bp"]].apply(concat_strs, axis=1)
    df.loc[isnan, "snp"] = newname
    
    # Give warning if there are duplicate IDs
    isdup = df.duplicated(subset="snp")
    if np.sum(isdup) > 0:
        dup_id = ", ".join(list(df.loc[isdup, "snp"]))
        message = (" Warning: Duplicate SNP IDs found: " + dup_id + 
                   "\n Warning: Dropping one of each duplicate.")
        print(message)
        df = df.drop_duplicates(subset="snp")

    
    # Set snp as index
    df.index = df["snp"]

    # Create genomic intervals
    df["interval"] = df.apply(make_GI, axis=1)

    return df

def make_GI(row):
    """ Given a pandas row, will create a genomic interval using chr and bp.
    """
    gi = GenomicInterval(str(row["chr"]), row["bp"], row["bp"] + 1, ".")
    return(gi)

def concat_strs(row, sep=":"):
    """ Concat strs together.
    """
    strs = [str(x) for x in list(row)]
    return sep.join(strs)

def parse_arguments():
    """ Parses command line arguments.
    """
    # Create top level parser.
    parser = argparse.ArgumentParser(description="GWAS stats enrichment")

    # Add files
    parser.add_argument('refstats', metavar="<reference sum stats>",
        help=('GWAS summary stats for selecting top SNPS and creating null.'),
        type=str)
    parser.add_argument('teststats', metavar="<test sum stats>",
        help=('GWAS summary stats to test for enrichment.'),
        type=str)
    parser.add_argument('ldscstats', metavar="<ld score file>",
        help=('List of LD score files to be merged. LD scores can be downloaded'
              ' from https://github.com/bulik/ldsc/#where-can-i-get-ld-scores'),
        type=str,
        nargs="+")
    parser.add_argument('--out', metavar="<out prefix>",
        help=('Output prefix. (default: out)'),
        type=str,
        default="out")
    parser.add_argument('--log',
        help=('Store log rather than print to stdout.'),
        action='store_true',
        default=False)
    # Add test options
    parser.add_argument('--test', choices=["lower", "upper", "two"],
        help=('Select one-tailed [upper/lower] or two-tailed test. '
              '(default: lower)'),
        type=str,
        default="lower")
    parser.add_argument('--topthresh', metavar="<float>",
        help=('SNPs below this p-value threshold will be used for test. '
              '(default: 1e-5)'),
        type=float,
        default=1e-5)
    parser.add_argument('--nullsize', metavar="<int>",
        help=('Number of SNPs to sample for null distribution.'
             ' (default: 100,000)'),
        type=int,
        default=100000)
    
    # Null profile options
    parser.add_argument('--mafrange', metavar="<float>",
        help=('Max MAF deviation from reference SNP for inclusion in null.'
             ' (default: 0.02)'),
        type=float,
        default=0.02)
    parser.add_argument('--posrange', metavar="<int>",
        help=('Null sample SNPs must be further than this distance from'
              ' reference SNP. (default: 1,000,000)'),
        type=int,
        default=1000000)
    parser.add_argument('--ldscrange', metavar="<int>",
        help=('Max LD score deviation from reference SNP for inclusion in null.'
              ' reference SNP. (default: 2)'),
        type=float,
        default=15) ################################# Needs changing ~~~~~~~~~~~~
    
    # Add column name args
    parser.add_argument('--testcol', metavar="<str>",
        help='Column with test-statistic. (default: frequentist_add_pvalue)',
        type=str,
        default="frequentist_add_pvalue")
    parser.add_argument('--snpcol', metavar="<str>",
        help='Column with SNP id. (default: rsid)',
        type=str,
        default="rsid")
    parser.add_argument('--bpcol', metavar="<str>",
        help='Column with SNP position. (default: position)',
        type=str,
        default="position")
    parser.add_argument('--chrcol', metavar="<str>",
        help='Column with chromosome number. (default: alternate_ids)',
        type=str,
        default="alternate_ids")
    parser.add_argument('--pcol', metavar="<str>",
        help='Column with chromosome number. (default: frequentist_add_pvalue)',
        type=str,
        default="frequentist_add_pvalue")
    parser.add_argument('--mafcol', metavar="<str>",
        help='Column with MAF. (default: all_maf)',
        type=str,
        default="all_maf")
    parser.add_argument('--sep', metavar="<str>",
        help='Column separater. (default: tab)',
        type=str,
        default="\t")
    parser.add_argument('--missing', metavar="<str>",
        help='List of values to use as missing. (default: . "" NA)',
        type=str,
        nargs="+",
        default=[".", "", "NA"])
    
    # Parse the arguments
    args = parser.parse_args()

    # Parse the arguments
    return args


if __name__ == '__main__':

    main()
