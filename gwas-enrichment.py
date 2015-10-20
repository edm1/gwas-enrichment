#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# The MIT License (MIT)
#
# Copyright (c) 2015 Edward Mountjoy
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import sys
import os
import argparse
import pandas as pd
import numpy as np
from HTSeq import GenomicInterval
import gzip
from operator import itemgetter
from statsmodels.sandbox.stats.multicomp import multipletests

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
    print_args(args)
    if args.features == None:
        args.features = []

    # Set up log
    if args.log == True:
        outname = outname = args.out + "_log.txt"
        log_h = open(outname, "w")
        sys.stdout = log_h

    # Set seed
    if not args.randseed:
        np.random.seed(123)

    #
    # Load input
    #

    # Load reference df
    print("Loading reference sum stats...")
    refdf = load_sumstats(args.refstats)
    print(" done")

    # Load additional genomic features
    print("Loading addtional genomic features...")
    feat_cols = []
    for infile in args.features:
        refdf, colname = bind_features(infile, refdf)
        feat_cols.append(colname)

    # Merge LD scores to reference df
    print("Loading LD scores...")
    refdf = load_merge_ldscores(refdf, args.ldscstats)
    print(" done")

    # Split refdf into test snps vs null snps
    print("Splitting test and null snps...")
    istestsnp = refdf["p"] < args.topthresh
    refdf_test = refdf.loc[istestsnp, :]
    refdf_null = refdf.loc[~istestsnp, :]
    print(" {0} test snps with p-value < {1} found".format(refdf_test.shape[0],
                                                           args.topthresh))
    if refdf_test.shape[0] == 0:
        sys.exit("Exiting: No test SNPs found.")

    #
    # Exclude test snps that are similar to each other
    #

    # Exclude test SNPs that are close to each other
    refdf_posexcluded = None
    if args.testposrange:
        print("Excluding test snps close to each other in genome...")
        refdf_test, refdf_posexcluded = exclude_test_snps_range(refdf_test)
        print(" {0} test snps excluded".format(refdf_posexcluded.shape[0]))
        print(" {0} test snps kept".format(refdf_test.shape[0]))

    # Exclude test SNPs that are in high LD
    refdf_r2excluded = None
    if args.ldmap:
        print("Excluding high r2 test snps...")
        refdf_test, refdf_r2excluded = exclude_test_snps_LD(refdf_test)
        print(" {0} test snps excluded".format(refdf_r2excluded.shape[0]))
        print(" {0} test snps kept".format(refdf_test.shape[0]))

    # Give warning if nullsize is greater than pool of null snps
    if args.topthresh > refdf_null.shape[0]:
        message = ("Warning: --nullsize is greater than the total number of"
                   "null snps.")
        print(message)

    #
    # Sample null distribution match each of the test SNPs
    #

    # For each test snp sample null snps with a similar profile
    sampled_df = sample_null_dist(refdf_test, refdf_null, feat_cols)

    # Output test snps and sample dataframes
    print("\nSaving copy of sampled SNPs and test SNPs...")
    outname = args.out + "_test-snps.tsv"
    refdf_test.to_csv(outname, sep="\t", header=True, index=False)
    if not refdf_posexcluded is None:
        outname = args.out + "_pos-excluded-test-snps.tsv"
        refdf_posexcluded.to_csv(outname, sep="\t", header=True, index=False)
    if not refdf_r2excluded is None:
        outname = args.out + "_r2-excluded-test-snps.tsv"
        refdf_r2excluded.to_csv(outname, sep="\t", header=True, index=False)
    if not args.supress:
        outname = args.out + "_null-dist-snps.tsv"
        sampled_df.to_csv(outname, sep="\t", header=True, index=False)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load second phenotype and compare test SNPs to null distribution
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Make empty results folder
    header = ["ref_file", "test_file", "snp", "chr", "pos", "maf", "stat",
              "enrich_pval", "null_size"]
    results_df = pd.DataFrame(columns=header)
    
    for testfile in args.teststats:

        # Name to give ref and test set
        uniq_names = uniq_name_path([args.refstats, testfile])
        # Skip if None (ref and test were the same file)
        if uniq_names is None:
            continue
        refname, testname = uniq_names

        # Load data
        print(" Loading test dataset {0}...".format(testname))
        testdf = load_sumstats(testfile)
        # Get test snps
        testdf_test = testdf.loc[refdf_test.index, :]

        # For each test SNP compare stat is null distribution
        results = compare_tests_to_null(testdf_test, sampled_df)

        # Prepare results file
        print(" Appending results...")

        # Add ref and testname to each result row
        results = [[refname, testname] + x for x in results]
        
        # Append rows to results df
        temp_df = pd.DataFrame(results, columns=header)
        results_df = results_df.append(temp_df)

    print("Writing output...")
    # Sort by pval
    results_df = results_df.sort("enrich_pval")
    # Do multiple testing corrections
    mt_correct = multipletests(results_df["enrich_pval"], method="fdr_bh",
                               alpha=0.05)[1]
    results_df["enrich_adjusted_pval"] = mt_correct
    # Write file
    outname = args.out + "_enrichment.tsv"
    results_df.to_csv(outname, sep="\t", header=True, index=False)

    # Close log
    if args.log == True:
        log_h.close()

    return 0

def exclude_test_snps_range(refdf_test):
    """ Test snps are sorted by P-value then snps further down the list are
        excluded if within genomc range: option --testposrange
    """
    # Sort by P-value so that we are always keeping the best
    refdf_test = refdf_test.sort("p")

    # Find exclusions
    print(" making exclusions...")
    toexclude = set([])
    for i in range(refdf_test.shape[0]):
        # Get test snp info
        row = refdf_test.iloc[i, :]
        snp = row["snp"]
        # Skip if already excluded
        if snp in toexclude:
            continue
        # Expand genomic interval
        test_chrom = row["chr"]
        test_bp = row["bp"]
        test_interval_exp = GenomicInterval(str(test_chrom),
                                            test_bp - args.testposrange * 1000,
                                            test_bp + args.testposrange * 1000,
                                            ".")
        # Exclude current snp from analysis
        refdf_test_loo = refdf_test.loc[~refdf_test.index.isin([snp]), :]
        # Find other test SNPs in interval
        iscontained = refdf_test_loo.loc[:, "interval"].apply(
            lambda x: x.is_contained_in(test_interval_exp))
        containedindex = refdf_test_loo.loc[iscontained, ].index
        toexclude = toexclude.union(list(containedindex))

    # Make exclusions
    refdf_excluded = refdf_test.loc[refdf_test.index.isin(toexclude), :]
    refdf_keep = refdf_test.loc[~refdf_test.index.isin(toexclude), :]

    return refdf_keep, refdf_excluded

def print_args(args):
    """ Print args to log.
    """
    print("Run arguments:")
    arg_dict = vars(args)
    for k in arg_dict:
        print(" {0}: {1}".format(k, arg_dict[k]))
    print("")

def bind_features(infile, refdf):
    """ Binds a column stating which each SNP is part of the feature set or not.
    """
    colname = os.path.split(infile)[1].rsplit(".", 1)[0]

    # Load a set of snps contained in feature
    isfeature = set([])
    with get_handle(infile, "r") as in_h:
        for line in in_h:
            snp = line.rstrip()
            isfeature.add(snp)

    # Make new column
    refdf[colname] = refdf["snp"].apply(lambda x: x in isfeature)
    
    print(" " + colname)
    print(refdf[colname].value_counts())

    return refdf, colname

def exclude_test_snps_LD(refdf_test):
    """ Removes test SNPs that are in clode LD with each other.
    """
    # Sort by P-value so that we are always keeping the best
    refdf_test = refdf_test.sort("p")

    # Load LD information
    print(" loading LD map...")
    ld_df = load_ldmap(args.ldmap, refdf_test["snp"])
    if ld_df.shape[0] == 0:
        return refdf_test

    # Make into a dict
    ld_dict = {}
    pairs = zip(ld_df["snpA"], ld_df["snpB"])
    for snpA, snpB in pairs:
        # Add snpB to snpA
        try:
            ld_dict[snpA].add(snpB)
        except KeyError:
            ld_dict[snpA] = set([snpB])
        # Recipricol
        try:
            ld_dict[snpB].add(snpA)
        except KeyError:
            ld_dict[snpB] = set([snpA])

    # Find exclusions
    print(" making exclusions...")
    toexclude = set([])
    for i in range(refdf_test.shape[0]):
        # Get test snp info
        row = refdf_test.iloc[i, :]
        snp = row["snp"]
        # Skip if already excluded or not in ld
        if snp in toexclude or snp not in ld_dict:
            continue
        # Exclude LD snps
        toexclude = toexclude.union(ld_dict[snp])

    # Make exclusions
    refdf_excluded = refdf_test.loc[refdf_test.index.isin(toexclude), :]
    refdf_keep = refdf_test.loc[~refdf_test.index.isin(toexclude), :]

    return refdf_keep, refdf_excluded
    

def load_ldmap(filen, snplist):
    """ Parses plink ld file to find snps in LD with those provided in snplist.
        Returns pandas df.
    """
    snpset = set(snplist)
    missingvals = set(args.missing)
    data = []
    with get_handle(filen) as in_h:
        # Get header
        header = in_h.readline()
        # Iterate over lines
        for line in in_h:
            line = line.rstrip().split()
            chrA, bpA, snpA, chrB, bpB, snpB, r2 = line
            # If SNP is missing, make new snp name
            if snpA in missingvals:
                snpA = ":".join(chromA, bpA)
            if snpB in missingvals:
                snpB = ":".join(chromB, bpB)
            # Add to rows if in refdf
            if snpA in snpset and snpB in snpset:
                data.append([snpA, snpB, float(r2)])
    # Make df
    df = pd.DataFrame(data, columns=["snpA", "snpB", "r2"])
    df = df.loc[df["r2"] > args.maxr2]
    return df

def compare_tests_to_null(testdf_test, sampled_df):
    """ For each SNP in the test dataset, compare it to the null and produce
        results table
    """
    print("Comparing test SNPs to null distribution...")
    results = []
    for i in range(testdf_test.shape[0]):

        # Get test snp info
        test_row = testdf_test.iloc[i, ]
        test_name = test_row["snp"]
        test_maf = test_row["maf"]
        test_chrom = test_row["chr"]
        test_bp = test_row["bp"]
        test_stat = test_row["stat"]

        # Get null dist
        null_dist_df = sampled_df.loc[sampled_df["test_snp"] == test_name, :]
        null_dist = np.array(null_dist_df["stat"])
        null_size = len(null_dist)
        if null_size == 0: ############# TODO make this better
            continue

        # Count how many times test stat is more extreme than null dist
        if args.test == "upper":
            counts = test_stat > null_dist
        elif args.test == "lower":
            counts = test_stat < null_dist
        elif args.test == "two":
            counts = abs(test_stat) < np.absolute(null_dist)
        # Covert to pval
        pval = (float(sum(counts)) + 1) / (null_size + 1)
        if args.test == "two":
            pval = pval / 2
        pval = 1 - pval
        # # Disallow Pvalue of 0.0
        # if pval == 0.0:
        #     pval = 1.0 / null_size

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

def sample_null_dist(refdf_test, refdf_null, feat_cols):
    """ For each SNP in test set, finds SNPs with similar profile in null set
        and samples them to produce a null distribution.
    """
    print("Sampling distribution of null for each test snp...")
    firstiter = True
    for i in range(refdf_test.shape[0]):

        nullsamplesize = args.nullsize
        
        # Get profile of test snp
        test_row = refdf_test.iloc[i, ]
        test_name = test_row["snp"]
        test_maf = test_row["maf"]
        test_chrom = test_row["chr"]
        test_bp = test_row["bp"]
        test_p = -np.log10(test_row["p"])
        test_l2bin = test_row["l2_bin"]
        print("\n sampling for snp: {0}".format(test_name))

        # Increase test interval by posrange
        test_interval = test_row["interval"]
        test_interval_exp = GenomicInterval(str(test_chrom),
                                            test_bp - args.posrange * 1000,
                                            test_bp + args.posrange * 1000,
                                            ".")

        # Copy df of candidate SNPs
        candidate_df = refdf_null.copy()
        print(" total candidates: {0}".format(candidate_df.shape[0]))

        # Keep SNPs where maf in range
        mafinrange = candidate_df.loc[:, "maf"].apply(maf_in_range,
                                                      test_maf=test_maf,
                                                      mafrange=args.mafrange)
        candidate_df = candidate_df.loc[mafinrange, :]
        print(" after maf match: {0}".format(candidate_df.shape[0]))

        # Keep SNPs where l2 score is in range
        l2same = candidate_df["l2_bin"] == test_l2bin
        candidate_df = candidate_df.loc[l2same, :]
        print(" after l2 match: {0}".format(candidate_df.shape[0]))

        # Keep SNPs where interval outside test_interval_exp
        iscontained = candidate_df.loc[:, "interval"].apply(
            lambda x: x.is_contained_in(test_interval_exp))
        candidate_df = candidate_df.loc[~iscontained, :]
        print(" after position exclusion: {0}".format(candidate_df.shape[0]))

        # For each feature column, keep snps that match test feature
        for col in feat_cols:
            test_featbool = test_row[col]
            matches = candidate_df.loc[:, col].apply(
                                                 lambda x: x == test_featbool)
            candidate_df = candidate_df.loc[matches, :]
            print(" after {0} feature match: {1}".format(col,
                                                         candidate_df.shape[0]))

        # Warnings
        cand_size = candidate_df.shape[0]
        if cand_size == 0:
            # Warn and skip if no snps in candidate
            message =  (" Warning: No null candidates matching {0} "
                        "profile".format(test_name))
            print(message)
            # Skip
            continue
        elif cand_size < nullsamplesize:
            # Warn if sample n > candidate_df n and set n to cand size

            message = [" Warning: test snp {0}".format(test_name),
            " Warning: --nullsize {0} > candidate snps {1}".format(args.nullsize,
                                                                   cand_size),
            " Warning: setting --nullsize to candidate size"]
            print("\n".join(message))
            nullsamplesize = candidate_df.shape[0]

        # Take random sample of candidate_df snps
        rand_idx = np.random.choice(candidate_df.shape[0], nullsamplesize,
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
                   " Warning: snps with missing LD scores will be dropped. "]
        print("\n".join(message))
        with open(outname, "w") as out_h:
            for snp in missing:
                out_h.write(snp + "\n")

    # Merge l2 scores to refdf
    refdf["l2"] = refdf.loc[:, "snp"].apply(get_l2, l2_dict=l2_dict)

    # Make percentile bins
    bins = np.arange(0, 100, 100 / args.ldscbins)
    perc = [0] + list(np.percentile(refdf["l2"], bins)) + [max(refdf["l2"]) + 1]
    refdf["l2_bin"] = pd.cut(refdf["l2"], perc, include_lowest=True)
    
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

    # Add input files
    parser.add_argument('refstats', metavar="<reference sum stats>",
        help=('GWAS summary stats for selecting top SNPS and creating null.'),
        type=str)
    parser.add_argument('teststats', metavar="<test sum stats>",
        help=('GWAS summary stat files to test for enrichment. Multiple test '
              'files can be tested at once against the same reference null.'),
        type=str,
        nargs="+")
    parser.add_argument('--ldscstats', metavar="<ld score file>",
        help=('List of LD score files to be merged. LD scores can be generated'
              ' using ldscr from https://github.com/bulik/ldsc/'),
        type=str,
        nargs="+",
        required=True)
    parser.add_argument('--features', metavar="<txt file>",
        help=('Genomic features to match profile on. Each file should be a'
              ' list of SNP names which share a feature.'),
        nargs='*')

    # Add output options
    parser.add_argument('--out', metavar="<out prefix>",
        help=('Output prefix. (default: out)'),
        type=str,
        default="out")
    parser.add_argument('--supress',
        help=("Suppress null dist output which is very large. (default: False)"),
        action="store_true",
        default=False)
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
              '(default: 1e-6)'),
        type=float,
        default=1e-6)
    parser.add_argument('--nullsize', metavar="<int>",
        help=('Number of SNPs to sample for null distribution.'
             ' (default: 100,000)'),
        type=int,
        default=100000)

    # Exclude similar test snps
    parser.add_argument('--testposrange', metavar="<int kb>",
        help=('Test SNPs are excluded if there is another test SNP with '
              'higher p-value within this genomic range. Useful to reduce '
              'multiple testing burden. (default: None)'),
        type=int,
        default=None)
    parser.add_argument('--ldmap', metavar="<plink ld file>",
        help=('LD map output from plink --r2. Used to exclude test SNPs in '
              'high LD. Useful to reduce multiple testing burden.'),
        type=str)
    parser.add_argument('--maxr2', metavar="<float>",
        help=('Test SNPs with r2 higher than this will be excluded. '
              'Only used in conjunction with --ldmap. (default: 0.7)'),
        type=float,
        default=0.7)
    
    # Null profile options
    parser.add_argument('--mafrange', metavar="<float>",
        help=('Max MAF deviation from reference SNP for inclusion in null.'
             ' (default: 0.02)'),
        type=float,
        default=0.02)
    parser.add_argument('--posrange',metavar="<int kb>",
        help=('Null sample SNPs must be further than this distance from'
              ' reference SNP. (default: 1000 kb)'),
        type=int,
        default=1000)
    parser.add_argument('--ldscbins', metavar="<int>",
        help=('Number of bins to make for LD scores. (default: 10)'),
        type=float,
        default=10)
    parser.add_argument('--randseed',
        help=('Use a random seed instead of 123.'),
        action="store_true",
        default=False)
    
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
        help=('Column with p-value significance. '
              '(default: frequentist_add_pvalue)'),
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

def uniq_name_path(path_list):
    """ Returns the shortest path suffix that identifies a file uniquely. (I'm
        sure there must be an easier way to do this)
    """
    # Convert paths to absolutes and normalise
    normpath_list = [os.path.normpath(os.path.abspath(x)) for x in path_list]
    # Check that no two paths are identical
    if len(set(normpath_list)) != len(path_list):
        return None
    # Find shortest suffix
    normpath_lol = [x.split(os.sep) for x in normpath_list]
    unzipped = zip(*normpath_lol)
    for i in xrange(len(unzipped) - 1, -1, -1):
        j = i
        if len(set(unzipped[i])) == len(path_list):
            break
    shortestsuffix = zip(*unzipped[j:])
    # Coverts back to path
    uniqpaths = [reduce(os.path.join, x) for x in shortestsuffix]
    return uniqpaths

if __name__ == '__main__':

    main()
