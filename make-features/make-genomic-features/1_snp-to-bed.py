#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# Reads stdin and writes snps in bed format to stdout
#

import sys

def main():

    # Set column names
    chr_col = "alternate_ids"
    pos_col = "position"
    snp_col = "rsid"
    # Set column separator
    sep = "\t"
    # Set possible missing values
    missing_vals = set(["", ".", "NA", "NaN"])

    # Get header and header positions
    header = sys.stdin.readline().rstrip().split("\t")
    chr_index = header.index(chr_col)
    pos_index = header.index(pos_col)
    snp_index = header.index(snp_col)

    # For each line convert to bed format
    for line in sys.stdin:
        # Get parts
        parts = line.rstrip().split("\t")
        chrom = parts[chr_index]
        pos = int(parts[pos_index])
        snp = parts[snp_index]
        # Make snp name if missing
        if snp in missing_vals:
            snp = "{0}:{1}".format(chrom, pos)
        # Print bed
        outline = [chrom, pos, pos + 1, snp, ".", "."]
        outline = [str(x) for x in outline]
        print("\t".join(outline))

if __name__ == '__main__':

    main()
