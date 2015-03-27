ldsc=ldscores/2*.l2.ldscore.gz

time python gwas-ct-enrich.py \
data/iq_subtest_sum_ex4sd_unrel1_warped_f8_chr20-21.tsv.gz \
data/visacutiy_best_acuity_ex4sd_unrel1_warped_f7_chr20-21.tsv.gz \
$ldsc \
--nullsize 10000 --log
