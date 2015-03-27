ldsc=data/ldscores/2*.l2.ldscore.gz
t1=data/iq_subtest_sum_ex4sd_unrel1_warped_f8_chr20-21.tsv.gz
t2=data/visacutiy_best_acuity_ex4sd_unrel1_warped_f7_chr20-21.tsv.gz
ldmap=data/ldmap/unrelated-children_maf0.01_info0.8.ld.gz
out=results/test1

time python gwas-ct-enrich.py $t1 $t2 $ldsc --nullsize 10000 --log --out $out
