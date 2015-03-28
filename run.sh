ldsc=data/ldscores/2*.l2.ldscore.gz
t1=data/iq_subtest_sum_ex4sd_unrel1_warped_f8_chr20-21.tsv.gz
t2=data/visacutiy_best_acuity_ex4sd_unrel1_warped_f7_chr20-21.tsv.gz
ldmap=data/ldmap/unrelated-children_maf0.01_info0.8.ld.gz
intronfeat=make-features/make-exons/snpinintron.txt
exonfeat=make-features/make-exons/snpinexon.txt
out=results/test1

time python gwas-ct-enrich.py $t1 $t2 $ldsc $ldmap \
--nullsize 10000 --out $out --maxr2 0.5 --mafrange 0.05 \
--features $intronfeat $exonfeat
