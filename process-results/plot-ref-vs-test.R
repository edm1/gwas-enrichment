library("ggplot2")
library("RColorBrewer")
library("reshape2")

if (interactive()) {
  setwd("/home/em13383/gwas-enrichment/process-results/")
  test_file = "data/pos-exclusion/visacutiy_best_acuity_V_iq_subtest_sum_enrichment.tsv"
  #test_file = "data/no-exclusion/res_visac-v-iq_enrichment.tsv"
  ref_file = "data/pos-exclusion/visacutiy_best_acuity_V_iq_subtest_sum_test-snps.tsv"
  #ref_file = "data/no-exclusion/res_visac-v-iq_test-snps.tsv"
  title = "test_title"
  outpref = "test_plot"
} else {
  args = commandArgs(trailingOnly=T)
  test_file = args[1]
  ref_file = args[2]
  title = args[3]
  outfile = args[4]
}

# Load test results
df.test = read.table(test_file, sep="\t", header=T)
colnames(df.test) = paste0("test_", colnames(df.test))
head(df.test)

# Load reference info
df.ref = read.table(ref_file, sep="\t", header=T)
colnames(df.ref) = paste0("ref_", colnames(df.ref))
head(df.ref)

# Merge
colnames(df.test)
colnames(df.ref)
df.data = merge(df.test, df.ref, by.x="test_snp", by.y="ref_snp")
df.data = df.data[order(df.data$test_enrich_pval), ]
head(df.data)

# Plot scatter
ggplot(df.data, aes(x=-log10(ref_p), y=-log10(test_enrich_pval), colour=factor(test_chr))) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), colour="blue", linetype="dashed", alpha=0.7) +
  geom_vline(xintercept=-log10(5e-8), colour="red", linetype="dashed", alpha=0.7) +
  ggtitle(title)

ggsave(filename=paste0(outpref, "_points.png"),
       height=10,
       width=12)

#
# Make a manhattan plot with two Y axes
#

# Melt data
head(df.data)
df.m = melt(df.data, c("test_chr", "test_pos", "test_snp"), c("ref_p", "test_enrich_pval"))
levels(df.m$variable) = c("Association p-value in reference", "SNP enrichment in test set")

# Sort by chr then pos then add index
df.m = df.m[order(df.m$test_chr, df.m$test_pos), ]
df.m["index"] = 1:nrow(df.m)

# Make a colour palette
nchr = length(unique(df.m$test_chr))
colpal = brewer.pal(nchr, "Dark2")
col = rep(colpal, nchr)

# Plot
ggplot(df.m, aes(x=index, y=-log10(value))) +
  # Add points, set colours and remove the legend
  geom_point(aes(colour=col[test_chr])) + 
  facet_grid(variable ~ ., scales="free_y") +
  scale_colour_manual(values=colpal) +
  theme(legend.position="none") +
  # Change labels and title
  labs(x="Index", y="-log10(P-value)") +
  ggtitle(title) +
  # Change theme
  theme(axis.text.x=element_text(colour="black", size=12),
        axis.text.y=element_text(colour="black", size=12),
        plot.title=element_text(face="bold"))

# Save plot
ggsave(filename=paste0(outpref, "_manhat.png"),
       width = 12,
       height = 6,
       dpi = 300)
