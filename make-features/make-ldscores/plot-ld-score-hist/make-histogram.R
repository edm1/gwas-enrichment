library("ggplot2")

# Load data
filen = "/home/em13383/gwas-enrichment/make-features/make-ldscores/plot-ld-score-hist/sample-ld-scores.txt.gz"
df = read.table(gzfile(filen))
head(df)

# Plot histo
ggplot(df, aes(x=V1)) + geom_histogram()
