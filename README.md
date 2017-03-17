# *Coenagrion puella* differential expression

- [Transcript quantification with `kallisto`](#transcript-quantification-with-kallisto)
- [Read in data and create sample descriptions](#read-in-data-and-create-sample-descriptions)
- [Run `DESeq` with different dispersion fits](#run-deseq-with-different-dispersion-fits)
  - [Choose a dispersion fit](#choose-a-dispersion-fit)
- [Regularized log transformation and PCA](#regularized-log-transformation-and-pca)

## Transcript quantification with `kallisto`

```sh
bunzip2 c_puella.fasta.bz2
# prep reference
align_and_estimate_abundance.pl --transcripts c_puella.fasta --est_method kallisto --trinity_mode --prep_reference --output_dir ./
# align library
align_and_estimate_abundance.pl --transcripts c_puella.fasta --seqType fq --left cm_c_1_S1_R1_001.fastq.gz --right cm_c_1_S1_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir cm_c_1
align_and_estimate_abundance.pl --transcripts c_puella.fasta --seqType fq --left cm_c_2_S9_R1_001.fastq.gz --right cm_c_2_S9_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir cm_c_2
align_and_estimate_abundance.pl --transcripts c_puella.fasta --seqType fq --left cm_t_1_S2_R1_001.fastq.gz --right cm_t_1_S2_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir cm_t_1
align_and_estimate_abundance.pl --transcripts c_puella.fasta --seqType fq --left cm_t_2_S10_R1_001.fastq.gz --right cm_t_2_S10_R2_001.fastq.gz --est_method kallisto --trinity_mode --output_dir cm_t_2
# switch the columns and add a header
awk '{ print $2 " " $1}' c_puella.fasta.gene_trans_map | sed '1 i\TXNAME\tGENEID' | sed 's/ /\t/g' > tx2gene.tsv
```

## Read in data and create sample descriptions

```r
library(DESeq2)
library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(tibble)
# this assumes you are the directory containing the kallisto output
dir <- "./"
run <- c("full_1", "full_2", "wnd_1", "wnd_2", "imm_1", "imm_2")
files <- file.path(dir, run, "abundance.tsv")
names(files) <- run
all(file.exists(files))
tx2gene <- read_tsv("tx2gene.tsv")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
# this derives the sample names from the count dataframe - safer than doing it manually
sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("treatment", "replicate"), sep = "_", remove = FALSE, extra = "drop")
rownames(sampleTable) <- sampleTable$library
# drop variables that won't be fitted
sampleTable <- dplyr::select(sampleTable, -library, -replicate)
```

## Run `DESeq` with different dispersion fits 

```r
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)
dds <- DESeq(dds, parallel = FALSE, fitType = "local")
ddp <- DESeq(dds, parallel = FALSE, fitType = "parametric")
```

### Choose a dispersion fit

In most cases `fitType = "local"` looks best. There is a third option `fitType = "mean"` which is not appropriate here. Compare the two types of fit:

```r
png("plots/disp_ests.png", width = 2 * 480)
par(mfrow = c(1, 2))
plotDispEsts(dds)
title("local")
plotDispEsts(ddp)
title("parametric")
dev.off()
```

![](https://github.com/Perugolate/puella_RNAseq/blob/master/plots/disp_ests.png)

## Regularized log transformation and PCA

Went with `fitType = "local"`.

```r
# use regularized log transformation prior to PCA
rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
# scree plot of PCs
png("plots/scree.png")
plot(pca)
dev.off()
```

![](https://github.com/Perugolate/puella_RNAseq/blob/master/plots/scree.png)

First three PCs.

```r
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
# this is super hacky
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>%
  .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))
png("plots/pca.png", width = 3 * 480)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()
```

![](https://github.com/Perugolate/puella_RNAseq/blob/master/plots/pca.png)


```r
# extract loadings with rownames as a variable
rot <- as.data.frame(pca$rotation) %>% rownames_to_column
png("plots/lod.png", width = 3 * 480)
par(mfrow = c(1, 3))
plot(rot$PC1, xlab = "", xaxt = "n", ylab = "PC1 loading")
abline(h = 0)
title("A. PC1")
plot(rot$PC2, xlab = "", xaxt = "n", ylab = "PC2 loading")
abline(h = 0)
title("B. PC2")
plot(rot$PC3, xlab = "", xaxt = "n", ylab = "PC3 loading")
title("C. PC3")
abline(h = 0)
dev.off()
```

![](https://github.com/Perugolate/puella_RNAseq/blob/master/plots/lod.png)

## Results

Extract results. Since the only variable is a factor with two levels it is straightforward to extract the results. For more complex models (like the one you will fit for *N. castaneus*) there will be additional arguments to the `results` function.

```r
# adjust ylim as needed
plotMA(results(dds, alpha = 0.05, addMLE = TRUE), ylim = c(-8, 8))
# convention filters on > 2-fold change
res <- results(dds, alpha = 0.05, addMLE = TRUE) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
```
