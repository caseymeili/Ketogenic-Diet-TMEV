# differential abundance using DESeq2
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(readxl)
library(apeglm)

# import data and build phyloseq object
otu_mat <- read_excel("/Users/caseymeili/Desktop/28160R/Fastq/KD/kd_phyloseq_nopbs.xls", sheet="OTU")
tax_mat <- read_excel("/Users/caseymeili/Desktop/28160R/Fastq/KD/kd_phyloseq_nopbs.xls", sheet="taxon")
Meta <- read_excel("/Users/caseymeili/Desktop/28160R/Fastq/KD/kd_phyloseq_nopbs.xls", sheet="Samples")
Meta <- Meta %>% tibble::column_to_rownames("Sample")
otu_mat <- otu_mat %>% tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>% tibble::column_to_rownames("#OTU ID")
OTU <- otu_table(as.matrix(otu_mat), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_mat))
samples <- sample_data(Meta)
Physeq <- phyloseq(OTU, TAX, samples)

# subset to 7 dpi samples
Physeq_7dpi <- subset_samples(Physeq, DietDay %in% c("Control 7 dpi", "KD 7 dpi"))
Physeq_7dpi <- prune_taxa(taxa_sums(Physeq_7dpi) > 0, Physeq_7dpi)

# agglomerate to genus level
Physeq_genus <- tax_glom(Physeq_7dpi, taxrank = "Genus", NArm = FALSE)

# clean genus labels
tax_df <- as.data.frame(tax_table(Physeq_genus))

tax_df$Genus <- ifelse(
  is.na(tax_df$Genus) | tax_df$Genus == "",
  paste0("Unclassified_", tax_df$Family),
  tax_df$Genus)

tax_table(Physeq_genus) <- tax_table(as.matrix(tax_df))

# filter low abundance taxa
Physeq_genus_filt <- filter_taxa(Physeq_genus, function(x) sum(x) > 10, TRUE)

# set factor order
sample_data(Physeq_genus_filt)$DietDay <- factor(
  sample_data(Physeq_genus_filt)$DietDay,
  levels = c("Control 7 dpi", "KD 7 dpi"))

# run DESeq2
dds <- phyloseq_to_deseq2(Physeq_genus_filt, ~ DietDay)
dds <- DESeq(dds)

# raw DESeq2 results for export/stats
res <- results(dds, contrast = c("DietDay", "KD 7 dpi", "Control 7 dpi"))

tax_join <- as.data.frame(tax_table(Physeq_genus_filt)) %>%
  rownames_to_column("OTU")

res_df <- as.data.frame(res) %>%
  rownames_to_column("OTU") %>%
  filter(!is.na(padj))

res_df <- left_join(res_df, tax_join, by = "OTU")
res_df$Genus <- as.character(res_df$Genus)

# annotate (raw statistics only)
lfc_threshold <- 1

res_df <- res_df %>%
  mutate(significant = padj < 0.05,
         strong = padj < 0.05 & abs(log2FoldChange) > lfc_threshold,
         enriched = ifelse(log2FoldChange > 0, "KD", "Control"))

# export raw DESeq2
sig_taxa <- res_df %>%
  filter(padj < 0.05) %>%
  distinct(OTU, .keep_all = TRUE) %>%
  arrange(padj)

write.csv(
  sig_taxa,
  file = "/Users/caseymeili/Downloads/KD-DESeq2_significant_genera_KD_vs_Control_7dpi.csv",
  row.names = FALSE)

# shrinkage (for visualization only)
res_shrunk <- lfcShrink(
  dds,
  coef = "DietDay_KD.7.dpi_vs_Control.7.dpi",
  type = "apeglm")

plot_df <- plot_df %>%
  mutate(
    enriched = ifelse(log2FoldChange > 0, "KD", "Control"),
    significant = padj < 0.05)

plot_df$Genus <- as.character(plot_df$Genus)

# select top taxa for plotting (shrunk values only)
top_bar <- plot_df %>%
  filter(padj < 0.05) %>%
  slice_max(order_by = abs(log2FoldChange), n = 10) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = factor(Genus, levels = Genus))

# define figure colors
col_kd <- "#DF7557"
col_control <- "#E6539A"

# plot
ggplot(top_bar, aes(x = log2FoldChange, y = Genus, fill = enriched)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("KD" = col_kd, "Control" = col_control)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 12) +
  labs(x = "log2 Fold Change (KD vs Control, 7 dpi)",
       y = NULL,
       fill = NULL,
       title = "Differentially Abundant Genera") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10))

# export plot
ggsave("/Users/caseymeili/Downloads/diff-kd-abund.jpg",
       width = 5,
       height = 5,
       dpi = 300)
