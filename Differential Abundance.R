# differential abundance using DESeq2
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(readxl)
library(apeglm)

# import data and build phyloseq object
otu_mat <-read_excel("/Users/caseymeili/Desktop/28160R/Fastq/KD/kd_phyloseq_nopbs.xls", sheet="OTU")
tax_mat<- read_excel("/Users/caseymeili/Desktop/28160R/Fastq/KD/kd_phyloseq_nopbs.xls", sheet="taxon")
Meta <-read_excel("/Users/caseymeili/Desktop/28160R/Fastq/KD/kd_phyloseq_nopbs.xls", sheet="Samples")
Meta <- Meta %>%
  tibble::column_to_rownames("Sample")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Meta)
samples
Physeq <-phyloseq(OTU, TAX, samples)

# subset to 7 dpi samples
Physeq_7dpi <- subset_samples(Physeq, DietDay %in% c("Control 7 dpi", "KD 7 dpi"))
Physeq_7dpi <- prune_taxa(taxa_sums(Physeq_7dpi) > 0, Physeq_7dpi)

# agglomerate to genus level
Physeq_genus <- tax_glom(Physeq_7dpi, taxrank = "Genus", NArm = FALSE)

# clean genus labels - replace NA or empty genus with "Unclassified_<Family>"
tax_df <- as.data.frame(tax_table(Physeq_genus))

tax_df$Genus <- ifelse(
  is.na(tax_df$Genus) | tax_df$Genus == "",
  paste0("Unclassified_", tax_df$Family),
  tax_df$Genus)

tax_table(Physeq_genus) <- tax_table(as.matrix(tax_df))

# filter low-abundance genera (helps reduce noise and multiple testing burden)
Physeq_genus_filt <- filter_taxa(Physeq_genus, function(x) sum(x) > 10, TRUE)

# set factor levels for DESeq2
sample_data(Physeq_genus_filt)$DietDay <- factor(
  sample_data(Physeq_genus_filt)$DietDay,
  levels = c("Control 7 dpi", "KD 7 dpi"))

# run DESeq2 at genus level
dds <- phyloseq_to_deseq2(Physeq_genus_filt, ~ DietDay)
dds <- DESeq(dds)

# log fold change shrinkage (for visualization)
res <- lfcShrink(
  dds,
  coef = "DietDay_KD.7.dpi_vs_Control.7.dpi",
  type = "apeglm")

# convert results to dataframe and merge taxonomy
# extract taxonomy with genus names
tax_df <- as.data.frame(tax_table(Physeq_genus_filt)) %>%
  rownames_to_column("OTU")

# Convert DESeq2 results
res_df <- as.data.frame(res) %>%
  rownames_to_column("OTU") %>%
  filter(!is.na(padj))

# Merge taxonomy into results
res_df <- left_join(res_df, tax_df, by = "OTU")

# Use genus as label
res_df$Genus <- as.character(res_df$Genus)

# annotate significance and direction
lfc_threshold <- 1

res_df <- res_df %>%
  mutate(
    significant = padj < 0.05,
    strong = padj < 0.05 & abs(log2FoldChange) > lfc_threshold,
    enriched = ifelse(log2FoldChange > 0, "KD", "Control")
  )

# select top genera for plotting
top_bar <- res_df %>%
  filter(strong) %>%
  slice_max(order_by = abs(log2FoldChange), n = 10) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = factor(Genus, levels = Genus))

# plot (LEfSe-style bar plot)
col_kd <- "#DF7557"
col_control <- "#E6539A"

ggplot(top_bar, aes(x = log2FoldChange, y = Genus, fill = enriched)) +
  geom_bar(stat = "identity") +
  
  scale_fill_manual(values = c(
    "KD" = col_kd,
    "Control" = col_control
  )) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  theme_minimal(base_size = 12) +
  labs(
    x = "log2 Fold Change (KD vs Control, 7 dpi)",
    y = NULL,
    fill = NULL,
    title = "Differentially Abundant Genera"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )

# save plot
ggsave("/Users/caseymeili/Downloads/diff-kd-abund.jpg",
       width = 5,
       height = 5,
       dpi = 300)

# export significant taxa
sig_taxa <- res_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::distinct(OTU, .keep_all = TRUE) %>%
  dplyr::arrange(padj)

write.csv(
  sig_taxa,
  file = "/Users/caseymeili/Downloads/DESeq2_significant_genera_KD_vs_Control_7dpi.csv",
  row.names = FALSE
)
