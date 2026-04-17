library(readxl)
library(phyloseq)
library(vegan)
library(tidyverse)
library(scales)
library(grid)

# import phyloseq object
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

# ensure factors are clean
sample_data(Physeq)$DietDay <- factor(sample_data(Physeq)$DietDay)

# bray curtis
pcoa_bc <- ordinate(Physeq, method = "PCoA", distance = "bray")

plot_ordination(Physeq, pcoa_bc, color = "DietDay") +
  geom_point(size = 3) +
  stat_ellipse(aes(color = DietDay), type = "t", linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Control Day -1" = "#085361",
      "KD Day -1" = "#21BFC9",
      "Control 7 dpi" = "#E6539A",
      "KD 7 dpi" = "#DF7557"
    ),
    name = NULL
  ) +
  theme_minimal()

# permanova
bc_dist <- phyloseq::distance(Physeq, method = "bray")

adonis2(bc_dist ~ Diet,
        data = as(sample_data(Physeq), "data.frame"))

# homogeneity of dispersion
disp <- betadisper(bc_dist, sample_data(Physeq)$Diet)
anova(disp)
permutest(disp)

# save plot
ggsave("/Users/caseymeili/Downloads/KD-all-infected-new.jpg",
       height = 5, width = 5, dpi = 300)




# make stacked bar chart
# replace NA / empty genus
tax_mat <- as(tax_table(Physeq), "matrix")
tax_mat[, "Genus"] <- ifelse(is.na(tax_mat[, "Genus"]) | tax_mat[, "Genus"] == "",
                             "Unknown",
                             tax_mat[, "Genus"])
tax_table(Physeq) <- tax_table(tax_mat)

# melt phyloseq object
df_long <- psmelt(Physeq)
df_long$Abundance <- as.numeric(df_long$Abundance)

# define ordered x axis
df_long$DietDay <- factor(
  df_long$DietDay,
  levels = c("Control Day -1",
             "KD Day -1",
             "Control 7 dpi",
             "KD 7 dpi"))

# identify top 10 genera
top_taxa <- df_long %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# group genera outside of top 10 as "Other"
df_long$Genus_plot <- ifelse(df_long$Genus %in% top_taxa,
                             df_long$Genus,
                             "Other")

# reorder genus factor
df_long$Genus_plot <- factor(
  df_long$Genus_plot,
  levels = c(setdiff(unique(df_long$Genus_plot), "Other"), "Other")
)

# normalize within sample
df_long <- df_long %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

# keep sample order
df_long <- df_long %>%
  group_by(DietDay) %>%
  mutate(Sample_order = factor(Sample, levels = unique(Sample))) %>%
  ungroup()

# color palette
custom_colors <- c("#DF7557", "#085361", 
                   "#21BFC9", "#E6539A", 
                   "#011993", "#FFD54F", 
                   "#A02B93", "#1B9E77", 
                   "#554345", "#985DC7",
                    "Other" = "#A3A3A3")

names(custom_colors) <- levels(df_long$Genus_plot)

# plot
ggplot(df_long, aes(x = Sample_order, y = RelAbundance, fill = Genus_plot)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~ DietDay, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing.x = unit(0.6, "lines"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = "Relative Abundance (%)",
    fill = "Genus"
  )

# save
ggsave("/Users/caseymeili/Downloads/four_group_taxa_barplot.jpg",
       width = 8,
       height = 4,
       dpi = 300)



