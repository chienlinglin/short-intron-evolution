library(dplyr)
library(readxl)
library(tidyverse)
library(ggplot2)
#### import gene types ####
protein_lncRNA_GRCh38 <- readRDS("protein_lncRNA_GRCh38.rds") # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R
gene_type <- protein_lncRNA_GRCh38 %>%
  distinct(gene_name,intron_class) %>%
  dplyr::rename(Intron_subtype = intron_class)

#### import embryo signature gene set ####
## resource - Capturing totipotency in human cells through spliceosomal repression
# https://ars.els-cdn.com/content/image/1-s2.0-S0092867424005191-mmc2.xlsx
embryo_genes <- read_excel("ZGA-gene-list.xlsx", 
                           sheet = "human embryo geneset", skip = 1)

bulk_info_embryo <- gene_type %>%
  distinct() %>%
  left_join(
    embryo_genes %>% dplyr::select(Gene, gene_set),
    by = c("gene_name" = "Gene")
  ) %>%
  mutate(gene_set = if_else(is.na(gene_set), "non-stage specific", gene_set)) %>%
  filter(gene_set != "non-stage specific")

table(bulk_info_embryo$gene_set)

gene_level_embryo <- bulk_info_embryo %>%
  distinct(gene_name, Intron_subtype, gene_set)


contingency_table <- table(gene_level_embryo$Intron_subtype, gene_level_embryo$gene_set)
contingency_table
# chi-square test: H₀ = Intron type is irrelevant to developmental stage  
# if p < 0.05 → refuse H₀
chisq_result <- chisq.test(contingency_table)
chisq_result$residuals

residual_df <- as.data.frame(as.table(chisq_result$residuals))

colnames(residual_df) <- c("Intron_type", "Stage", "Residual")

#### plot fig4c ####
embryoStage_genes <- residual_df %>%
  mutate(
    Stage = factor(Stage,
                   levels = c("Oocyte", "Zygote-2","4-cell","8-morula",
                              "Blastocyst","ESC","non-stage specific")),
    Intron_type = factor(Intron_type,
                         levels = rev(c("short_GC-rich","short_mixed","short_AT-rich", "long/no_intron_gene"))),
    sig = abs(Residual) > 2 ) %>%
  ggplot(aes(x = Stage,
             y = Intron_type,
             fill = Residual)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Residual)),
            size = 3) +
  scale_fill_gradientn(
    colours = c("steelblue4", "white","white","white", "firebrick"),
    values = scales::rescale(c(-4.5, -1, 0, 1, 4.5)),
    limits = c(-4.5, 4.5)
  ) +
  theme_classic(base_size = 14)+
  labs(x = "Human embryo stage-specific gene set\n(Li et al.)", y = NULL,
       subtitle = paste0("Pearson's χ² test\ndf = ",
                         chisq_result$parameter,
                         "\nχ² =",
                         round(chisq_result$statistic,digits = 1),
                         "\nP-value = 2.4e-12")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(fill = NA, colour = "black"),
    plot.subtitle = element_text(size = 8),
    axis.title = element_text(size = 10)
  ) 

embryoStage_genes

ggsave(filename = "embryoStage_genes.pdf",embryoStage_genes,
       width = 6, height = 3.3)