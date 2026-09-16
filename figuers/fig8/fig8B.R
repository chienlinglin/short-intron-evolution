library(readxl)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
## hg19 was used for TSA seq
## download from doi: 10.1083/jcb.201807108. Epub 2018 Aug 28.
TSA_score_JCB_201807108 <- read_excel("TSA_score_JCB_201807108.xlsx") %>%
  dplyr::select(1:13)


#### with hg38 ####
gene_intronType <- readRDS("protein_lncRNA_GRCh38.rds") %>% # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R
  transmute(gene_id,gene_name,
            intron_class = case_when(
              has_intron == "FALSE" ~ "no\nintron", 
              intron_class == "short_intron_gene" ~ "contain\nshort intron",
              TRUE ~ "only\nlong intron"
            )
  ) 

combine_tsa_score_gene_status <- merge(TSA_score_JCB_201807108,gene_intronType,
                                       by=c("gene_id","gene_name"))

### plot
library(ggplot2)
library(ggpubr)
my_comparisons <- list(
  c("no\nintron", "only\nlong intron"),
  c("no\nintron", "contain\nshort intron"),
  c("only\nlong intron", "contain\nshort intron")
)
sample_counts <- combine_tsa_score_gene_status %>%
  mutate(intron_class = factor(intron_class,
                               levels = c("no\nintron", 
                                          "only\nlong intron", 
                                          "contain\nshort intron"))) %>%
  count(intron_class)
tsa_score_genes_plot <- combine_tsa_score_gene_status %>%
  mutate(intron_class = factor(intron_class,
                               levels = c("no\nintron", 
                                          "only\nlong intron", 
                                          "contain\nshort intron"))) %>%
  ggplot(aes(x = intron_class, y = tsa_value, fill = intron_class)) + 
  
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.4,fill="white") +  
  
  stat_compare_means(
    comparisons = my_comparisons,
    label.y = c(3.2, 3.5, 3.2),  
    method = "wilcox.test"
  ) +
  
  scale_fill_manual(values = c("no\nintron" = "gray20",
                               "only\nlong intron" = "royalblue",
                               "contain\nshort intron" = "firebrick")) +
  
  geom_text(
    data = sample_counts,
    aes(x = intron_class, y = 0, label = paste0("n = ", n)),
    nudge_y = -2.2,  
    size = 4.5,
    inherit.aes = FALSE
  ) +
  
  theme_linedraw(base_size = 15) +
  xlab("Gene types")+
  ylab("TSA score")+
  theme(legend.position = "none",
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.05),  
        panel.grid.minor = element_line(colour = "grey70", linewidth = 0.05),  
        axis.text.x = element_text(
          face = "bold"
        ))

ggsave(tsa_score_genes_plot,
       filename = "tsa_score_genes_plot.pdf",
       width = 4,height = 5,units = "in")

