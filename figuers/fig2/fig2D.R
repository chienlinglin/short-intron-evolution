library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
###### all human proteins/genes in the orthoFinder  #####
ortholog_gene_max_age <- read.csv("ortholog_gene_max_age.csv")

plot_df <- ortholog_gene_max_age %>%
  count(relativeAge, transcriptType) %>%
  dplyr::rename("Gene type"="transcriptType")

plot_df$`Gene type` <- factor(plot_df$`Gene type`,
                              levels = c("noIntron",
                                         "longIntron",
                                         "shortIntronContaining"))
gene_type_colors <- c(
  noIntron = "#A3A3A3",               
  longIntron = "#6B9AC4",             
  shortIntronContaining = "#D46A6A"  
)

gene_stackedBar_ageScore <- ggplot(plot_df, aes(x = factor(relativeAge), y = n, fill = `Gene type`)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Relative Age",
    y = "Number of Genes"
  ) +
  scale_fill_manual(values = gene_type_colors)+
  theme_linedraw(base_size = 14)+ 
  theme(legend.position = c(0.35, 0.85),
        legend.background = element_rect(fill = "transparent", color = NA) 
        )

gene_proportion_ageScore <- ggplot(plot_df, aes(x = factor(relativeAge), y = n, fill = `Gene type`)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    x = "Relative Age",
    y = "Proportion"
  ) +
  scale_fill_manual(values = gene_type_colors)+
  theme_linedraw(base_size = 14)+ 
  theme(legend.position = "none")

combine <- gene_stackedBar_ageScore + gene_proportion_ageScore

combine

ggsave(filename = "gene_ageScore.pdf",plot =combine ,
       device = "pdf",width = 6,height = 6,units = "in",dpi = 300)

