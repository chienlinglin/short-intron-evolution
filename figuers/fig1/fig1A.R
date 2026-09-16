library(ggplot2)
library(dplyr)


#### new version of short intron distribution ####
# from gencode.v38.annotation.gtf"
gencode_intron_df <- readRDS("intron_df_GC_GRCh38.rds") # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R


#### Fig1A ####
## 01. intronLen_hist_allIntron
intronLen_hist_allIntron <- gencode_intron_df %>% 
  filter(intron_length>=30) %>%
  ggplot(aes(x=log10(intron_length))) + 
  annotate(
    "rect",
    xmin =  -Inf,
    xmax = log10(100),
    ymin = -Inf,
    ymax = Inf,
    fill = "red",
    alpha = 0.15
  ) +
  annotate(
    "rect",
    xmin = log10(100),
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = "steelblue",
    alpha = 0.15
  ) + 
  geom_histogram(fill="#458B74",
                 alpha=0.9,
                 bins=300) +
  annotate(
    "text", x = mean(c(log10(30), log10(100))),y = Inf, label = "Short\nintrons", vjust = 2, colour = "red4", size = 5) +
  annotate(
    "text",x = 2.8,y = Inf, label = "Long\nintrons", vjust = 2, colour = "steelblue4",size = 5 ) +
  theme_linedraw(base_size=14)+
  #geom_vline(xintercept=log10(100),linetype=5,size = 0.6,colour = "gray80") +
  xlab(expression(log[10]*"(Intron length)"))+
  ylab("Count") +
  theme(      
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)
  ) 

intronLen_hist_allIntron

ggsave(filename = "intronLen_hist_allIntron.pdf",
       plot = intronLen_hist_allIntron, width = 7, height = 4)