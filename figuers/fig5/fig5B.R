library(dplyr)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(ggplot2)
library(ggpubr)
library(VarCon)

unique_intronPOS <- read.csv("unique_intronPOS.csv") #  # derived from short-intron-evolution/data_processing/human_introns_examined.R
### 5'SS region (-3 ~ +6)
unique_intronPOS$fivePrime <- 
  as.character(getSeq(genome, 
                      unique_intronPOS$seqnames,
                      start = ifelse(unique_intronPOS$strand=="+",
                                     unique_intronPOS$start - 3,
                                     unique_intronPOS$end -5
                      ), 
                      end = ifelse(unique_intronPOS$strand=="+",
                                   unique_intronPOS$start + 5,
                                   unique_intronPOS$end + 3),
                      strand = unique_intronPOS$strand)
  )
## 3'SS region (-20 ~ +3)
unique_intronPOS$threePrime <- 
  as.character(getSeq(genome, 
                      unique_intronPOS$seqnames,
                      start = ifelse(unique_intronPOS$strand=="+",
                                     unique_intronPOS$end - 19,
                                     unique_intronPOS$start -3
                      ), 
                      end = ifelse(unique_intronPOS$strand=="+",
                                   unique_intronPOS$end + 3,
                                   unique_intronPOS$start + 19),
                      strand = unique_intronPOS$strand)
  )

## maxent
unique_intronPOS$maxent_5ss <- calculateMaxEntScanScore(unique_intronPOS$fivePrime, 5)
unique_intronPOS$maxent_3ss <- calculateMaxEntScanScore(unique_intronPOS$threePrime, 3)
final_df <- unique_intronPOS %>%
  mutate(intron_type = ifelse(intron_length>100,
                              "Long",
                              ifelse(
                                GC_content >= 0.45,"GC-rich\nshort","AT-rich\nshort"
                              ))
  ) %>%
  pivot_longer(cols = c("maxent_5ss", "maxent_3ss"), 
               names_to = "Splice_site", 
               values_to = "MaxEntScore") %>%
  mutate(MaxEntScore = as.numeric(MaxEntScore))


#### plot fig5B ####
# for compare
my_comparisons <- list(
  c("GC-rich\nshort", "AT-rich\nshort"),
  c("AT-rich\nshort", "Long"),
  c("GC-rich\nshort", "Long")
)
# for hline
long_median <- final_df %>%
  mutate(Splice_site = factor(ifelse(Splice_site == "maxent_5ss", "5'SS", "3'SS"),
                              levels = c("5'SS", "3'SS")),
         intron_type = factor(intron_type,
                              levels = c("GC-rich\nshort", "AT-rich\nshort", "Long"))
  ) %>%
  filter(intron_type == "Long") %>%
  group_by(Splice_site) %>%
  summarise(median_score = median(MaxEntScore, na.rm = TRUE), .groups = "drop")

# main plot
maxent_intron <- final_df %>%  
  mutate(Splice_site = factor(ifelse(Splice_site=="maxent_5ss","5'SS","3'SS"),
                              levels = c("5'SS","3'SS")),
         intron_type = factor(intron_type,
                              levels = c("GC-rich\nshort",
                                         "AT-rich\nshort",
                                         "Long"))
  ) %>%
  ggplot(aes(x=intron_type, y=MaxEntScore,
             fill=intron_type)) + 
  geom_boxplot(outliers = F, alpha = 0.5)+
  facet_grid( ~ Splice_site,
              scales = "free_y") +
  geom_hline(data     = long_median,
             aes(yintercept = median_score),
             color    = "grey88",
             linewidth = 0.6,
             linetype = "dashed") +
  xlab("Introns")+
  ylab("Splice site strength\n(MaxEntScore)")+
  stat_compare_means(
    comparisons   = my_comparisons,
    method        = "wilcox.test",
    size = 3,
    label.y = c(12,12,13)) +
  theme_linedraw(base_size = 14.5) +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, color = "black", face = "bold")
  )

maxent_intron 

ggsave(filename = "maxent_GCATintrons.pdf", plot = maxent_intron, 
       width = 6, height = 6)