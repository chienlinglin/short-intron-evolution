library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggseqlogo)
library(cowplot)

unique_intronPOS <- read.csv("unique_intronPOS.csv") # derived from short-intron-evolution/data_processing/human_introns_examined.R
## the BP strength is predicted by svm-bpfinder (full length of intron sequences)
## bash for SVM 
## . SVM_bp.sh inputSVM_intronSeq.fa Hsap 100 > outSVM_intronBP.fa

## read in BP/PPT prediction result
## choose the maximum svm_scr within the same intron
all_max_scr <- read.csv("all_max_scr.csv")

## plot df
plot_df <- all_max_scr %>%
  separate(col = seq_id, 
           into = c("intron_pos","intron_length","intron_GC","intron_GTAG"),
           sep = "\\|"
  ) %>%
  mutate(group = factor(ifelse(intron_length>100,
                               "Long",
                               ifelse(
                                 intron_GC >= 0.45,"GC-rich\nshort","AT-rich\nshort"
                               )),
                        levels = c("GC-rich\nshort","AT-rich\nshort",
                                   "Long")
                        
  )
  )

# bp_scr (svm_scr)
my_comparisons <- list(
  c("GC-rich\nshort", "AT-rich\nshort"),
  c("GC-rich\nshort", "Long"),
  c("AT-rich\nshort", "Long")
)
bp_scr_plot <- ggplot(plot_df, aes(x = group, y = bp_scr, fill = group)) +
  geom_boxplot(width=0.7 ,alpha = 0.7) +
  labs(x = "", y = "Branch point scores") +
  theme_linedraw(base_size = 14) +
  stat_compare_means(
    comparisons   = my_comparisons,
    method        = "wilcox.test",
    size = 3,
    label.y = c(4.5,5.1,4.5)) + 
  scale_y_continuous(limits = c(NA, 5.7)) + 
  theme(legend.position = "none" ,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank()
  )

bp_scr_plot

# seq
shortGC_seqs <- as.character(plot_df$bp_seq[which(plot_df$group=="GC-rich\nshort")]) %>% toupper()
shortAT_seqs <- as.character(plot_df$bp_seq[which(plot_df$group=="AT-rich\nshort")]) %>% toupper()
long_seqs  <- as.character(plot_df$bp_seq[which(plot_df$group=="Long")]) %>% toupper()
seqs_list <- list(`GC-rich\nshort` = shortGC_seqs, 
                  `AT-rich\nshort` = shortAT_seqs, 
                  `Long\nintrons` = long_seqs)

bp_seq_logo <- ggseqlogo(seqs_list, method = "bits", ncol = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, color = "black", face = "bold")
  )

bp_scr_plot <- bp_scr_plot +
  theme(plot.margin = margin(t = 5, r = 5, b = 0, l = 5))  
bp_seq_logo <- bp_seq_logo +
  theme(plot.margin = margin(t = 0, r = 5, b = 5, l = 5))

#### plot fig5C ####
combine_bp <- plot_grid(bp_scr_plot, bp_seq_logo,  
                        ncol = 1, align = 'v',
                        rel_heights = c(3, 2))

ggsave(combine_bp,filename = "combine_introns_bp.pdf",
       width = 4.5,height = 5.5)