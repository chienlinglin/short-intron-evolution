library(dplyr)
library(Rsamtools)
library(GenomicRanges)
library(msa)
library(ggplot2)
library(tidyr)
library(forcats)

#### high-GC example - GET4 ####
ageScore4_highGC_GET4 <- readRDS("ageScore4_highGC_GET4.rds")

## length and GC 
gc_df <- ageScore4_highGC_GET4$gc_df
pdf("ageScore4_highGC_GET4.pdf",width = 6,height = 5)
ggplot(gc_df, aes(x = fct_rev(species), y = log10_intronLength, fill = gc)) +
  geom_col(width = 0.6) +
  scale_fill_gradient2(
    high = "brown",     
    mid = "#f4f1ee",     
    low = "cadetblue",    
    midpoint = 0.5,
    limits = c(0, 1),
    oob = scales::squish
  )+
  theme_linedraw() + 
  coord_flip() +
  labs(
    x = "",
    y = expression(log[10]*"(intron length)"),
    fill = "Intron GC-content",
    title = expression(paste("Human ", italic(GET4), " intron 7 conserved junction"))
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "inside",
    legend.background = element_blank(),
    legend.position.inside = c(0.85,0.8)
  ) +
  ylim(c(0,4.5))

dev.off()

## alignemtn 3' and 5' regions (MSA)
alignment_3ss <- ageScore4_highGC_GET4$alignment_3ss

reordered_DNA <- DNAStringSet(apply(alignment_3ss[, 1:(ncol(alignment_3ss) - 2)], 1, paste0, collapse = ""))
names(reordered_DNA) <- alignment_3ss$original_name
alignment3ss_reordered <- msa(reordered_DNA, method = "ClustalW", order = "input")
msaPrettyPrint(
  alignment3ss_reordered,
  output = "pdf",
  file = "ageScore4_highGC_GET4_3ss.pdf",
  showNames = "left",
  showLogo = "top",
  logoColors = "standard area",
  askForOverwrite = FALSE
)

alignment_5ss <- ageScore4_highGC_GET4$alignment_5ss

reordered_DNA <- DNAStringSet(apply(alignment_5ss[, 1:(ncol(alignment_5ss) - 2)], 1, paste0, collapse = ""))
names(reordered_DNA) <- alignment_5ss$original_name
alignment5ss_reordered <- msa(reordered_DNA, method = "ClustalW", order = "input")

msaPrettyPrint(
  alignment5ss_reordered,
  output = "pdf",
  file = "ageScore4_highGC_GET4_5ss.pdf",
  showNames = "left",
  showLogo = "top",
  logoColors = "standard area",
  askForOverwrite = FALSE
)


#### low-GC example - UBA3 ####
ageScore4_highAT_UBA3 <- readRDS("ageScore4_highAT_UBA3.rds")

## length and GC 
gc_df <- ageScore4_highAT_UBA3$gc_df
pdf("ageScore4_LowGC_UBA3.pdf",width = 6,height = 5)
ggplot(gc_df, aes(x = fct_rev(species), y = log10_intronLength, fill = gc)) +
  geom_col(width = 0.6) +
  scale_fill_gradient2(
    high = "brown",     
    mid = "#f4f1ee",     
    low = "cadetblue",    
    midpoint = 0.5,
    limits = c(0, 1),
    oob = scales::squish
  )+
  theme_linedraw() + 
  coord_flip() +
  labs(
    x = "",
    y = expression(log[10]*"(intron length)"),
    fill = "Intron GC-content",
    title = expression(paste("Human ", italic(GET4), " intron 7 conserved junction"))
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "inside",
    legend.background = element_blank(),
    legend.position.inside = c(0.85,0.8)
  ) +
  ylim(c(0,4.5))

dev.off()

## alignemtn 3' and 5' regions (MSA)
alignment_3ss <- ageScore4_highAT_UBA3$alignment_3ss

reordered_DNA <- DNAStringSet(apply(alignment_3ss[, 1:(ncol(alignment_3ss) - 2)], 1, paste0, collapse = ""))
names(reordered_DNA) <- alignment_3ss$original_name
alignment3ss_reordered <- msa(reordered_DNA, method = "ClustalW", order = "input")
msaPrettyPrint(
  alignment3ss_reordered,
  output = "pdf",
  file = "ageScore4_LowGC_UBA3_3ss.pdf",
  showNames = "left",
  showLogo = "top",
  logoColors = "standard area",
  askForOverwrite = FALSE
)

alignment_5ss <- ageScore4_highAT_UBA3$alignment_5ss

reordered_DNA <- DNAStringSet(apply(alignment_5ss[, 1:(ncol(alignment_5ss) - 2)], 1, paste0, collapse = ""))
names(reordered_DNA) <- alignment_5ss$original_name
alignment5ss_reordered <- msa(reordered_DNA, method = "ClustalW", order = "input")

msaPrettyPrint(
  alignment5ss_reordered,
  output = "pdf",
  file = "ageScore4_LowGC_UBA3_5ss.pdf",
  showNames = "left",
  showLogo = "top",
  logoColors = "standard area",
  askForOverwrite = FALSE
)




