library(stringr)
library(dplyr)
library(tidyr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
library(ggpubr)
library(ggplot2)
library(purrr)
library(patchwork)

#### RI after KD-SON ####
SON_RI.MATS.JC <- read.delim("RI.MATS.JC.txt") 
# Calculate intron coordinates and annotation attributes
SON_RI.MATS.JC <- SON_RI.MATS.JC %>%
  mutate(
    alternativeIS = upstreamEE + 1,  # Convert 0-based to 1-based start
    alternativeIE = downstreamES,   # Already 0-based, compatible with Biostrings
    alt_intronLen = alternativeIE - alternativeIS + 1,  # Compute intron length
    Significance = ifelse(FDR <= 0.05, "sig", "non-sig"),  # Mark significance
    Direction = ifelse(IncLevelDifference < 0, "RI\n(KD SON)", "SI\n(KD SON)")  # Retained or Skipped intron
  )
# Extract alternative intron sequences from the genome
alt_intron_seq <- getSeq(
  genome,
  SON_RI.MATS.JC$chr,
  start = SON_RI.MATS.JC$alternativeIS,
  end = SON_RI.MATS.JC$alternativeIE,
  strand = SON_RI.MATS.JC$strand
)

# Calculate GC content of each intron sequence
alt_intron_gc_content <- (str_count(as.character(alt_intron_seq), "G") + str_count(as.character(alt_intron_seq), "C")) / 
  width(as.character(alt_intron_seq))

# Add sequence and GC content information back to the data frame
SON_RI.MATS.JC <- SON_RI.MATS.JC %>%
  mutate(
    alt_intron_seq = as.character(alt_intron_seq),
    alt_intron_gc_content = alt_intron_gc_content
  )

fill_colors <- c("sig" = "#D87C7C", "non-sig" = "#7C9FD8", "all" = "gray80")
my_comparisons <- list(c("all", "non-sig"), c("non-sig", "sig"), c("all", "sig"))


# Create a combined dataframe including background group ("all")
ri_plot_df <- bind_rows(
  SON_RI.MATS.JC,
  SON_RI.MATS.JC %>% mutate(Significance = "all")
) %>%
  mutate(log10_intronLen = log10(alt_intronLen)) %>%
  dplyr::select(Significance, Direction, log10_intronLen, alt_intron_gc_content)

# Reshape into long format for facet plotting
ri_long_df <- ri_plot_df %>%
  pivot_longer(
    cols = c(log10_intronLen, alt_intron_gc_content),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = recode(Metric,
                         log10_intronLen = "Alternative Intron Length (log10)",
                         alt_intron_gc_content = "Alternative Intron GC-content"))

# Generate the combined violin + boxplot with facet by metric and direction
ri_boxplot_p1 <- ri_long_df %>% filter(Metric=="Alternative Intron Length (log10)") %>%
  ggplot(aes(x = Significance, y = Value, fill = Significance)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width=0.3,outlier.size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = "log10(Length)") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4.2, 4.2, 4.5),  # Adjust if needed
                     method = "wilcox.test") +
  facet_grid(Direction~Metric, scales = "free_y")

ri_boxplot_p2 <- ri_long_df %>% filter(Metric=="Alternative Intron GC-content") %>%
  ggplot(aes(x = Significance, y = Value, fill = Significance)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width=0.3,outlier.size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = "GC-content") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.9, 0.9, 1.0),  # Adjust if needed
                     method = "wilcox.test") +
  facet_grid(Direction~Metric, scales = "free_y")

ri_boxplot <- ri_boxplot_p1/ri_boxplot_p2


#### SE after KD-SON ####
SON_SE.MATS.JC <- read.delim("SE.MATS.JC.txt")
# Calculate intron coordinates and annotation attributes
SON_SE.MATS.JC <- SON_SE.MATS.JC %>%
  mutate(
    upstreamIS = upstreamEE+1,
    upstreamIE = exonStart_0base, # Already 0-based, compatible with Biostrings
    downstreamIS = exonEnd+1,
    downstreamIE = downstreamES, # Already 0-based, compatible with Biostrings
    altExonLen = exonEnd-exonStart_0base,  
    upstreamIL = upstreamIE - upstreamIS + 1,   
    downstreamIL = downstreamIE - downstreamIS + 1,  
    Significance = ifelse(FDR <= 0.05, "sig", "non-sig"),  # Mark significance
    Direction = ifelse(IncLevelDifference < 0, "RE\n(KD SON)", "SE\n(KD SON)")  # Retained or Skipped intron
  )
# Extract alternative intron sequences from the genome
alt_exon_seq <- getSeq(
  genome,
  SON_SE.MATS.JC$chr,
  start = SON_SE.MATS.JC$exonStart_0base+1,
  end = SON_SE.MATS.JC$exonEnd,
  strand = SON_SE.MATS.JC$strand
)
alt_up_intron_seq <- getSeq(
  genome,
  SON_SE.MATS.JC$chr,
  start = SON_SE.MATS.JC$upstreamIS,
  end = SON_SE.MATS.JC$upstreamIE,
  strand = SON_SE.MATS.JC$strand
)
alt_down_intron_seq <- getSeq(
  genome,
  SON_SE.MATS.JC$chr,
  start = SON_SE.MATS.JC$downstreamIS,
  end = SON_SE.MATS.JC$downstreamIE,
  strand = SON_SE.MATS.JC$strand
)
# Calculate GC content of each intron sequence
alt_exon_gc_content <- (str_count(as.character(alt_exon_seq), "G") + 
                          str_count(as.character(alt_exon_seq), "C")) / 
  width(as.character(alt_exon_seq))

alt_up_intron_gc_content <- (str_count(as.character(alt_up_intron_seq), "G") + 
                               str_count(as.character(alt_up_intron_seq), "C")) / 
  width(as.character(alt_up_intron_seq))

alt_down_intron_gc_content <- (str_count(as.character(alt_down_intron_seq), "G") + 
                                 str_count(as.character(alt_down_intron_seq), "C")) / 
  width(as.character(alt_down_intron_seq))

SON_SE.MATS.JC <- SON_SE.MATS.JC %>%
  mutate(
    alt_exon_seq = as.character(alt_exon_seq),
    alt_up_intron_seq = as.character(alt_up_intron_seq),
    alt_down_intron_seq = as.character(alt_down_intron_seq),
    alt_exon_gc_content = alt_exon_gc_content,
    alt_up_intron_gc_content = alt_up_intron_gc_content,
    alt_down_intron_gc_content = alt_down_intron_gc_content
  )
# Create a combined dataframe including background group ("all")
se_plot_df <- bind_rows(
  SON_SE.MATS.JC,
  SON_SE.MATS.JC %>% mutate(Significance = "all")
) %>%
  mutate(log10_altUP_intronLen = log10(nchar(alt_up_intron_seq)),
         log10_altDOWN_intronLen = log10(nchar(alt_down_intron_seq)),
         log10_alt_exonLen = log10(nchar(alt_exon_seq)),
  ) %>%
  dplyr::select(Significance, Direction, 
         log10_altUP_intronLen, log10_altDOWN_intronLen, log10_alt_exonLen,
         alt_up_intron_gc_content,alt_down_intron_gc_content, alt_exon_gc_content)
# Reshape into long format for facet plotting
se_long_df <- se_plot_df %>%
  pivot_longer(
    cols = c(log10_altUP_intronLen, log10_altDOWN_intronLen, log10_alt_exonLen,
             alt_up_intron_gc_content,alt_down_intron_gc_content, alt_exon_gc_content),
    names_to = "Metric",
    values_to = "Value"
  ) 

se_long_df$Metric <- factor(se_long_df$Metric,
                            levels = c("log10_altUP_intronLen","log10_alt_exonLen","log10_altDOWN_intronLen",
                                       "alt_up_intron_gc_content","alt_exon_gc_content","alt_down_intron_gc_content"))
# Generate the combined boxplot with facet by metric and direction
se_boxplot_p1 <- se_long_df %>% filter(Metric=="log10_altUP_intronLen"|
                                         Metric=="log10_alt_exonLen"|
                                         Metric=="log10_altDOWN_intronLen") %>%
  filter(Direction=="SE\n(KD SON)") %>%
  ggplot(aes(x = Significance, y = Value, fill = Significance)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width=0.3,outlier.size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = "log10(Length)") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(5.8, 5.8, 6.1),  # Adjust if needed
                     method = "wilcox.test") +
  facet_grid(Direction~Metric, scales = "free_y")
se_boxplot_p2 <- se_long_df %>% filter(Metric=="alt_up_intron_gc_content"|
                                         Metric=="alt_exon_gc_content"|
                                         Metric=="alt_down_intron_gc_content") %>%
  filter(Direction=="SE\n(KD SON)") %>%
  ggplot(aes(x = Significance, y = Value, fill = Significance)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width=0.3,outlier.size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = "GC-content") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(1, 1, 1.05),  # Adjust if needed
                     method = "wilcox.test") +
  facet_grid(Direction~Metric, scales = "free_y")

se_boxplot <- se_boxplot_p1/se_boxplot_p2


### combine ####
combine <- (ri_boxplot | se_boxplot) +
  plot_layout(widths = c(0.7,2))


ggsave(filename = "combine.pdf",plot = combine,
       device = "pdf",width = 19,height = 9)
