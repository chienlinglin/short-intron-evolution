library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

#### import the differential expression transcripts analysis by limma-voom package ####
limma_DET_SON <- read_csv("limma_DET_SON.csv") %>%
  dplyr::rename("transcriptID" = "...1") %>%
  mutate(significance = case_when(`adj.P.Val`<=0.05 & logFC>=1 ~ "upregulated\n(n=1390)",
                                  `adj.P.Val`<=0.05 & logFC<=-1 ~ "downregulated\n(n=505)",
                                  TRUE ~ "non-sig\n(n=10257)"))
limma_DET_SON$significance <- factor(limma_DET_SON$significance,
                                     levels = c("non-sig\n(n=10257)",
                                                "downregulated\n(n=505)",
                                                "upregulated\n(n=1390)"))

#### extract the transcript information ####
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
txdb <- txdbmaker::makeTxDbFromGFF("gencode.v38.annotation.gtf")
transcript_ids <- limma_DET_SON$transcriptID

#### Function to calculate GC content ####
calc_gc <- function(seq_list) {
  lapply(seq_list, function(dna_set) {
    if (length(dna_set) == 0) return(NA)
    gc_values <- sapply(dna_set, function(exon) { 
      sum(letterFrequency(exon, letters = "GC")) / nchar(exon)
    })
    return(gc_values)
  })
}
#### Extract exons ####
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
exon_info <- exons[names(exons) %in% transcript_ids]
exon_lengths <- lapply(exon_info, width)  # calculate length
exon_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, exon_info)
class(exon_seqs)
exon_gc <- calc_gc(exon_seqs)
exon_df <- do.call(rbind, Map(function(tx, len, gc) {
  data.frame(
    transcriptID = tx,
    exon_length = unlist(len),
    exon_gc = unlist(gc),
    order = seq_along(len)
  )
}, names(exon_lengths), exon_lengths, exon_gc))

row.names(exon_df) <- NULL

exon_counts <- as.data.frame(exons) %>%
  group_by(group_name) %>%
  summarise(n_exons = max(exon_rank), .groups = "drop")

exon_df <- exon_df %>%
  left_join(exon_counts, by = c("transcriptID" = "group_name")) %>%
  group_by(transcriptID) %>%
  mutate(mRNA_length = sum(exon_length),
         median_exon_length = median(exon_length),
         mean_exon_gc = mean(exon_gc)) %>%
  ungroup()

#### intron ####
introns <- intronsByTranscript(txdb, use.names = TRUE)
intron_info <- introns[names(introns) %in% transcript_ids]
intron_counts <- sapply(intron_info, length)
intron_count_df <- data.frame(
  transcriptID = names(intron_counts),
  intron_count = ifelse(intron_counts == 0, NA, intron_counts)
)
intron_lengths <- lapply(intron_info, width)
intron_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, intron_info)
intron_gc <- calc_gc(intron_seqs)
intron_df <- do.call(rbind, Map(function(tx, len, gc) {
  
  # no intron（ single-exon transcript）
  if (length(len) == 0 || length(gc) == 0) {
    return(data.frame(
      transcriptID = tx,
      intron_length = NA,
      intron_gc = NA,
      order = NA
    ))
  }
  
  data.frame(
    transcriptID = tx,
    intron_length = unlist(len),
    intron_gc = unlist(gc),
    order = seq_along(len)
  )
  
}, names(intron_lengths), intron_lengths, intron_gc))

row.names(intron_df) <- NULL
intron_df <- intron_df %>%
  left_join(intron_count_df,by="transcriptID") %>%
  group_by(transcriptID) %>%
  mutate(median_intron_length = median(intron_length),
         mean_intron_gc = mean(intron_gc)) %>%
  ungroup()

#### summary the transcript info ####
enrich_NS <- limma_DET_SON$transcriptID[which(limma_DET_SON$significance=="upregulated\n(n=1390)")]
depleted_NS <- limma_DET_SON$transcriptID[which(limma_DET_SON$significance=="downregulated\n(n=505)")]
summary_df <- inner_join(
  intron_df %>% dplyr::distinct(transcriptID,median_intron_length,mean_intron_gc),
  exon_df %>% dplyr::distinct(transcriptID,mRNA_length,median_exon_length,mean_exon_gc,n_exons),
  by = "transcriptID") %>%
  filter(transcriptID %in% c(enrich_NS,depleted_NS)) %>%
  mutate(NS_enrich = ifelse(transcriptID %in% enrich_NS, 1, 0))

# remonve no introns 
summary_df_w_introns <- summary_df %>%
  filter(!is.na(median_intron_length)) %>%
  mutate(mRNA_length_log = log10(mRNA_length),
         median_intron_length_log = log10(median_intron_length),
         median_exon_length_log = log10(median_exon_length))


#### logistic modeling ####
# Z-score standardization
df_scaled <- summary_df_w_introns %>% 
  mutate(across( c(mean_intron_gc, 
                   mean_exon_gc,
                   n_exons, 
                   mRNA_length_log, 
                   median_intron_length_log, median_exon_length_log),
                 ~ as.numeric(scale(.)) )) 

model <- glm(NS_enrich ~ mean_intron_gc * median_intron_length_log +
               mean_exon_gc + n_exons + mRNA_length_log + 
               median_exon_length_log, 
             data = df_scaled, 
             family = binomial)
summary(model)

library(broom)
tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = dplyr::recode(term,
                         mean_intron_gc = "Intron GC-content",
                         median_intron_length_log =  "Intron length\n(log10)",
                         mean_exon_gc = "Exon GC-content",
                         median_exon_length_log = "Exon length\n(log10)",
                         mRNA_length_log = "mRNA length\n(log10)",
                         n_exons = "Number of exons",
                         `mean_intron_gc:median_intron_length_log` = "Intron GC-content\nx Intron length" 
    ),
    significance_level = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    effect_direction = case_when(
      estimate > 1 ~ "positive",
      estimate < 1 ~ "negative",
      TRUE ~ "neutral"
    ),
    color_group = case_when(
      significance_level == "ns" ~ "Not significant",
      effect_direction == "positive" ~ "Positive effect",
      effect_direction == "negative" ~ "Negative effect"
    )
  )

model_odds_NSgene <- ggplot(tidy_model, aes(x = reorder(term, estimate),
                                            y = estimate,
                                            color = color_group)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, size = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10() +
  scale_color_manual(values = c(
    "Positive effect" = "#D55E00",   
    "Negative effect" = "#0072B2",   
    "Not significant" = "grey70"
  )) +
  geom_text(aes(label = significance_level),
            vjust = -0.5,
            hjust = -0.3,
            size = 4) +
  labs(
    x = "Transcript features",
    y = "Odds Ratio\n(logistic regression)",
    color = ""
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text.y = element_text(size=15))

model_odds_NSgene

ggsave(filename = "model_odds_NSgene.pdf",
       plot = model_odds_NSgene,
       width = 5, height = 4
)