library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

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
limma_DET_SRRM2 <- read_csv("rScriptsAndData/limma_DET_SRRM2.csv") %>%
  dplyr::rename("transcriptID" = "...1") %>%
  mutate(significance = case_when(`adj.P.Val`<=0.05 & logFC>=1 ~ "upregulated\n(n=617)",
                                  `adj.P.Val`<=0.05 & logFC<=-1 ~ "downregulated\n(n=160)",
                                  TRUE ~ "non-sig\n(n=11375)"))
limma_DET_SRRM2$significance <- factor(limma_DET_SRRM2$significance,
                                       levels = c("non-sig\n(n=11375)",
                                                  "downregulated\n(n=160)",
                                                  "upregulated\n(n=617)"))
#### extract the transcript information ####
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
txdb <- makeTxDbFromGFF("gencode.v38.annotation.gtf")

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
    value = unlist(len),
    type = "exon_length",
    order = seq_along(len)
  ) %>%
    bind_rows(data.frame(
      transcriptID = tx,
      value = unlist(gc),
      type = "exon_gc",
      order = seq_along(gc)
    ))
}, names(exon_lengths), exon_lengths, exon_gc))
row.names(exon_df) <- 1:nrow(exon_df)
# saveRDS(exon_df,"exon_df.rds")
#### Extract introns ####
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
  if (length(len) == 0 || length(gc) == 0) {
    return(data.frame(
      transcriptID = tx,
      value = NA,  # fill with NA
      type = c("intron_length", "intron_gc"),
      order = NA  # fill with NA
    ))
  }
  
  data.frame(
    transcriptID = tx,
    value = unlist(len),
    type = "intron_length",
    order = seq_along(len)
  ) %>%
    bind_rows(data.frame(
      transcriptID = tx,
      value = unlist(gc),
      type = "intron_gc",
      order = seq_along(gc)
    ))
}, names(intron_lengths), intron_lengths, intron_gc))
row.names(intron_df) <- 1:nrow(intron_df)
# saveRDS(intron_df,"intron_df.rds")

##### combine exon/intron information #####
son_trans <- limma_DET_SON %>% dplyr::select(transcriptID,logFC,adj.P.Val,significance)
srrm2_trans <- limma_DET_SRRM2 %>% dplyr::select(transcriptID,logFC,adj.P.Val,significance)

#### plot ####
plot_exon_intron_analysis <- function(input, comparisons, target_protein) {
  test_ex <- merge(input, exon_df, by = "transcriptID")
  test_in <- merge(input, intron_df, by = "transcriptID")
  
  # Define the function
  plot_boxplot <- function(data, type_filter, y_transform = identity, y_limits) {
    data %>%
      filter(type == type_filter & !is.na(value)) %>%
      group_by(transcriptID, significance) %>% 
      summarise(mean_value = mean(value), .groups = "drop") %>%
      ggplot(aes(x = as.factor(significance), y = y_transform(mean_value), fill = as.factor(significance))) +
      geom_violin(width = 1, alpha = 0.4, color = NA, trim = FALSE) +
      geom_boxplot(width = 0.2) + 
      scale_fill_manual(values=c("gray", "lightskyblue","tomato")) +
      stat_compare_means(comparisons = comparisons, 
                         label.y = if (type_filter %in% c("exon_length", "intron_length")) {
                           c(6, 6, 6.4) 
                         } else {
                           c(0.86, 0.86, 0.91)
                         })+
      labs(y = ifelse(type_filter %in% c("exon_length", "intron_length"), 
                      "log10(Length)", "GC-content"), 
           x = target_protein) +
      ylim(y_limits) + 
      theme_minimal() +
      theme(text = element_text(size = 18),legend.position = "none")  # Remove legend
  }
  
  # Generate plots with fixed y-axis scales
  p1 <- plot_boxplot(test_ex, "exon_length", log10, c(0, 7.1))
  p2 <- plot_boxplot(test_ex, "exon_gc", identity, c(0, 1.0))
  p3 <- plot_boxplot(test_in, "intron_length", log10, c(0, 7.1))
  p4 <- plot_boxplot(test_in, "intron_gc", identity, c(0, 1.0))
  
  # Show all plots
  ggarrange(p1, p2, p3, p4, 
            labels = c("Exon", "Exon", "Intron", "Intron"), 
            ncol = 2, nrow = 2,
            label.x = 0.1)
}


son_comparisons <- list(c("non-sig\n(n=10257)", "downregulated\n(n=505)"), 
                        c("downregulated\n(n=505)", "upregulated\n(n=1390)"),
                        c("non-sig\n(n=10257)", "upregulated\n(n=1390)"))

srrm2_comparisons <- list(c("non-sig\n(n=11375)", "downregulated\n(n=160)"), 
                          c("downregulated\n(n=160)", "upregulated\n(n=617)"),
                          c("non-sig\n(n=11375)", "upregulated\n(n=617)"))


son_plots <- plot_exon_intron_analysis(son_trans, son_comparisons, "SON")

pdf("transLenGC_son.pdf",width = 9.5,height = 9)
print(son_plots)
dev.off()
