library(rtracklayer)
library(dplyr)
library(tidyverse)
library(purrr)
library(GenomicRanges)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)

#### build GRCm38 introns df ####
build_introns_from_gtf <- function(gtf_path) {
  
  # ---- Load GTF ----
  message("Loading GTF...")
  gtf <- import(gtf_path)
  
  # ---- Extract exons ----
  message("Extracting exons...")
  exons <- gtf[gtf$type == "exon"]
  exon_df <- as.data.frame(exons)
  
  # ---- Remove version numbers ----
  exon_df$transcript_id <- sub("\\..*$", "", exon_df$transcript_id)
  exon_df$gene_id       <- sub("\\..*$", "", exon_df$gene_id)
  
  # ---- Build introns ----
  message("Building intron coordinates...")
  intron_df <- exon_df %>%
    arrange(seqnames, transcript_id, start) %>%
    group_by(transcript_id) %>%
    filter(n() > 1) %>%                      # must have >1 exon
    mutate(
      exon_number = as.numeric(exon_number),
      intron_start = lag(end) + 1,
      intron_end   = start - 1
    ) %>%
    filter(!is.na(intron_start) & intron_end > intron_start) %>%
    transmute(
      seqnames = as.character(seqnames),
      intron_start,
      intron_end,
      strand = strand,
      gene_id,
      transcript_id,
      exon_number,
      intron_number = ifelse(strand == "+",
                             exon_number - 1,
                             exon_number),
      intron_length = intron_end - intron_start + 1
    ) %>%
    ungroup()
  
  # ---- Build complete transcript list ----
  message("Preparing full transcript table...")
  all_tx <- exon_df %>%
    distinct(transcript_id, gene_id, seqnames, strand)
  
  intron_df_full <- all_tx %>%
    left_join(intron_df,
              by = c("transcript_id","gene_id","seqnames","strand"))
  
  message("Done.")
  
  # Output everything
  return(list(
    exon_df = exon_df,
    intron_df = intron_df,
    intron_df_full = intron_df_full
  ))
}
mouse_gtf_path <- "gencode.vM21.annotation.gtf"
mouse_res <- build_introns_from_gtf(mouse_gtf_path)

mouse_intron <- mouse_res$intron_df

### intron GC ####
intron_gr <- GRanges(
  seqnames = mouse_intron$seqnames,
  ranges   = IRanges(start = mouse_intron$intron_start,
                     end   = mouse_intron$intron_end),
  strand   = mouse_intron$strand
)

intron_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, intron_gr)

gc_content <- letterFrequency(intron_seqs, letters = c("G", "C"),
                              as.prob = TRUE) %>% rowSums()

mouse_intron$intron_GC <- gc_content

mouse_intron %>%
  mutate(intron_type = factor(ifelse(
    intron_length <= 100, "short","long" 
  ),levels = c("short","long"))) %>%
  ggplot(aes(x = intron_GC, fill = intron_type)) +
  geom_density(alpha = .3) + 
  geom_vline(xintercept = 0.45,linetype = 2) +
  theme_linedraw()+
  theme(legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.85),
        strip.text = element_blank()
  ) +
  labs(x = "Intron GC Content", y = "Density", fill = "Intron Type")

#### intron_type ####  
gene_intron_class <-  mouse_GRCm38_res$intron_df_full %>%
  left_join(gene,
            by = c("gene_id")) %>%
  left_join(mouse_intron, by =c("transcript_id","gene_id","seqnames","strand","intron_start","intron_end",
                                "exon_number","intron_number","intron_length")) %>%
  mutate(
    GC_types = case_when(
      is.na(intron_GC)      ~ NA_character_,
      intron_GC >= 0.45     ~ "GC-rich",
      TRUE                  ~ "AT-rich"
    ),
    is_short = intron_length <= 100 & intron_length >= 30
  ) %>%
  group_by(gene_id, gene_name,gene_type) %>%
  summarise(
    has_intron      = any(!is.na(intron_length)),
    any_short       = any(is_short, na.rm = TRUE),
    all_long        = if (any(!is.na(intron_length))) all(intron_length > 100, na.rm = TRUE) else FALSE,
    # only look at GC_types among short introns
    short_GC_types  = list(unique(GC_types[is_short])),
    .groups         = "drop"
  ) %>%
  mutate(
    intron_class = case_when(
      any_short & sapply(short_GC_types, function(x) all(c("GC-rich", "AT-rich") %in% x)) ~ "short_mixed",
      any_short & sapply(short_GC_types, function(x) "GC-rich" %in% x)     ~ "short_GC-rich",
      any_short & sapply(short_GC_types, function(x) "AT-rich" %in% x)     ~ "short_AT-rich",
      TRUE                                                                   ~ "long/no_intron_gene"
    )
  ) %>%
  dplyr::select(-short_GC_types)


saveRDS(gene_intron_class,"GRCm38_mm10_mouse_geneAnnotation.rds")
