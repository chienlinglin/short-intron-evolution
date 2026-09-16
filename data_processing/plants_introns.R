library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(purrr)
library(Biostrings)
library(GenomicRanges)
#### 01_extract_introns ####
setwd("/data2/Data/ang_chu_data/plants")
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
Arabidopsis_intron <- build_introns_from_gtf("/data2/Data/ang_chu_data/plants/gtf/Arabidopsis_thaliana.TAIR10.62.gtf")
Oryza_intron <- build_introns_from_gtf("/data2/Data/ang_chu_data/plants/gtf/Oryza_sativa.IRGSP-1.0.62.gtf")
Physcomitrium_intron <- build_introns_from_gtf("/data2/Data/ang_chu_data/plants/gtf/Physcomitrium_patens.Phypa_V3.62.gtf")


combined_df <- bind_rows(
  Arabidopsis_intron$intron_df |> mutate(species = "Arabidopsis thaliana"),
  Oryza_intron$intron_df      |> mutate(species = "Oryza sativa"),
  Physcomitrium_intron$intron_df |> mutate(species = "Physcomitrium patens")
) %>%
  mutate(log_intron_width = log10(intron_length))
saveRDS(object = combined_df,file = "outputs/combined_df.rds")

#### 02_extract_sequece_GC ####
combined_df <- readRDS("/data2/Data/ang_chu_data/plants/outputs/combined_df.rds")
fasta_dir <- "/data2/Data/ang_chu_data/plants/dna_sequence/"

genomes <- list(
  Arabidopsis_thaliana = readDNAStringSet(paste0(fasta_dir, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")),
  Oryza_sativa         = readDNAStringSet(paste0(fasta_dir, "Oryza_sativa.IRGSP-1.0.dna.toplevel.fa")),
  Physcomitrium_patens = readDNAStringSet(paste0(fasta_dir, "Physcomitrium_patens.Phypa_V3.dna.toplevel.fa"))
)
lapply(genomes, function(g) head(names(g), 3))

get_gc <- function(df, genome) {
  names(genome) <- sub(" .*", "", names(genome))
  
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = df$intron_start, end = df$intron_end),
    strand   = df$strand
  )
  
  seqs <- getSeq(genome, gr)
  df$intron_seqs <- as.character(seqs)
  df$intron_GC <- as.numeric(letterFrequency(seqs, letters = "GC", as.prob = TRUE))
  df
}

combined_gc <- bind_rows(
  get_gc(combined_df %>% filter(species == "Arabidopsis thaliana"), genomes$Arabidopsis_thaliana),
  get_gc(combined_df %>% filter(species == "Oryza sativa"),         genomes$Oryza_sativa),
  get_gc(combined_df %>% filter(species == "Physcomitrium patens"), genomes$Physcomitrium_patens)
)

saveRDS(object = combined_gc,file = "outputs/combined_gc.rds")
