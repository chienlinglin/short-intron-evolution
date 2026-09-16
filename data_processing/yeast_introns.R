library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3

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
yeast_intron <- build_introns_from_gtf("/data2/Data/ang_chu_data/short_intron/yeast/Saccharomyces_cerevisiae.R64-1-1.62.gtf")

yeast_intron_df <- yeast_intron$intron_df %>%
  mutate(seqnames = paste0("chr",seqnames),
         seqnames = ifelse(seqnames=="chrMito","chrM",seqnames)
  )

yeast_intron_df$intron_seqs <- getSeq(
  genome,
  yeast_intron_df$seqnames,
  start = yeast_intron_df$intron_start,
  end = yeast_intron_df$intron_end,
  strand = yeast_intron_df$strand) %>% 
  as.character()

calculate_gc_content <- function(seq) {
  # Calculate GC content for each sequence independently
  gc_ratios <- sapply(seq, function(s) {
    bases <- unlist(strsplit(as.character(s), ""))
    gc_count <- sum(bases %in% c("G", "C"))
    seq_length <- nchar(as.character(s))
    gc_ratio <- gc_count / seq_length
    return(gc_ratio)
  })
  return(gc_ratios)
}
yeast_intron_df$intron_gc <- calculate_gc_content(yeast_intron_df$intron_seqs)

yeast_intron_df$intron_type <- ifelse(yeast_intron_df$intron_width<=100,"<= 100",">100") # 116
saveRDS(yeast_intron_df,file = "/data2/Data/ang_chu_data/short_intron/yeast/yeast_intron_df.rds")