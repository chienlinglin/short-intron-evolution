# ============================================================
# Purpose:
#   Derive transcript-level intron coordinates from the GENCODE
#   GRCh38 annotation, calculate intron GC content and canonical
#   GT-AG splice-site status, and generate gene-level intron
#   classifications.
#
# Input:
#   - GENCODE v38 GTF annotation
#   - BSgenome.Hsapiens.UCSC.hg38
#
# Output:
#   - exon_df_GRCh38.rds
#   - intron_df_full_GRCh38.rds
#   - intron_df_GC_GRCh38.rds / .csv
#   - all_gene_GRCh38.rds / .csv
#   - protein_lncRNA_GRCh38.rds
#
# Definitions:
#   - Short intron: 30-100 nt
#   - GC-rich: GC content >= 0.45
#   - AT-rich: GC content < 0.45
#   - Gene classes:
#       short_GC-rich
#       short_AT-rich
#       short_mixed
#       long/no_intron_gene
# ============================================================
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
genome <- BSgenome.Hsapiens.UCSC.hg38
GRCh38_gtf_path <- "gencode.v38.annotation.gtf"

gtf <- import(GRCh38_gtf_path)

exons <- gtf[gtf$type == "exon"]

exon_df <- as.data.frame(exons)

exon_df$transcript_id <- sapply(strsplit(exon_df$transcript_id, "\\."), `[`, 1)  
exon_df$gene_id <- sapply(strsplit(exon_df$gene_id, "\\."), `[`, 1)
saveRDS(exon_df,"exon_df_GRCh38.rds")


intron_df <- exon_df %>%
  arrange(seqnames, transcript_id, start) %>%
  group_by(transcript_id) %>%
  filter(n() > 1) %>%
  mutate(
    intron_start = lag(end) + 1,
    intron_end = start - 1
  ) %>%
  filter(!is.na(intron_start) & intron_end > intron_start) %>%
  transmute(
    seqnames = as.character(seqnames),
    start = intron_start,
    end = intron_end,
    strand = strand,
    gene_id = gene_id,
    transcript_id = transcript_id,
    intron_length = intron_end - intron_start + 1
  )

all_transcripts <- exon_df %>%
  distinct(transcript_id, gene_id, seqnames, strand)

intron_df_full <- all_transcripts %>%
  left_join(intron_df, by = c("transcript_id", "gene_id", "seqnames", "strand"))


saveRDS(intron_df_full,"intron_df_full_GRCh38.rds")

#### get GC-content of introns ####
intron_df_only <- intron_df_full %>%
  filter(!is.na(end))

intron_gr <- GRanges(
  seqnames = intron_df_only$seqnames,
  ranges   = IRanges(start = intron_df_only$start,
                     end   = intron_df_only$end),
  strand   = intron_df_only$strand
)

intron_seqs <- getSeq(genome, intron_gr)

gc_content <- letterFrequency(intron_seqs, letters = c("G", "C"),
                              as.prob = TRUE) %>% rowSums()

intron_df_gc <- intron_df_only %>%
  mutate(
    intron_seq = as.character(intron_seqs),
    GC_content = gc_content,
    GTAG_ss = ifelse(
      substr(intron_seq, 1, 2) == "GT" &
        substr(intron_seq, nchar(intron_seq) - 1, nchar(intron_seq)) == "AG",
      "yes",
      "no"
    )
  )

intron_df_gc <- intron_df_gc %>%
  mutate(GTAG_ss = ifelse(
    substr(intron_seq, 1, 2) == "GT" &
      substr(intron_seq, nchar(intron_seq) - 1, nchar(intron_seq)) == "AG",
    "yes",
    "no"
  ))

intron_df_gc_noSeq <- intron_df_gc %>%
  dplyr::select(-intron_seq)

saveRDS(intron_df_gc_noSeq,"intron_df_GC_GRCh38.rds")

write.csv(
  intron_df_gc_noSeq,
  file = "intron_df_GC_GRCh38.csv",
  row.names = FALSE
)
#### gene-level annotations ####
intron_df_GC_GRCh38 <- readRDS("intron_df_GC_GRCh38.rds") 
intron_df_full <- readRDS("intron_df_full_GRCh38.rds")

gene_intron_class <- intron_df_full %>%
  left_join(intron_df_GC_GRCh38,
            by = c("transcript_id", "gene_id", "seqnames", "strand", "start", "end", "intron_length")) %>%
  mutate(
    GC_types = case_when(
      is.na(GC_content)      ~ NA_character_,
      GC_content >= 0.45     ~ "GC-rich",
      TRUE                  ~ "AT-rich"
    ),
    is_short = intron_length <= 100 & intron_length >= 30
  ) %>%
  group_by(gene_id) %>%
  summarise(
    has_intron      = any(!is.na(intron_length)),
    any_short       = any(is_short, na.rm = TRUE),
    has_long        = if (any(!is.na(intron_length))) all(intron_length > 100, na.rm = TRUE) else FALSE,
    # only look at GC_types among short introns
    short_GC_types  = list(unique(GC_types[is_short])),
    .groups         = "drop"
  ) %>%
  mutate(
    intron_class = case_when(
      any_short & sapply(short_GC_types, function(x) length(x) > 1)        ~ "short_mixed",
      any_short & sapply(short_GC_types, function(x) "GC-rich" %in% x)     ~ "short_GC-rich",
      any_short & sapply(short_GC_types, function(x) "AT-rich" %in% x)     ~ "short_AT-rich",
      TRUE                                                                   ~ "long/no_intron_gene"
    )
  ) %>%
  select(-short_GC_types)

intron_df_annotated <- intron_df_full %>%
  left_join(gene_intron_class, by = "gene_id")


gene_coordinates <- exon_df %>%
  group_by(gene_id,transcript_type,gene_name) %>%
  summarise(
    seqnames = unique(as.character(seqnames))[1], 
    start = min(start), # minimum of the start
    end = max(end), # maximum of the end
    strand = unique(as.character(strand))[1]
  )

gene_info_final <- gene_intron_class %>%
  left_join(gene_coordinates, by = "gene_id")

allGene <- gene_info_final %>% 
  mutate(name = paste(gene_id, gene_name, transcript_type, intron_class,sep = "|"),
         transcript_type = ifelse(transcript_type == "protein_coding"|
                                    transcript_type == "lncRNA",
                                  transcript_type, "other_non_coding"
         )) %>%
  distinct()

saveRDS(allGene,"all_gene_GRCh38.rds")

write.csv(
  allGene,
  file = "all_gene_GRCh38.csv",
  row.names = FALSE
)
protein_lncRNA <- gene_info_final %>% 
  filter(transcript_type=="protein_coding"|transcript_type=="lncRNA") %>%
  mutate(name = paste(gene_id, gene_name, transcript_type, intron_class,sep = "|"))

saveRDS(protein_lncRNA,"/data2/Data/ang_chu_data/short_intron/shortIntronGene_bedFiles/protein_lncRNA_GRCh38.rds")