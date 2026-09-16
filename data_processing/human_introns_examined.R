# ============================================================
# Purpose:
#   Generate a non-redundant set of human introns from the
#   GENCODE v38 intron annotation and export short- and long-
#   intron BED files for downstream splicing-signal analyses.
#
#   The resulting intron set is defined by unique genomic
#   coordinates and strand, without restricting to canonical
#   GT-AG splice sites. Only introns >= 30 nt are retained.
#
#   This intron set is used for downstream analyses including:
#     1. 5' and 3' splice-site sequence logos
#     2. Splice-site strength analysis using MaxEntScan
#     3. Branch-point strength analysis using SVM-BPfinder
#     4. CLIP-seq enrichment analysis
#
# Input:
#   - intron_df_GC_GRCh38.rds
#     Generated from the GENCODE v38 annotation. (getIntrons_GRCh38_gtf.R)
#     Contains transcript-level intron coordinates, intron length,
#     GC content, and GT-AG splice-site annotation.
#
# Output:
#   - unique_intronPOS.rds
#   - unique_intronPOS.csv
#       Non-redundant human introns >= 30 nt, defined by unique
#       genomic coordinates and strand.
#
#   - short_introns.bed
#       BED6 file containing short introns (30-100 nt).
#
#   - long_introns.bed
#       BED6 file containing long introns (>100 nt).
# ============================================================
library(dplyr)
setwd("/data2/Data/ang_chu_data/short_intron/shortIntron_splicing_signals")
### extract from gtf gencode v38
gencode_intron_df <- readRDS("/data2/Data/ang_chu_data/short_intron/shortIntronGene_bedFiles/intron_df_GC_GRCh38.rds") 

unique_intronPOS <- gencode_intron_df %>%
  filter(intron_length >= 30) %>%
  distinct(seqnames, start, end, strand, 
           intron_length, GC_content, GTAG_ss)
saveRDS(unique_intronPOS,file = "unique_intronPOS.rds")
write.csv(unique_intronPOS,"unique_intronPOS.csv",row.names = FALSE)
#
unique_intronPOS

bed_short <- unique_intronPOS %>%
  filter(intron_length <= 100) %>%
  dplyr::mutate(
    chrom = seqnames,
    chromStart = start - 1, # bed file is zero-based
    chromEnd = end,
    name = paste0(seqnames, ":", start, "-", end,"_short"),
    score = 0
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, name, score, strand)

write.table(
  bed_short,
  "short_introns.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

bed_long <- unique_intronPOS %>%
  filter(intron_length > 100) %>%
  dplyr::mutate(
    chrom = seqnames,
    chromStart = start - 1,
    chromEnd = end,
    name = paste0(seqnames, ":", start, "-", end,"_long"),
    score = 0
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, name, score, strand)

write.table(
  bed_long,
  "long_introns.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)