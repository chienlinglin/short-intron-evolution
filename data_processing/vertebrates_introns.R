library(rtracklayer)
library(dplyr)
library(ggplot2)
setwd("/data2/Data/ang_chu_data/short_intron/TOGA/NcbiRefgene_geneAnnotationGTF/")

## exon/intron start/end width ##
process_gtf_file <- function(gtf_file, valid_chr_file) {
  
  # Read valid chromosome file
  valid_chr <- read.delim(valid_chr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
    dplyr::filter(!is.na(GC.Percent)) %>% 
    dplyr::select(UCSC.style.name, RefSeq.seq.accession)
  
  # Extract chromosome names from valid_chr
  valid_chr_names <- unique(c(valid_chr$UCSC.style.name, valid_chr$RefSeq.seq.accession))
  print("lodaed valid_chr")
  # Import GTF file and filter by valid chromosomes during import
  gtf_data <- import(gtf_file, format = "gtf")
  print("imported gtf_data")
  gtf_df <- as.data.frame(gtf_data) %>%
    # Filter based on valid chromosome names
    dplyr::filter(seqnames %in% valid_chr_names) %>%
    dplyr::select(-c(score, phase, exon_id, source))
  print("gtf_df selelcted and filtered valid chromosomes")
  # Filter only transcript type
  gtf_df_transcript_only <- gtf_df %>%
    filter(type == "transcript") %>%
    dplyr::rename("pre_mRNA_width" = "width") %>%
    dplyr::select(-exon_number) %>%
    distinct()
  print("filtered gtf_df - transcript")
  # Filter only exon type
  gtf_df_exon_only <- gtf_df %>%
    dplyr::filter(type == "exon") %>%
    dplyr::rename("exon_start" = "start",
                  "exon_end" = "end",
                  "exon_width" = "width")
  print("filtered gtf_df - exon")
  # Create intron data
  gtf_df_intron_only <- gtf_df_exon_only %>%
    group_by(transcript_id) %>%
    mutate(exon_number = as.numeric(exon_number)) %>%
    mutate(intron_start = exon_end + 1,
           intron_end = lead(exon_start - 1),
           type = "intron",
           intron_number = ifelse(strand == "+", exon_number, exon_number - 1)) %>%
    ungroup() %>%
    mutate(intron_width = ifelse(is.na(intron_end), NA, intron_end - intron_start + 1)) %>%
    dplyr::select(seqnames, intron_start, intron_end, intron_width, strand, type, gene_id, transcript_id, gene_name, intron_number) %>%
    na.omit()
  print("filtered gtf_df - intron")
  # Summarize mRNA width
  mRNA_width <- gtf_df_exon_only %>%
    group_by(gene_id, transcript_id) %>%
    summarise(mRNA_width = sum(exon_width, na.rm = TRUE))
  
  gtf_df_transcript_only <- gtf_df_transcript_only %>%
    left_join(mRNA_width, by = c("gene_id", "transcript_id"))
  print("mRNA width done")
  
  return(list(transcript_data = gtf_df_transcript_only,
              exon_data = gtf_df_exon_only,
              intron_data = gtf_df_intron_only))
}
##### process gtf files #####
hg38_gtf <- "/data2/Data/ang_chu_data/from_home/short_intron/hg38.refGene.gtf"
valid_chr_hg38 <- "/data2/Data/ang_chu_data/short_intron/TOGA/chr_refgene_trnasfer/valid_chrom/hg38_sequence_report.tsv"
hg38_list <- process_gtf_file(hg38_gtf,valid_chr_hg38)
saveRDS(hg38_list,"/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information/hg38_list.rds")

mm10_gtf <- "mm10.ncbiRefSeq.gtf"
valid_chr_mm10 <- "/data2/Data/ang_chu_data/short_intron/TOGA/chr_refgene_trnasfer/valid_chrom/mm10_sequence_report.tsv"
mm10_list <- process_gtf_file(mm10_gtf,valid_chr_mm10)
saveRDS(mm10_list,"/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information/mm10_list.rds")

galGal6_gtf <- "galGal6_ncbiRefSeq.gtf"
valid_chr_galGal6 <- "/data2/Data/ang_chu_data/short_intron/TOGA/chr_refgene_trnasfer/valid_chrom/galGal6_sequence_report.tsv"
galGal6_list <- process_gtf_file(galGal6_gtf,valid_chr_galGal6)
saveRDS(galGal6_list,"/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information/galGal6_list.rds")

turtle_gtf <- "GCF_011386835.1_ASM1138683v2.ncbiRefSeq.gtf"
valid_chr_turtle <- "/data2/Data/ang_chu_data/short_intron/TOGA/chr_refgene_trnasfer/valid_chrom/GCF_011386835_sequence_report.tsv"
turtle_list <- process_gtf_file(turtle_gtf,valid_chr_turtle)
saveRDS(turtle_list,"/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information/turtle_list.rds")

xenTro10_gtf <- "xenTro10_ncbiRefSeq.gtf"
valid_chr_xenTro10 <- "/data2/Data/ang_chu_data/short_intron/TOGA/chr_refgene_trnasfer/valid_chrom/xenTro10_sequence_report.tsv"
xenTro10_list <- process_gtf_file(xenTro10_gtf,valid_chr_xenTro10)
saveRDS(xenTro10_list,"/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information/xenTro10_list.rds")

danRer11_gtf <- "danRer11_ncbiRefSeq.gtf"
valid_chr_danRer11 <- "/data2/Data/ang_chu_data/short_intron/TOGA/chr_refgene_trnasfer/valid_chrom/danRer11_sequence_report.tsv"
danRer11_list <- process_gtf_file(danRer11_gtf,valid_chr_danRer11)
saveRDS(danRer11_list,"/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information/danRer11_list.rds")
#### plot transcript/intron length ####
transcript_width_df <- data.frame(
  species = factor(rep(c("hg38", "mm10", "galGal6", "turtle", "xenTro10", "danRer1"), 
                       times = c(length(hg38_list$transcript_data$pre_mRNA_width),
                                 length(mm10_list$transcript_data$pre_mRNA_width),
                                 length(galGal6_list$transcript_data$pre_mRNA_width),
                                 length(turtle_list$transcript_data$pre_mRNA_width),
                                 length(xenTro10_list$transcript_data$pre_mRNA_width),
                                 length(danRer1_list$transcript_data$pre_mRNA_width))),
                   levels = c("hg38", "mm10", "galGal6", "turtle", "xenTro10", "danRer1")),
  pre_mRNA_width = c(hg38_list$transcript_data$pre_mRNA_width,
                     mm10_list$transcript_data$pre_mRNA_width,
                     galGal6_list$transcript_data$pre_mRNA_width,
                     turtle_list$transcript_data$pre_mRNA_width,
                     xenTro10_list$transcript_data$pre_mRNA_width,
                     danRer1_list$transcript_data$pre_mRNA_width),
  mRNA_width = c(hg38_list$transcript_data$mRNA_width,
                 mm10_list$transcript_data$mRNA_width,
                 galGal6_list$transcript_data$mRNA_width,
                 turtle_list$transcript_data$mRNA_width,
                 xenTro10_list$transcript_data$mRNA_width,
                 danRer1_list$transcript_data$mRNA_width)
)


ggplot(transcript_width_df, aes(x = species, y = log2(mRNA_width),fill=species)) +
  geom_boxplot() +
  labs(title = "pre_mRNA Width across Species",
       x = "Species",
       y = "log2(mRNA Width)")+ theme(legend.position = "none")

intron_width_df <- data.frame(
  species = factor(rep(c("hg38", "mm10", "galGal6", "turtle", "xenTro10", "danRer1"), 
                       times = c(length(hg38_list$intron_data$intron_width),
                                 length(mm10_list$intron_data$intron_width),
                                 length(galGal6_list$intron_data$intron_width),
                                 length(turtle_list$intron_data$intron_width),
                                 length(xenTro10_list$intron_data$intron_width),
                                 length(danRer1_list$intron_data$intron_width))),
                   levels = c("hg38", "mm10", "galGal6", "turtle", "xenTro10", "danRer1")),
  intron_width = c(hg38_list$intron_data$intron_width,
                   mm10_list$intron_data$intron_width,
                   galGal6_list$intron_data$intron_width,
                   turtle_list$intron_data$intron_width,
                   xenTro10_list$intron_data$intron_width,
                   danRer1_list$intron_data$intron_width)
)

ggplot(intron_width_df, aes(x = species, y = log10(intron_width),fill=species)) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  labs(title = "Intron Width across Species",
       x = "Species",
       y = "log10(intron Width)")+ theme(legend.position = "none")


#### exon/intron GC-content ####
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Ggallus.UCSC.galGal6)
## turtle and xenTro10 sequences were obtained from fasta files. ##
## refer to /data2/Data/ang_chu_data/short_intron/TOGA/fasta_files/
library(BSgenome.Drerio.UCSC.danRer11)
hg38_genome <- BSgenome.Hsapiens.UCSC.hg38
mm10_genome <- BSgenome.Mmusculus.UCSC.mm10
galGal6_genome <- BSgenome.Ggallus.UCSC.galGal6
danRer11_genome <- BSgenome.Drerio.UCSC.danRer11

#### exon ####
exon_hg38 <- hg38_list$exon_data
exon_mm10 <- mm10_list$exon_data
exon_galGal6 <- galGal6_list$exon_data
exon_danRer11 <- danRer11_list$exon_data

exon_hg38$exon_seq <- getSeq(hg38_genome, exon_hg38$seqnames,
                             start=exon_hg38$exon_start, 
                             end=exon_hg38$exon_end,
                             strand=exon_hg38$strand) %>% as.character()

exon_mm10$exon_seq <- getSeq(mm10_genome, exon_mm10$seqnames,
                             start=exon_mm10$exon_start, 
                             end=exon_mm10$exon_end,
                             strand=exon_mm10$strand) %>% as.character()

exon_galGal6$exon_seq <- getSeq(galGal6_genome, exon_galGal6$seqnames,
                                start=exon_galGal6$exon_start, 
                                end=exon_galGal6$exon_end,
                                strand=exon_galGal6$strand) %>% as.character()

exon_danRer11$exon_seq <- getSeq(danRer11_genome, exon_danRer11$seqnames,
                                 start=exon_danRer11$exon_start, 
                                 end=exon_danRer11$exon_end,
                                 strand=exon_danRer11$strand) %>% as.character()

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
exon_hg38$exon_gc <- calculate_gc_content(exon_hg38$exon_seq)
exon_mm10$exon_gc <- calculate_gc_content(exon_mm10$exon_seq)
exon_galGal6$exon_gc <- calculate_gc_content(exon_galGal6$exon_seq)
exon_danRer11$exon_gc <- calculate_gc_content(exon_danRer11$exon_seq)
exon_turtle <- readRDS("/data2/Data/ang_chu_data/short_intron/TOGA/fasta_files/exon_turtle.rds")
exon_xenTro10 <- readRDS("/data2/Data/ang_chu_data/short_intron/TOGA/fasta_files/exon_xenTro10.rds")
setwd("/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information")
all_exon_list <- list(
  exon_hg38 = exon_hg38,
  exon_mm10 = exon_mm10,
  exon_galGal6 = exon_galGal6,
  exon_turtle = exon_turtle,
  exon_xenTro10 = exon_xenTro10,
  exon_danRer11 = exon_danRer11
)
saveRDS(all_exon_list,"all_exon_list.rds")
#### intron ####
intron_hg38 <- hg38_list$intron_data
intron_mm10 <- mm10_list$intron_data
intron_galGal6 <- galGal6_list$intron_data
intron_danRer11 <- danRer11_list$intron_data

intron_hg38$intron_seq <- getSeq(hg38_genome, intron_hg38$seqnames,
                                 start=intron_hg38$intron_start, 
                                 end=intron_hg38$intron_end,
                                 strand=intron_hg38$strand) %>% as.character()

intron_mm10$intron_seq <- getSeq(mm10_genome, intron_mm10$seqnames,
                                 start=intron_mm10$intron_start, 
                                 end=intron_mm10$intron_end,
                                 strand=intron_mm10$strand) %>% as.character()

intron_galGal6$intron_seq <- getSeq(galGal6_genome, intron_galGal6$seqnames,
                                    start=intron_galGal6$intron_start, 
                                    end=intron_galGal6$intron_end,
                                    strand=intron_galGal6$strand) %>% as.character()

intron_danRer11$intron_seq <- getSeq(danRer11_genome, intron_danRer11$seqnames,
                                     start=intron_danRer11$intron_start, 
                                     end=intron_danRer11$intron_end,
                                     strand=intron_danRer11$strand) %>% as.character()

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
intron_hg38$intron_gc <- calculate_gc_content(intron_hg38$intron_seq)
intron_mm10$intron_gc <- calculate_gc_content(intron_mm10$intron_seq)
intron_galGal6$intron_gc <- calculate_gc_content(intron_galGal6$intron_seq)
intron_danRer11$intron_gc <- calculate_gc_content(intron_danRer11$intron_seq)
intron_turtle <- readRDS("/data2/Data/ang_chu_data/short_intron/TOGA/fasta_files/intron_turtle.rds")
intron_xenTro10 <- readRDS("/data2/Data/ang_chu_data/short_intron/TOGA/fasta_files/intron_xenTro10.rds")
setwd("/data2/Data/ang_chu_data/short_intron/TOGA/r_scripts_for_6_species/00_species_information")
all_intron_list <- list(
  intron_hg38 = intron_hg38,
  intron_mm10 = intron_mm10,
  intron_galGal6 = intron_galGal6,
  intron_turtle = intron_turtle,
  intron_xenTro10 = intron_xenTro10,
  intron_danRer11 = intron_danRer11
)
saveRDS(all_intron_list,"all_intron_list.rds")