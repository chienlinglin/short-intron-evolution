library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
# RNA expression download from GTEx v10
# Referenve genome version: GENCODE Version 39
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
sample.df <- read.delim("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", as.is=TRUE, header=TRUE)
sampleID_tissue <- sample.df %>% dplyr::select(SAMPID,SMTS) %>% 
  mutate(SAMPID = str_replace_all(string = SAMPID ,pattern = "-",replacement = "."))

subject.df <- read.delim("GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt", as.is=TRUE, header=TRUE)
sampid_cols <- colnames(subject.df)[3:ncol(subject.df)]
valid_sampid <- intersect(sampid_cols, sampleID_tissue$SAMPID)

sampid_to_tissue <- setNames(sampleID_tissue$SMTS, sampleID_tissue$SAMPID)
sampid_by_tissue <- split(valid_sampid, sampid_to_tissue[valid_sampid])

subjectMeanSD <- lapply(seq_along(sampid_by_tissue), function(i) {
  sampids <- sampid_by_tissue[[i]]
  tissue_name <- names(sampid_by_tissue)[i]
  
  df <- subject.df[, c("transcript_id", "gene_id", sampids), drop = FALSE]
  df$rowMean <- rowMeans(df[, -c(1, 2)], na.rm = TRUE)
  df$rowSD <- apply(df[, -c(1, 2)], 1, function(x) sd(x, na.rm = TRUE))
  
  df <- df %>% 
    select(transcript_id, gene_id, rowMean, rowSD) %>% 
    mutate(tissue = tissue_name)
  
  return(df)
})

names(subjectMeanSD) <- names(sampid_by_tissue)

#### read different intron types of transcripts ####
read_and_process_gtf <- function(file_path) {
  data <- read.delim(file_path, header = FALSE, comment.char = "#", sep = "\t") 
  # Process the data
  data_1 <- separate(data, V9, into = paste0("V9_", 1:7), sep = "; ", extra = "drop")
  colnames(data_1) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "transcript_id", "gene_type", "gene_name", "transcript_type", "transcript_name", "exon_number")
  data_2 <- data_1 %>%
    mutate(
      gene_id = str_replace(gene_id, pattern = "gene_id ", replacement = ""),
      transcript_id = str_replace(transcript_id, pattern = "transcript_id ", replacement = ""),
      transcript_name = str_replace(transcript_name, pattern = "transcript_name ", replacement = ""),
      exon_number = str_replace(exon_number, pattern = "exon_number ", replacement = ""),
      gene_type = str_replace(gene_type, pattern = "gene_type ", replacement = ""),
      gene_name = str_replace(gene_name, pattern = "gene_name ", replacement = ""),
      transcript_type = str_replace(transcript_type, pattern = "transcript_type ", replacement = "")
    )
  return(data_2)
}
human_ensembl_gtf_df <- read_and_process_gtf("gencode.v39.annotation.gtf")

calculate_introns <- function(df) {
  # Group by transcript_id and sort by start position
  df <- df %>% mutate(exon_number=as.numeric(exon_number))
  exon_df <-  df %>% filter(feature == "exon") %>% arrange(exon_number) %>% group_by(transcript_id) %>%
    mutate(
      intron_start =  ifelse(strand=="+",end + 1,lead(end)+1),
      intron_end = ifelse(strand=="+",lead(start) - 1,start-1),
      intron_width = intron_end - intron_start+1
    )
  return(exon_df)
}
human_ensembl_gtf_w_intron <- calculate_introns(human_ensembl_gtf_df)

transcriptIntronType <- transcriptIntron_info %>% 
  group_by(gene_id,transcript_type, transcript_id) %>%
  summarise(
    intron_status = ifelse(
      any(intronType == "shortIntron"),
      "transcripts with short introns",
      "transcripts without short introns")) %>%
  ungroup()


transcriptIntronType <- transcriptIntron_info %>% 
  group_by(gene_id,transcript_type, transcript_id) %>%
  summarise(
    intron_status = ifelse(
      any(intronType == "shortIntron"),
      "transcripts with short introns",
      "transcripts without short introns")) %>%
  ungroup()

transcriptIntronType # 244939
selectedTranscriptType <- c("Mt_rRNA","Mt_tRNA","miRNA","misc_RNA","rRNA","scRNA","snRNA","snoRNA","ribozyme","sRNA","scaRNA", # Non-coding RNA predicted using sequences from Rfam and miRBase
                            "lncRNA", # Generic long non-coding RNA biotype that replaced the following biotypes: 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic and sense_overlapping.
                            "retained_intron", # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants.
                            "protein_coding" # Contains an open reading frame (ORF).
)
selected_transcriptIntronType <- transcriptIntronType[which(transcriptIntronType$transcript_type %in% selectedTranscriptType),]

#### combine the all tissue data ####
allTissue <- do.call(rbind,subjectMeanSD)
View(allTissue)
allTissue_transcriptType <- left_join(selected_transcriptIntronType,allTissue,by=c("gene_id","transcript_id"))
View(allTissue_transcriptType)
allTissue_transcriptType <- allTissue_transcriptType %>%
  mutate(log_rowMean = log2(rowMean+1))

#### only protein coding transcripts ####
protein_coding_transcripts <- human_ensembl_gtf_w_intron %>% filter(gene_type=="protein_coding") # 164578

transcriptIntronType <- protein_coding_transcripts %>% 
  mutate(intronID = paste0(seqname,":",intron_start,"-",intron_end),
         intronType = case_when(intron_width <= 100 ~ "shortIntron",
                                intron_width > 100 ~ "longIntron",
                                TRUE ~ "withoutIntron"))%>%
  dplyr::select(gene_id,transcript_id,intronID,intron_width,intronType) %>%
  group_by(gene_id,transcript_id) %>%
  summarise(
    intron_status = ifelse(
      any(intronType == "shortIntron"),
      "transcripts with short introns",
      "transcripts without short introns")) %>%
  ungroup()
allTissue <- do.call(rbind,subjectMeanSD)

allTissue_transcriptType <- left_join(transcriptIntronType,allTissue,by=c("gene_id","transcript_id"))
saveRDS(allTissue_transcriptType,"allTissue_transcriptType_proteinCoding.rds")

#### selected_tissues ####
selected_tissues <- c("Muscle","Brain","Nerve","Heart","Lung","Testis","Liver","Kidney","Bladder","Salivary Gland","Adrenal Gland")
allTissue_transcriptType <- allTissue_transcriptType %>%
  mutate(log_rowMean_10 = log10(rowMean+1))

selected_tissues_plot <- allTissue_transcriptType %>% filter(tissue %in% selected_tissues) %>%
  ggplot(aes(x=tissue, y=log_rowMean_10,fill=intron_status)) + 
  geom_boxplot(outliers = FALSE) +
  ylab(expression(log[10]*"(TPM+1)"))+
  xlab("")+ 
  stat_compare_means(aes(group = intron_status), label = "p.signif", label.y = 1.5) +
  scale_fill_discrete(name = "Protein-coding transcripts", 
                      labels = c("with\nshort introns\n(n=30305)", "without\nshort introns\n(n=134273)"))+
  theme_linedraw()+
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        text = element_text(size = 14),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position="bottom")

selected_tissues_plot

ggsave(filename = "selected_tissues_plot.pdf",plot = selected_tissues_plot,
       width = 7,height = 6,units = "in")