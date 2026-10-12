library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
# human embryo download from
# https://petropoulos-lanner-labs.clintec.ki.se/resource/human_embryo_stem_datasets/psd.R3.6.em.seurat.ob.rds
# doi:10.1038/s41592-024-02493-2
#### 01. import and process human embryo ####
psd.R3.6.em.seurat.ob <- readRDS("psd.R3.6.em.seurat.ob.rds")
psd.R3.6.em.seurat.ob <- UpdateSeuratObject(psd.R3.6.em.seurat.ob)

####  Extract RAW counts
#### keep this object untouched (never rename its colnames) 
#### so it always stays a reliable lookup table by original sample ID.

counts_mat <- GetAssayData(psd.R3.6.em.seurat.ob, assay = "RNA", layer = "counts")
head(colnames(counts_mat))


#### Pull per-cell/sample metadata
sample_info <- psd.R3.6.em.seurat.ob@meta.data
devTime     <- sample_info$devTime
publication <- sample_info$publication

## Hard check: sample_info rows and counts_mat columns MUST be
## in the exact same order. If this is FALSE, stop immediately
## instead of silently producing wrong groupings downstream.
stopifnot(identical(rownames(sample_info), colnames(counts_mat)))

#### Build a mapping table
## original_id  -> the untouched sample ID (matches counts_mat colnames)
## group_label  -> publication + devTime, used for pseudo-bulk grouping

group_label <- paste(publication, devTime, sep = "_")

sample_mapping <- data.frame(
  original_id      = colnames(counts_mat),   # always derived from the untouched matrix
  sample_id        = rownames(sample_info),
  publication       = publication,
  devTime           = devTime,
  group_label       = group_label,
  stringsAsFactors  = FALSE
)

#### Pseudo-bulk aggregation (SUM of RAW counts per group)
##    Averaging is done on raw counts, NOT on log-normalized data

groups <- unique(sample_mapping$group_label)

pseudobulk_counts <- sapply(groups, function(g) {
  cols <- sample_mapping$original_id[sample_mapping$group_label == g]
  
  # Defensive check: every column name must exist in counts_mat.
  # If this fails, it almost always means counts_mat or
  # sample_mapping got renamed/reordered somewhere upstream.
  stopifnot(all(cols %in% colnames(counts_mat)))
  
  if (length(cols) == 1) {
    counts_mat[, cols]
  } else {
    Matrix::rowSums(counts_mat[, cols, drop = FALSE])
  }
})

colnames(pseudobulk_counts) <- groups
dim(pseudobulk_counts)
saveRDS(pseudobulk_counts,"pseudobulk_counts.rds")

#### CPM normalization of the pseudo-bulk matrix
##    (library-size normalization only; does NOT correct for
##    gene length -- this is CPM, not TPM/RPKM/FPKM)

libsize <- colSums(pseudobulk_counts)

pseudobulk_cpm <- sweep(pseudobulk_counts, 2, libsize, FUN = "/") * 1e6
pseudobulk_cpm <- as.data.frame(pseudobulk_cpm)
pseudobulk_cpm$gene_name <- row.names(pseudobulk_cpm)
row.names(pseudobulk_cpm) <- NULL
dim(pseudobulk_cpm)

#### Gene types
protein_lncRNA_GRCh38 <- readRDS("protein_lncRNA_GRCh38.rds") # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R

pseudobulk_cpm_info <- pseudobulk_cpm %>%
  inner_join(protein_lncRNA_GRCh38 %>%
               dplyr::select(gene_name,intron_class),
             by = "gene_name")

pseudobulk_long <- pseudobulk_cpm_info %>%
  pivot_longer(
    cols      = -c(gene_name, intron_class),
    names_to  = "group_label",
    values_to = "gene_value"
  )

pseudobulk_long <- pseudobulk_long %>%
  separate(col = group_label,sep = "_",into = c("publication","devTime"), remove = FALSE)

devTime_order <- c("Zygote", "E1", "E2", "E3", "E4", "E5", "E6", "E7", 
                   "E8", "E9", "E10", "E12", "E14", "CS7")

pseudobulk_long$devTime <- factor(pseudobulk_long$devTime, levels = devTime_order)


pseudobulk_long <- pseudobulk_long %>%
  mutate(log_value = log10(gene_value + 1),
         publication = factor(publication,
                              levels = rev(c("Tyser et al","Xiang et al","Yanagida et al", 
                                             "Meistermann et al","Petropoulos et al","Yan et al"))),
         intron_class = factor(intron_class,
                               levels = c("long/no_intron_gene","short_GC-rich","short_mixed","short_AT-rich"))
  )
saveRDS(pseudobulk_long,file = "pseudobulk_long.rds")


#### 02. plot Fig4A ####
pseudobulk_long <- readRDS("pseudobulk_long.rds")
pseudobulk_newType <- pseudobulk_long %>%
  mutate(intron_type = ifelse(intron_class == "long/no_intron_gene",
                              "long/no_intron_gene",
                              "short_intron_gene"))
human <- pseudobulk_newType %>%
  transmute(gene = gene_name,
            intron_type,
            group_label,
            source = "Human",
            gene_value) %>%
  distinct()
group_order <- c("Yan et al_Zygote", "Yan et al_E1","Yan et al_E2","Yan et al_E3","Yan et al_E4", "Yan et al_E6" ,
                 "Petropoulos et al_E3","Petropoulos et al_E4", "Petropoulos et al_E5","Petropoulos et al_E6", "Petropoulos et al_E7",
                 "Yanagida et al_E5", "Yanagida et al_E6","Yanagida et al_E7", 
                 "Meistermann et al_E4", "Meistermann et al_E5","Meistermann et al_E6",
                 "Xiang et al_E6","Xiang et al_E7","Xiang et al_E8","Xiang et al_E9","Xiang et al_E10","Xiang et al_E12","Xiang et al_E14",     
                 "Tyser et al_CS7" )

human_plot <- human %>%
  transmute(
    gene,
    sample = group_label,
    intron_type,
    expression = log10(gene_value + 1),
    source) %>%
  group_by(source, sample, intron_type) %>%
  summarise(
    median_exp = median(expression, na.rm=TRUE),
    .groups = "drop") %>%
  mutate(sample = factor(sample,
                         levels = c(devTime_order,
                                    group_order))
         )

paired_df <- human_plot %>%
  select(source, sample, intron_type, median_exp) %>%
  pivot_wider(
    names_from = intron_type,
    values_from = median_exp
  ) %>%
  filter(source == "Human") 

x_position <- c(
  short_intron = 1,
  long_no_intron = 2
)
x_labels <- c(
  short_intron = "Short intron genes",
  long_no_intron = "Long/no-intron genes"
)
#### fig4A
paired_dot <- ggplot(paired_df) +
  geom_segment(
    aes(x = x_position["short_intron"], xend = x_position["long_no_intron"],
        y = short_intron_gene, yend = `long/no_intron_gene`, colour = sample),
    linewidth = 1, alpha = 0.5) +
  geom_point(
    aes(x = x_position["short_intron"], y = short_intron_gene), size = 2, colour = "#E15759") +
  geom_point(
    aes(x = x_position["long_no_intron"], y = `long/no_intron_gene`),
    size = 2, colour = "#4E79A7") +
  scale_x_continuous(breaks = x_position, labels = x_labels[names(x_position)], limits = c(0.8, 2.2)) +
  labs(x = NULL, y = expression("Median "*log[10]*"(expression + 1)"), colour = "Sample") +
  theme_classic(base_size = 13)

paired_dot

ggsave(filename = "paired_dot.pdf",paired_dot,
       width = 7,height = 5)

#### 03. plot fig4B ####
### import GTEx 
protein_lncRNA_GRCh38 <- readRDS("protein_lncRNA_GRCh38.rds") # gene annotations

allTissue_transcriptType_proteinCoding <- readRDS("allTissue_transcriptType_proteinCoding.rds") # derived from short-intron-evolution/figuers/fig1/fig1C.R

allTissue_transcriptType_proteinCoding$gene_id <- sub("\\..*$", "", allTissue_transcriptType_proteinCoding$gene_id)

GTEx_expression_info <- allTissue_transcriptType_proteinCoding %>%
  inner_join(protein_lncRNA_GRCh38 %>%
               dplyr::distinct(gene_id, intron_class,gene_name),
             by = "gene_id"
  ) %>%
  dplyr::distinct(gene_name,gene_id, rowMean, tissue, intron_class)

multi_id_genes <- c("MATR3", "PINX1", "POLR2J3", "SIGLEC5", "TBCE", "TMSB15B")  

gene_level_df <- GTEx_expression_info %>%
  group_by(gene_name, gene_id, tissue, intron_class) %>%
  summarise(rowMean = sum(rowMean, na.rm = TRUE), .groups = "drop") %>% # summaries isoform expression with sum gene expression
  mutate(gene_name = if_else(
    gene_name %in% multi_id_genes,
    paste0(gene_name, "_", gene_id),   
    gene_name
  ))

#### combine GTEx and embryo
embryo_df <- pseudobulk_long %>%
  transmute(
    gene_name,
    sample = group_label,
    source = "Embryo",
    intron_class,
    expression = log10(gene_value + 1)
  )

gtex_df <- gene_level_df %>%
  transmute(
    gene_name,
    sample = tissue,
    source = "Adult (GTEx)",
    intron_class,
    expression = log10(rowMean + 1)
  )

combined_df <- bind_rows(
  embryo_df,
  gtex_df) %>%
  mutate(source = factor(source,
                         levels = c("Embryo","Adult (GTEx)")))

heatmap_df <- combined_df %>%
  group_by(source, sample, intron_class) %>%
  summarise(
    median_exp = median(expression, na.rm=TRUE),
    .groups = "drop"
  )

heatmap_df <- heatmap_df %>%
  group_by(source, sample) %>%
  mutate(
    zscore = scale(median_exp)[,1]
  ) %>%
  ungroup()

remove_list <- c("Meistermann",
                 "Tyser", 
                 "Xiang",
                 "Yan", 
                 "Yanagida")
remove_pattern <- str_c(remove_list, collapse = "|")

wide_df <- heatmap_df %>%
  filter(!str_detect(string = sample,pattern = remove_pattern)) %>%
  dplyr::select(source, sample, intron_class, zscore) %>%
  pivot_wider(names_from = intron_class, values_from = zscore) %>%
  dplyr::select(source, sample,
                `short_GC-rich`,
                short_mixed,
                `short_AT-rich`,
                `long/no_intron_gene`)


mat <- as.matrix(wide_df %>% select(-source, -sample))
rownames(mat) <- wide_df$sample


row_anno <- rowAnnotation(
  Source = wide_df$source,
  col = list(Source = c("Embryo" = "darkseagreen", "Adult (GTEx)" = "thistle4"))
)

max_abs <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs),
                      c("steelblue3", "white", "tomato2"))

#### fig4B
pdf("intron_class_heatmap.pdf", width = 5, height = 6)
Heatmap(
  mat,
  name = "Z-score",
  col = col_fun,
  left_annotation = row_anno,
  row_split = wide_df$source,      
  cluster_row_slices = FALSE,      
  cluster_rows = TRUE,             
  cluster_columns = FALSE,         
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8.2),
  column_title = "Relative abundance of gene classes within each sample\n(z-score computed within sample)"
)
dev.off()