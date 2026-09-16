library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(SummarizedExperiment)
library(ggpubr)

# mouse embryo dataset download from doi: 10.1038/s41556-026-01971-3. Epub 2026 Jun 1.
mouse_embryo_mean_expression_CPM <- read.csv("download.cncb.ac.cn/OMIX/OMIX014646/mouse_embryo_mean_expression_CPM.csv")
GRCm38_mm10_mouse_geneAnnotation <- readRDS("GRCm38_mm10_mouse_geneAnnotation.rds") # refer to mouse_gene_annotation.R
#### embryo stages ####
#### filter 0 counts in all stages ####
stage_cols <- setdiff(colnames(mouse_embryo_mean_expression_CPM), "gene")
mouse_embryo_mean_expression_CPM_filtered <- mouse_embryo_mean_expression_CPM[
  rowSums(mouse_embryo_mean_expression_CPM[, stage_cols]) > 0, ]

rnaExp_geneAnno <- mouse_embryo_mean_expression_CPM_filtered %>%
  inner_join(GRCm38_mm10_mouse_geneAnnotation,by=c("gene"="gene_name")) %>%
  filter(gene_type %in% c("protein_coding","lincRNA"))

rnaExp_geneAnno_long <- rnaExp_geneAnno %>%
  pivot_longer(
    cols      = -c(gene, gene_id, gene_type, 
                   has_intron, any_short, all_long,
                   intron_class),
    names_to  = "group_label",
    values_to = "gene_value"
  )


devTime_order <- c("E1.5", "E2.0",  "E2.5",  "E3.0",  "E3.5",  "E4.0",  "E4.5",  "E5.0",  "E5.5",
                   "E6.0",  "E6.5", "E7.0",  "E7.5",  "E7.75", "E8.0", "E8.5",  "E9.0",  "E9.5",
                   "E10.0", "E10.5", "E11.0", "E11.5", "E12.0", "E12.5", "E13.0", "E13.5", "E14.0", "E14.5",
                   "E15.0", "E15.5","E16.0", "E16.5", "E17.0", "E17.5", "E18.0", "E18.5", "E19.0")

rnaExp_geneAnno_long$group_label <- factor(rnaExp_geneAnno_long$group_label, levels = devTime_order)


rnaExp_geneAnno_long <- rnaExp_geneAnno_long %>%
  mutate(log_value = log10(gene_value + 1),
         intron_class = factor(intron_class,
                               levels = c("short_GC-rich","short_mixed","short_AT-rich","long/no_intron_gene")))
saveRDS(rnaExp_geneAnno_long,file = "mouse_embryo.rds")

#### only embryo ####
embryo_df <- rnaExp_geneAnno_long %>%
  dplyr::transmute(
    gene, intron_class, group_label, gene_value, log_value,
    source = "Mouse embryo\n(Cao et al.)"
  )


combine_df <- embryo_df %>%
  mutate(group_label = factor(group_label,
                              levels = c(devTime_order))
  )


summarise_df <- combine_df %>%
  group_by(group_label, intron_class) %>%
  summarise(
    median_exp = median(log_value, na.rm=TRUE),
    .groups = "drop"
  )

summarise_df <- summarise_df %>%
  group_by(group_label) %>%
  mutate(
    zscore = scale(median_exp)[,1]
  ) %>%
  ungroup()


wide_df <- summarise_df %>%
  dplyr::select(group_label, intron_class, zscore) %>%
  pivot_wider(names_from = intron_class, values_from = zscore)

mat <- as.matrix(wide_df %>% dplyr::select(-group_label))
rownames(mat) <- wide_df$group_label


max_abs <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("steelblue3", "white", "tomato2"))

source_lookup <- combine_df %>%
  distinct(group_label, source) 

row_source <- source_lookup$source[match(rownames(mat), source_lookup$group_label)]

row_anno <- rowAnnotation(
  Source = row_source,
  col = list(Source = c("Mouse embryo\n(Cao et al.)" = "darkseagreen"))
)


pdf("mouse_embryo_only_heatmap.pdf", width = 5.5, height = 7)
Heatmap(
  mat,
  name = "RNA expression\n(Z−score by row)",
  col = col_fun,
  #left_annotation = row_anno,
  row_split = row_source,
  cluster_row_slices = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.1),
  row_names_gp = grid::gpar(fontsize = 8.2)
)
dev.off()
