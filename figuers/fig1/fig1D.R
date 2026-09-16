library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(biomaRt)
library(org.Hs.eg.db)
# RNA stability in HeLa download from GSE117299
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE117299&format=file&file=GSE117299%5FStability%2Etxt%2Egz
HeLaStability <- read.delim("GSE117299_Stability.txt") %>%
  filter(Sample=="HeLa")

# transform symbol with ensg
mart <- useMart("ensembl", 
                dataset = "hsapiens_gene_ensembl",
                host = "https://mart.ensembl.org")


gene_id <- str_replace(unique(all_predicted_hit_endo_m6a$gene_ids), "\\..*", "")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = gene_id,  
  mart = mart
)

simplifed_gene_info <- gene_info %>% dplyr::select(ensembl_gene_id,gene_biotype,hgnc_symbol)

endo_m6a_status <- left_join(all_predicted_hit_endo_m6a,simplifed_gene_info,
                             by=c("gene_ids_clean"="ensembl_gene_id"))

#### merge with half-lives
m6a_status_halfLives <- inner_join(HeLaStability,endo_m6a_status,
                                   by=c("Symbol"="hgnc_symbol"))

saveRDS(m6a_status_halfLives,"half-lives/m6a_status_halfLives.rds")

gene_class <- m6a_status_halfLives %>%
  group_by(Symbol) %>%
  summarise(
    m6a_annotation = if_else(any(annotation == "endo_Met"), "endo_m6A_gene", "endo_unm6A_gene"),
    gene_annotation = if_else(any(intron_anno_light == "flanking\nwith\nshort\nintron(s)"), 
                              "short intron(s)\ncontaining", 
                              "no\nshort intron"),
    Half = dplyr::first(Half)
  )


# plot the half-lives vs gene types
table(gene_class$gene_annotation)

plot_df <- gene_class %>%
  mutate(gene_annotation=factor(ifelse(gene_annotation=="short intron(s)\ncontaining",
                                       "with\nshort introns\n(n=2966)",
                                       "without\nshort introns\n(n=8443)"),
                                levels = c("without\nshort introns\n(n=8443)",
                                           "with\nshort introns\n(n=2966)"))
  ) %>%
  filter(Half!=0) 

shortIntronGene_halfLives <- plot_df %>%
  ggplot(aes(x = gene_annotation, y = log2(Half), fill = gene_annotation)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(
    values = c(
      "with\nshort introns\n(n=2966)" = "#F8766D" ,
      "without\nshort introns\n(n=8443)" = "#00BFC4" 
    )
  ) +
  coord_flip() +
  stat_compare_means(method = "wilcox.test", label.x = 2.2, label.y = 6) +
  labs(x = "Protein-coding/lncRNA\ngenes",
       y = expression("mRNA half-life ("*log[2]*"(hours))"),) +
  theme_linedraw(base_size = 16) +
  geom_hline(
    aes(yintercept = median(log2(Half))),
    data = subset(plot_df, gene_annotation == "without\nshort introns\n(n=8443)"),
    color = "gray20", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3))

shortIntronGene_halfLives

ggsave(filename = "shortIntronGene_halfLives.pdf",plot = shortIntronGene_halfLives,
       device = "pdf",width = 8,height = 3)