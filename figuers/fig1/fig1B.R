library(stringi)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(webr)
### 1. import the OGEE dataset
OGEE <- read.delim("CSEGs_CEGs_R.txt",header = TRUE, sep = "\t") # download from OGEE V3
## select CEGs	Core-essential genes (Genes essential in 80% or more tested cell lines.)
OGEE_CEGs <- OGEE %>% filter(essentiality=="CEGs")

## 2. import the genes with/without short introns
all_gene_GRCh38 <- readRDS("all_gene_GRCh38.rds") # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R

# genes have protein_coding isoforms > lncRNA isoforms > other_non_coding isoforms
all_gene_filtered <- all_gene_GRCh38 %>%
  distinct(gene_id, gene_name,
           has_intron, intron_class,
           transcript_type ) %>%
  mutate(priority = case_when(
    transcript_type == "protein_coding" ~ 1,
    transcript_type == "lncRNA" ~ 2,
    TRUE ~ 3
  )) %>%
  group_by(gene_id, gene_name) %>%
  filter(priority == min(priority)) %>%
  dplyr::select(-priority) %>%
  ungroup()

## 3. combine with OGEE
gene_intronType_essential <- all_gene_filtered %>% 
  mutate(essentiality=ifelse(gene_id %in% OGEE_CEGs$ensembl,
                             "essential",NA),
         essentiality = ifelse(is.na(essentiality),
                               ifelse(gene_name %in% OGEE_CEGs$gene,
                                      "essential","no"),
                               essentiality
         )
  )

#### check the OGEE gene name ####
### All gene/ensembl name in OGEE_CEGs should be found in gene_intronType ###
mart <- useMart("ensembl", 
                dataset = "hsapiens_gene_ensembl",
                host = "https://mart.ensembl.org")

listFilters(mart)
listAttributes(mart)
lost_ensemblID <- OGEE_CEGs[which(!OGEE_CEGs$ensembl %in% gene_intronType_essential$gene_id),]$gene 

ensemblID <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'external_gene_name',
                   values = lost_ensemblID,
                   mart = mart)

# If still not found --> Manually finding
OGEE_CEGs <- 
  OGEE_CEGs %>% mutate( gene=ifelse(ensembl=="ENSG00000136628","EPRS1",gene),
                        gene=ifelse(ensembl=="ENSG00000276368","H2AC14",gene),
                        gene=ifelse(ensembl=="ENSG00000276410","H2BC3",gene),
                        gene=ifelse(ensembl=="ENSG00000180596","H2BC4",gene),
                        gene=ifelse(ensembl=="ENSG00000274290","H2BC6",gene),
                        gene=ifelse(ensembl=="ENSG00000277224","H2BC7",gene),
                        gene=ifelse(ensembl=="ENSG00000275713","H2BC9",gene),
                        gene=ifelse(ensembl=="ENSG00000278191; ENSG00000110619","CARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000262788; ENSG00000141076","UTP4",gene),
                        gene=ifelse(ensembl=="ENSG00000127054","INTS11",gene),
                        gene=ifelse(ensembl=="ENSG00000115866","DARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000106105","GARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000204628","RACK1",gene),
                        gene=ifelse(ensembl=="ENSG00000196305","IARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000065427","KARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000133706","LARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000166986","MARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000134440","NARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000100138","SNU13",gene),
                        gene=ifelse(ensembl=="ENSG00000149716","LTO1",gene),
                        gene=ifelse(ensembl=="ENSG00000172053","QARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000113643","RARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000022277","RTF2",gene),
                        gene=ifelse(ensembl=="ENSG00000031698","SARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000113407","TARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000154582","ELOC",gene),
                        gene=ifelse(ensembl=="ENSG00000070010","UFD1",gene),
                        gene=ifelse(ensembl=="ENSG00000231945; ENSG00000204394; ENSG00000226589; ENSG00000224264; ENSG00000096171; ENSG00000227686; ENSG00000231116","VARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000140105","WARS1",gene),
                        gene=ifelse(ensembl=="ENSG00000071462","BUD23",gene),
                        gene=ifelse(ensembl=="ENSG00000134684","YARS1",gene)) #


## 3.2 recombine with OGEE
gene_intronType_essential <- all_gene_filtered %>% 
  mutate(essentiality=ifelse(gene_id %in% OGEE_CEGs$ensembl,
                             "essential",NA),
         essentiality = ifelse(is.na(essentiality),
                               ifelse(gene_name %in% OGEE_CEGs$gene,
                                      "essential","no"),
                               essentiality
         )
  ) %>%
  filter(transcript_type != "other_non_coding") %>% # only keep "protein_coding" and "lncRNA" 
  mutate(
    intron_class = ifelse(intron_class == "long/no_intron_gene",
                          "long/no_intron_gene",
                          "short_intron_gene")) 
## 4. fisher exact test
table(gene_intronType_essential$intron_class,
      gene_intronType_essential$essentiality)
# essential    no
# long/no_intron_gene       332 29697
# short_intron_gene         273  6206

data_matrix <- matrix(c(273, 6206, 332, 29697), nrow = 2, byrow = TRUE,
                      dimnames = list(short_intron_gene = c("yes", "no"),
                                      essentiality = c("essential", "no")))


fisher <- fisher.test(data_matrix) # p-value < 2.2e-16; odds ratio 3.934626 
chisq <- chisq.test(data_matrix) # p-value < 2.2e-16


## 5. plot fig1B
# Summarize data for plotting
plot_data <- gene_intronType_essential %>%
  group_by(intron_class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(
    percentage = count / sum(count) * 100,
    lab.pos = cumsum(percentage) - 0.5 * percentage,
    intron_class = factor(
      intron_class,
      levels = c("short_intron_gene",
                 "long/no_intron_gene")
    )
  )


donut <- plot_data %>%
  mutate(intron_class = factor(ifelse(intron_class=="short_intron_gene",
                                      "Short intron\ngenes",
                                      "Long or no\nintron genes"),
                               levels=c("Short intron\ngenes",
                                        "Long or no\nintron genes"))) %>% 
  ggplot(aes(x = 2, y = percentage, fill = intron_class))+
  geom_bar(stat = "identity")+
  coord_polar("y", start = 200) +
  geom_text(aes(x = 1.5, y = lab.pos,
                label = paste(intron_class,"\n",round(percentage,digits = 1),
                              "%", 
                              sep = "")), 
            col = "black", size=6) +
  theme_void() +
  scale_fill_manual(values=c("#FF3030","#607B8B"))+
  xlim(.1,2.5)+
  theme(legend.position="none")
donut

## pie plots
Palette <- c("#FF7F50","#CDCDC1","#B2DFEE")


## short intron
short_intron_gene <- gene_intronType_essential %>% filter(intron_class == "short_intron_gene") %>%
  mutate(transcript_type = str_replace(transcript_type,
                                       pattern = "_",
                                       replacement = "\n"
  ),
  type = paste0(ifelse(essentiality=="essential",
                       "Essential\n",""),
                transcript_type)) %>%
  group_by(type) %>%
  summarise(count = n(), .groups = 'drop')%>%
  arrange(desc(type)) %>%
  mutate(prop = count / sum(count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

pie_short <- 
  ggplot(short_intron_gene, aes(x = "", y = prop, fill = type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=Palette)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none") +
  geom_text(aes(#x = 1.3,
    y = ypos,
    label = paste0(type,"\n",round(prop,digits = 1),"%")), 
    color = "black", size=6)

## long intron
long_intron_gene <- gene_intronType_essential %>% filter(intron_class != "short_intron_gene") %>%
  mutate(transcript_type = str_replace(transcript_type,
                                       pattern = "_",
                                       replacement = "\n"),
         type = paste0(ifelse(essentiality=="essential",
                              "Essential\n",""),
                       transcript_type)) %>%
  group_by(type) %>%
  summarise(count = n(), .groups = 'drop')%>%
  arrange(desc(type)) %>%
  mutate(prop = count / sum(count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

pie_long <- 
  ggplot(long_intron_gene, aes(x = "", y = prop, fill = type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=Palette)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = paste0(type,"\n",round(prop,digits = 1),"%")), color = "black", size=6)



new_plots <- ggarrange(pie_short, donut, pie_long,
                       ncol = 3, nrow = 1)+
  annotate("text", x = 0.5, y = 0.1, label = paste0("Chi-square test p.value < 2.2e-16\nOdds Ratio = ",
                                                    round(fisher$estimate,digits = 2)),
           size = 5,
           hjust = 0.5)

new_plots

ggsave(filename = "Gene_Essentiality.pdf", plot = new_plots,
       width = 8, height = 3.5)
