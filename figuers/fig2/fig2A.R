library(tidyr)
library(stringr)
library(ggtree)
library(ape)
library(aplot)
library(viridis)
library(ggplot2)
library(ggridges)
library(ggstance)
library(biomaRt)
library(dplyr)
library(Biostrings)
# yeast intron
yeast_intron_df <- readRDS("yeast_intron_df.rds") # derived from short-intron-evolution/data_processing/yeast_introns.R
yeast <- yeast_intron_df %>%
  dplyr::transmute(
    species="Saccharomyces cerevisiae",
    seqnames,intron_start,intron_end,intron_length,intron_seqs,intron_GC=intron_gc
  )
# plants intron
plants <- readRDS("combined_gc.rds") # derived from short-intron-evolution/data_processing/plants_introns.R
plants <- plants %>%
  dplyr::select(species,seqnames,intron_start,intron_end,intron_length,intron_seqs,intron_GC)

# animal intron 
all_intron_list <- readRDS("all_intron_list.rds") # derived from short-intron-evolution/data_processing/veterbrates_introns.R
names(all_intron_list) <- c("Homo sapiens","Mus musculus","Gallus gallus", "Chrysemys picta bellii","Xenopus tropicalis","Danio rerio")
animal <-  data.frame()
for (name in names(all_intron_list)) {
  data <- all_intron_list[[name]]
  temp_data <- data.frame(species = name,
                          seqnames = data$seqnames,
                          intron_start = as.numeric(data$intron_start),
                          intron_end = as.numeric(data$intron_end),
                          intron_seqs = data$intron_seq,
                          intron_length = as.numeric(data$intron_width),
                          intron_GC = data$intron_gc)
  animal <- rbind(animal, temp_data) %>% distinct()
}


### combine all species and filter fake introns ####
rbind_df <- bind_rows(yeast, animal, plants) %>%
  mutate(
    species = str_replace(species, " ", "\n"),  
    species = factor(species,
                     levels = c("Physcomitrium\npatens", "Oryza\nsativa", "Arabidopsis\nthaliana",
                                "Saccharomyces\ncerevisiae",
                                "Danio\nrerio", "Xenopus\ntropicalis", "Chrysemys\npicta bellii", "Gallus\ngallus", "Homo\nsapiens", "Mus\nmusculus"))
  )

filtered_intron <- rbind_df %>%
  filter(
    intron_length >= 4,
    substr(intron_seqs, 1, 2) %in% c("GT", "AT"),
    substr(intron_seqs, nchar(intron_seqs) - 1, nchar(intron_seqs)) %in% c("AG", "AC")
  ) %>%
  mutate(log_intron_length = log10(intron_length)) 

n_label <- filtered_intron %>%
  count(species) %>%
  mutate(label = paste0("n = ", scales::comma(n)))

#### density plot of intron length and GC  ####
species_colors <- c("palegreen4", "darkolivegreen3", "lavenderblush4")

length <- filtered_intron %>%
  ggplot(aes(x = log_intron_length, fill = species)) +
  geom_density(alpha = .4) +
  #geom_histogram(aes(y = after_stat(count)), bins = 50, alpha = .4)+
  facet_wrap(~species, nrow = 1) +
  geom_text(
    data = n_label,
    aes(label = label, x = Inf, y = Inf),
    hjust = 1.1, vjust = 1.5,
    size = 3, inherit.aes = FALSE
  ) +
  theme_linedraw()+
  #scale_fill_manual(values = species_colors) +
  theme(
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, color = "black", face = "bold.italic")
  ) +
  labs(x = expression(log[10]*"(Intron Width)"), y = "Density")

cbPalette <- c("#E69F00", "#56B4E9")
GC <- filtered_intron  %>%
  mutate(intron_type = ifelse(intron_length <= 100, "<=100", ">100")) %>%
  ggplot(aes(x = intron_GC, fill = intron_type)) +
  geom_density(alpha = .3) +
  facet_wrap(~species, nrow = 1) +
  scale_fill_manual(values = cbPalette) +
  theme_linedraw() +
  theme(legend.position = "inside",
        legend.background = element_blank(),
        legend.position.inside = c(0.94, 0.7),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_blank()
  ) +
  labs(x = "Intron GC Content", y = "Density", fill = "Intron Type")

#### plot tree ####
tree <- read.tree("tree_wt_plants.txt")

tree_wt_plants <- ggtree(tree) +
  geom_label(
    aes(label = label),
    data = function(x) dplyr::filter(x, !isTip, !is.na(label), label != ""),
    size = 4,
    label.size = 0.5,
    fill = "white"
  ) +
  layout_dendrogram()

#### combine tree + density ####
all_merged <- tree_wt_plants/length/GC
ggsave(filename = "all_merged.pdf",all_merged,
       width = 15, height = 7)