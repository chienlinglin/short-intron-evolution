library(dplyr)
library(tidyr)
library(stringr)
library(rtracklayer)
library(ggplot2)
library(ggpubr)
library(purrr)
library(furrr)
library(progressr)
library(future.apply)
# read the gene file
gene_intronType <- readRDS("protein_lncRNA_GRCh38.rds") %>% # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R
  mutate(seqnames = factor(seqnames,
                           levels = paste0("chr", c(1:22, "X", "Y", "M"))))
# define order index
compute_observed_counts <- function(data, target_class, collapse_threshold = 21L) {
  
  distance_df <- data %>%
    arrange(seqnames, start) %>%
    group_by(seqnames) %>%
    mutate(gene_order = row_number()) %>%
    ungroup() %>%
    filter(intron_class == target_class) %>%
    group_by(seqnames) %>%
    arrange(seqnames, gene_order) %>%
    mutate(distance_to_next = lead(gene_order) - gene_order) %>%
    ungroup()
  
  distance_df %>%
    select(seqnames, gene_id, distance_to_next, intron_class) %>%
    filter(!is.na(distance_to_next)) %>%
    group_by(seqnames, distance_to_next) %>%
    summarise(observed_count = n(), .groups = "drop") %>%
    mutate(
      distance_to_next = if_else(
        distance_to_next >= collapse_threshold,
        collapse_threshold,
        as.integer(distance_to_next)
      )
    ) %>%
    group_by(seqnames, distance_to_next) %>%
    summarise(observed_count = sum(observed_count), .groups = "drop") %>%
    arrange(seqnames, distance_to_next)
}
observed_counts_AT    <- compute_observed_counts(gene_intronType, "short_AT-rich")
observed_counts_GC    <- compute_observed_counts(gene_intronType, "short_GC-rich")

#### simulation ####
options(future.globals.maxSize = 4 * 1024^3)
options(future.wait.timeout = 600)
gene_intronType <- readRDS("protein_lncRNA_GRCh38.rds") %>% # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R
  mutate(seqnames = factor(seqnames,
                           levels = paste0("chr", c(1:22, "X", "Y", "M"))))

shortGC_gene_ids <- gene_intronType %>%
  filter(intron_class == "short_GC-rich") %>%
  pull(gene_id) %>% unique()
shortAT_gene_ids <- gene_intronType %>%
  filter(intron_class == "short_AT-rich") %>%
  pull(gene_id) %>% unique()

# all genes information
all_genes <- gene_intronType %>%
  select(seqnames, gene_id, start) %>% 
  arrange(seqnames, start) %>%
  group_by(seqnames) %>%
  mutate(original_order = row_number()) %>%
  ungroup()
# shuffle gene_order，
# calcualte the distance_to_next of short gene  
set.seed(12)
plan(multisession, workers = 2)
safe_simulate_distance <- function(i, all_genes, short_gene_ids) {
  shuffled <- all_genes %>%
    group_by(seqnames) %>%
    mutate(shuffled_order = sample(row_number())) %>%
    ungroup() %>%
    filter(gene_id %in% short_gene_ids) %>%
    arrange(seqnames, shuffled_order) %>%
    group_by(seqnames) %>%
    mutate(distance_to_next = lead(shuffled_order) - shuffled_order) %>%
    ungroup() %>%
    select(seqnames, gene_id, distance_to_next) %>%
    mutate(sim_id = i)
  
  return(shuffled)
}

handlers(global = TRUE)
handlers("txtprogressbar")
with_progress({
  simulated_GC_distances <- future_map_dfr(
    1:1000,
    ~ safe_simulate_distance(.x, all_genes = all_genes, short_gene_ids = shortGC_gene_ids),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE  
  )
})

with_progress({
  simulated_AT_distances <- future_map_dfr(
    1:1000,
    ~ safe_simulate_distance(.x, all_genes = all_genes, short_gene_ids = shortAT_gene_ids),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE  
  )
})
# saveRDS(simulated_GC_distances,"simulated_GC_distances.rds")
# saveRDS(simulated_AT_distances,"simulated_AT_distances.rds")

all_ids <- 1:1000
simGC_counts_per_run <- simulated_GC_distances %>%
  filter(!is.na(distance_to_next)) %>%
  group_by(sim_id, seqnames, distance_to_next) %>%
  summarise(count = n(), .groups = "drop")%>%
  group_by(seqnames, distance_to_next) %>%
  complete(sim_id = all_ids, fill = list(count = 0)) %>%
  ungroup()
simAT_counts_per_run <- simulated_AT_distances %>%
  filter(!is.na(distance_to_next)) %>%
  group_by(sim_id, seqnames, distance_to_next) %>%
  summarise(count = n(), .groups = "drop")%>%
  group_by(seqnames, distance_to_next) %>%
  complete(sim_id = all_ids, fill = list(count = 0)) %>%
  ungroup()

# calculate mean and sd across sim_id
simGC_counts_summary <- simGC_counts_per_run %>%
  group_by(seqnames, distance_to_next) %>%
  summarise(
    mean_count = round(mean(count), digits = 2),
    mean_sd = sd(count),
    .groups = "drop"
  )
simAT_counts_summary <- simAT_counts_per_run %>%
  group_by(seqnames, distance_to_next) %>%
  summarise(
    mean_count = round(mean(count), digits = 2),
    mean_sd = sd(count),
    .groups = "drop"
  )

# saveRDS(simGC_counts_summary,"simGC_counts_summary.rds")
# saveRDS(simAT_counts_summary,"simAT_counts_summary.rds")

simGC_counts_summary_collapsed <- simGC_counts_summary %>%
  mutate(distance_to_next = ifelse(distance_to_next >= 21, 21L, distance_to_next)) %>%
  group_by(seqnames, distance_to_next) %>%
  summarise(mean_count = sum(mean_count), .groups = "drop")
simAT_counts_summary_collapsed <- simAT_counts_summary %>%
  mutate(distance_to_next = ifelse(distance_to_next >= 21, 21L, distance_to_next)) %>%
  group_by(seqnames, distance_to_next) %>%
  summarise(mean_count = sum(mean_count), .groups = "drop")
#### combine ####
compareGC_counts <- full_join(
  observed_counts_GC,
  simGC_counts_summary_collapsed,
  by = c("seqnames", "distance_to_next")) %>%
  replace_na(list(
    observed_count = 0,
    mean_count = 0
  ))

compareAT_counts <- full_join(
  observed_counts_AT,
  simAT_counts_summary_collapsed,
  by = c("seqnames", "distance_to_next")) %>%
  replace_na(list(
    observed_count = 0,
    mean_count = 0
  ))
#### ks ####
run_ks_test <- function(data, p_adjust_method = "BH") {
  data %>%
    group_by(seqnames) %>%
    nest() %>%
    mutate(
      ks_p_value = map_dbl(data, function(df) {
        obs <- rep(df$distance_to_next, df$observed_count)
        rnd <- rep(df$distance_to_next, round(df$mean_count))
        if (length(obs) > 1 && length(rnd) > 1) {
          ks.test(obs, rnd)$p.value
        } else {
          NA_real_
        }
      })
    ) %>%
    ungroup() %>%
    mutate(
      ks_p_adj = p.adjust(ks_p_value, method = p_adjust_method),
      p_label  = paste0("KS adj.p = ", signif(ks_p_adj, 2))
    ) %>%
    select(seqnames, ks_p_value, ks_p_adj, p_label)
}
ks_AT    <- run_ks_test(compareAT_counts)
ks_GC    <- run_ks_test(compareGC_counts)

ks_both <- rbind(ks_GC %>% mutate(Gene_type = "GC-rich short"),
                 ks_AT %>% mutate(Gene_type = "AT-rich short")
) %>%
  mutate(Gene_type = factor(Gene_type,
                            levels = c("GC-rich short",
                                       "AT-rich short")))

ks_both$seqnames <- factor(ks_both$seqnames, 
                           levels = paste0("chr", c(1:22, "X", "Y", "M")))

#### plot ####
plot_GCdata <- compareGC_counts %>%
  left_join(ks_GC, by = "seqnames") %>%
  mutate(factor(seqnames, 
                levels = paste0("chr", c(1:22, "X", "Y", "M"))),
         Gene_type = "GC-rich short"
  )
plot_ATdata <- compareAT_counts %>%
  left_join(ks_AT, by = "seqnames") %>%
  mutate(factor(seqnames, 
                levels = paste0("chr", c(1:22, "X", "Y", "M"))),
         Gene_type = "AT-rich short"
  )

plot_df <- rbind(plot_GCdata, plot_ATdata) %>%
  mutate(Gene_type = factor(Gene_type,
                            levels = c("GC-rich short",
                                       "AT-rich short")))

plot_distance_distribution <- function(plot_data, ks_data,
                                       real_color     = "#ff7f0e",
                                       shuffle_color  = "#7D7D7D",
                                       plot_title = "") {
  ggplot(plot_data, aes(x = distance_to_next - 1)) +
    geom_col(aes(y = observed_count, fill = "Real order"),
             position = "identity", alpha = 0.6) +
    geom_col(aes(y = mean_count, fill = "Shuffled 1000\ntimes mean"),
             position = "identity", alpha = 0.4) +
    geom_text(
      data          = ks_data,
      aes(x = Inf, y = Inf, label = p_label),
      hjust         = 1.1,
      vjust         = 1.2,
      inherit.aes   = FALSE,
      size          = 3
    ) +
    facet_grid(Gene_type ~ seqnames, scales = "free_y") +
    scale_fill_manual(values = c(
      "Real order"             = real_color,
      "Shuffled 1000\ntimes mean" = shuffle_color
    )) +
    labs(
      x    = paste0(plot_title," gene ordinal distance"),
      y    = "Count",
      fill = ""
    ) +
    theme_linedraw() +
    theme(
      legend.position          = "top",
      text                     = element_text(size = 14),
      axis.text.x              = element_text(size = 10),
      panel.grid.major         = element_line(color = "grey95", linewidth = 0.4),
      panel.grid.minor         = element_blank(),
      strip.background         = element_blank(),
      strip.text               = element_text(size = 12, color = "black", face = "bold")
    )
}

selected_chrom <- c(paste0("chr", c(6,8,11,16,17,19)))

GC_AT_barplot_real_randon_sele <- plot_distance_distribution(plot_df %>%
                                                               filter(seqnames %in% selected_chrom), 
                                                             ks_both %>%
                                                               filter(seqnames %in% selected_chrom)
)

GC_AT_barplot_real_randon_other <- plot_distance_distribution(plot_df %>%
                                                                filter(!seqnames %in% selected_chrom),    
                                                              ks_both %>%
                                                                filter(!seqnames %in% selected_chrom))

ggsave(GC_AT_barplot_real_randon_sele,
       filename = "GC_AT_barplot_real_randon_sele.pdf",
       width = 12, height = 4.5)
ggsave(GC_AT_barplot_real_randon_other,
       filename = "GC_AT_barplot_real_randon_other.pdf",
       width = 20, height = 4.5)