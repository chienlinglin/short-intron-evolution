library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggplot2)
### genecode v38 version
# GC
sheets <- getSheetNames("human_ordersGC.xlsx")
human_ordersGC <- setNames(
  lapply(sheets, 
         function(x) {
           read.xlsx("human_ordersGC.xlsx", sheet = x)
           }), 
  sheets)

# AT
sheets <- getSheetNames("human_ordersAT.xlsx")
human_ordersAT <- setNames(
  lapply(sheets, 
         function(x) {
           read.xlsx("human_ordersAT.xlsx", sheet = x)
           }),
  sheets
)
### functions for processing df ####
make_cluster_plot_df <- function(orders, cutoff = 20) {
  
  if (!all(c("observed_raw", "shuffle_mean_df") %in% names(orders))) {
    stop("orders must contain observed_raw and shuffle_mean_df")
  }
  
  raw_df <- rbind(
    orders$shuffle_mean_df,
    orders$observed_raw %>%
      transmute(
        value = delta_order,
        type = "observed"
      )
  )
  
  raw_df_mod <- raw_df %>%
    mutate(
      distance_group = ifelse(round(value, 0) >= cutoff,
                              paste0("≥", cutoff),
                              as.character(round(value, 0))),
      distance_group = factor(
        distance_group,
        levels = c(as.character(1:(cutoff-1)), paste0("≥", cutoff))
      ),
      type = ifelse(type == "Shuffled_mean",
                    paste0("Shuffled background (n=1000)"),
                    "Observed short intron positions")
    )
  
  return(raw_df_mod)
}
make_species_df <- function(orders, species_name) {
  df <- make_cluster_plot_df(orders)
  df$species <- species_name
  df$ks_label <- colnames(orders$ks_text)
  return(df)
}
#### prepare for ploting ####
humanGC_df      <- make_species_df(human_ordersGC,     "Human")

humanAT_df      <- make_species_df(human_ordersAT,     "Human")

humanGC_df2 <- humanGC_df %>%
  mutate(intron_type = "GC-rich")
humanAT_df2 <- humanAT_df %>%
  mutate(intron_type = "AT-rich")

human_combined_df <- bind_rows(humanGC_df2, humanAT_df2) %>%
  mutate(intron_type = factor(intron_type, levels = c("GC-rich", "AT-rich")))

ks_combined_df <- human_combined_df %>%
  distinct(intron_type, species, ks_label) %>%
  mutate(x = Inf, y = Inf)



#### plots ####
human_combined_plot <- human_combined_df %>%
  filter(!is.na(distance_group)) %>%
  ggplot(aes(x = distance_group, fill = type)) +
  geom_bar(color = "black", alpha = 0.5, position = "identity", lwd = 0.1) +
  scale_fill_manual(values = c("#FF6A6A", "#2F4F4F")) +
  facet_grid(intron_type ~ ., scales = "free_y",
             labeller = labeller(intron_type = c(
               "GC-rich" = "GC-rich short",
               "AT-rich" = "AT-rich short"
             ))
             
  ) +
  labs(fill = "",
       x = "Number of exons to next GC-rich or AT-rich short intron",
       y = "Count"
  ) +
  geom_text(
    data = ks_combined_df,
    aes(label = ks_label, x = x, y = y),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    size = 4, color = "black"
  ) +
  theme_linedraw() +
  theme(
    text = element_text(size = 13),
    legend.position = "inside",
    legend.background = element_blank(),
    legend.position.inside = c(0.65, 0.75),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    panel.grid.major         = element_line(color = "grey95", linewidth = 0.4),
    panel.grid.minor         = element_blank(),
    strip.background         = element_blank(),
    strip.text               = element_text(size = 12, color = "black", face = "bold")
  )
#### save plots ####
ggsave(filename = "human_GCAT.pdf",
       human_combined_plot,
       width = 5, height = 5
)