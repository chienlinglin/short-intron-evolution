library(DBI)
library(RSQLite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(readr)

#### import result
result_all_tissues <- read_csv("aging_MSR_result_all_tissues.csv")

result_all_tissues <- result_all_tissues %>%
  mutate(
    tissue = tissue %>%
      sub("_misspliced$", "", .) %>%     
      tolower() %>%  
      tools::toTitleCase() 
  )

forest_df <- result_all_tissues %>%
  mutate(
    sig = p_adj < 0.05,
    p_adj_label = ifelse(
      p_adj < 2.2e-16,
      2.2e-16,
      p_adj),
    neglog10_p_adj = -log10(p_adj_label)
  )

#### plot fig4E ####
max_neglog10 <- max(forest_df$neglog10_p_adj)
min_neglog10 <- min(forest_df$neglog10_p_adj)
order_by_allSites <- c("Salivary_gland", "Blood", "Blood_vessel", "Nerve", "Pancreas",
                       "Muscle", "Small_intestine", "Stomach",
                       "Heart", "Esophagus", "Spleen", "Thyroid",
                       "Colon", "Adrenal_gland", "Lung", "Adipose_tissue",
                       "Brain", "Skin")

shortGCAT_long_increasingMSR_odds_plot <- forest_df %>%
  filter(comparison %in% c("long vs short AT-rich","long vs short GC-rich")) %>%
  mutate(tissue = factor(tissue,
                         levels = rev(order_by_allSites)),
         comparison = factor(comparison,
                             levels = c("long vs short GC-rich",
                                        "long vs short AT-rich"))) %>% 
  ggplot(aes(x = OR, y = tissue, color = sig, size = neglog10_p_adj)) +
  facet_wrap( ~ comparison, ncol = 2) +
  geom_point(alpha = 0.85) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    width = 0.2,alpha = 0.85,
    size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  scale_color_manual(
    values = c("TRUE"  = "#D55E00",
               "FALSE" = "grey60"),
    labels = c("TRUE"  = "Sig. (p.adj < 0.05)",
               "FALSE" = "Not sig. (p.adj ≥ 0.05)"),
    name   = ""
  ) +
  scale_size_continuous(
    name   = expression(-log[10]("p.adj")),
    range  = c(1, 4),
    breaks = c(min_neglog10, 1.30103, 5, max_neglog10),
    labels = c(round(min_neglog10,digits = 2), "1.3", "5" , "> 15.66")
  ) +
  labs(
    x = "Odds ratio (short vs long intron)",
    y = NULL,
    color = ""
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",    
        legend.box.just = "left")

shortGCAT_long_increasingMSR_odds_plot

ggsave(filename = "shortGCAT_long_increasingMSR_odds_plot.pdf",
       plot = shortGCAT_long_increasingMSR_odds_plot,width = 6, height = 5)