library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)

#### reference: Dvinge, Heidi et al. ####  
#### “RNA components of the spliceosome regulate tissue- and cancer-specific alternative splicing.” 
#### Genome research vol. 29,10 (2019): 1591-1604. doi:10.1101/gr.246678.118
### excel: Supplemental_Table_S1.xlsx
# Delta PSI values compared to scrambled control oligo for snRNA knockdown within each cell line.
# "-": Splicing events that were either not detected or not statistically significant.
#### they aligned RNA-seq with hg19 ####
#### import excel ####
alternative_events <- read_excel("Supplemental_Table_S1.xlsx",skip = 3) 

muf7_alt <- alternative_events %>%
  dplyr::select(splicing_category, gene_ID, gene_HGNC, event_location,
                starts_with("MCF7")) %>%
  distinct()

hela_alt <- alternative_events %>%
  dplyr::select(splicing_category, gene_ID, gene_HGNC, event_location,
                starts_with("HeLa")) %>%
  distinct()
#### focus on Retained intron
muf7_RI <- muf7_alt %>%
  filter(splicing_category=="Retained intron") %>% 
  distinct() %>%
  separate(col = "event_location",
           into = c("chr","start","end","strand"),
           sep = ":",remove = FALSE
  ) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end),
         intron_length = end-start+1,
         intron_group = case_when(
           intron_length <= 100 ~ "<=100",
           intron_length > 100 ~ ">100"),
         cell_type = "MCF7"
  )

muf7_df_long <- muf7_RI %>%
  pivot_longer(
    cols = c(MCF7.U1, MCF7.U2, MCF7.U4, MCF7.U6),
    names_to = "KD",
    values_to = "delta_PSI"
  ) %>%
  mutate(
    KD = str_remove(KD, "^MCF7\\."),
    delta_PSI = na_if(delta_PSI, "-"),
    delta_PSI = as.numeric(delta_PSI),  
    intron_group = factor(intron_group,
                          levels = c("<=100",">100"))
  )

hela_RI <- hela_alt %>%
  filter(splicing_category=="Retained intron") %>% 
  distinct() %>%
  separate(col = "event_location",
           into = c("chr","start","end","strand"),
           sep = ":",remove = FALSE
  ) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end),
         intron_length = end-start+1,
         intron_group = case_when(
           intron_length <= 100 ~ "<=100",
           intron_length > 100 ~ ">100"),
         cell_type = "HeLa"
  )

hela_df_long <- hela_RI %>%
  pivot_longer(
    cols = c(HeLa.U1, HeLa.U2, HeLa.U4, HeLa.U6),
    names_to = "KD",
    values_to = "delta_PSI"
  ) %>%
  mutate(
    KD = str_remove(KD, "^HeLa\\."),
    delta_PSI = na_if(delta_PSI, "-"),
    delta_PSI = as.numeric(delta_PSI),  
    intron_group = factor(intron_group,
                          levels = c("<=100",">100"))
  )
all_df <- bind_rows(hela_df_long,
                    muf7_df_long)
#### ECDF plot vs intron size ####
df_marked <- all_df %>%
  group_by(KD,cell_type) %>%
  mutate(
    retention = delta_PSI > 0,
    spliced = delta_PSI < 0
  ) %>%
  ungroup()
df_bg <- df_marked %>%
  mutate(group = "Background")
df_retention<- df_marked %>%
  filter(retention) %>%
  mutate(group = "Retention\n(delta PSI > 0)") 
df_spliced <- df_marked %>%
  filter(spliced) %>%
  mutate(group = "Spliced\n(delta PSI < 0)")

df_ecdf <- bind_rows(df_bg, 
                     df_retention,
                     df_spliced) %>%
  mutate(group = factor(group,
                        levels = c("Retention\n(delta PSI > 0)",
                                   "Spliced\n(delta PSI < 0)" ,
                                   "Background"
                        )))

ecdf_plot_len_splicing <- df_ecdf %>%
  filter(intron_length >= 30) %>%
  ggplot(aes(x = log10(intron_length), color = group)) +
  stat_ecdf(
    aes(fill = group),
    geom = "area",
    alpha = 0.1,
    linewidth = 0
  ) +
  stat_ecdf(
    linewidth = 0.7,
    alpha = 0.8
  ) +
  facet_grid(cell_type ~ KD) +
  theme_linedraw(base_size = 14) +
  labs(
    x = expression(log[10]*"(Length of alternative spliced introns)"),
    y = "Cumulative proportion (ECDF)",
    color = "",
    fill = ""
  ) +
  
  ## Publication-ready colors (Nature-friendly)
  scale_color_manual(
    values = c(
      "Retention\n(delta PSI > 0)" = "#D55E00",
      "Spliced\n(delta PSI < 0)"   = "#1E90FF",
      "Background" = "grey66"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Retention\n(delta PSI > 0)" = "#D55E00",
      "Spliced\n(delta PSI < 0)"   = "#1E90FF",
      "Background" = "grey66"
    )
  ) +
  
  theme(
    legend.position = "bottom",
    #legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

ecdf_plot_len_splicing



### ks ###
library(purrr)
library(broom)

ks_results <- df_ecdf %>%
  filter(intron_length >= 30) %>%
  group_by(cell_type, KD) %>%  
  summarise(
    # sample size
    n_ret = sum(group == "Retention\n(delta PSI > 0)"),
    n_spl = sum(group == "Spliced\n(delta PSI < 0)"),
    n_bg  = sum(group == "Background"),
    
    # KS tests
    ks_ret_bg = list(ks.test(
      intron_length[group == "Retention\n(delta PSI > 0)"],
      intron_length[group == "Background"]
    )),
    ks_spl_bg = list(ks.test(
      intron_length[group == "Spliced\n(delta PSI < 0)"],
      intron_length[group == "Background"]
    )),
    ks_res_spl = list(ks.test(
      intron_length[group == "Retention\n(delta PSI > 0)"],
      intron_length[group == "Spliced\n(delta PSI < 0)"]
    )),
    .groups = "drop"
  ) %>%
  mutate(
    # p-values
    p_ret_bg = map_dbl(ks_ret_bg, ~ .x$p.value),
    p_spl_bg = map_dbl(ks_spl_bg, ~ .x$p.value),
    p_ret_spl = map_dbl(ks_res_spl, ~ .x$p.value),
    
    # KS D statistics
    D_ret_bg = map_dbl(ks_ret_bg, ~ .x$statistic),
    D_spl_bg = map_dbl(ks_spl_bg, ~ .x$statistic),
    D_ret_spl = map_dbl(ks_res_spl, ~ .x$statistic)
  ) %>%
  select(-ks_ret_bg, -ks_spl_bg, -ks_res_spl)

df_pvals <- ks_results %>%
  mutate(label = paste0(
    "R=", n_ret, ", S=", n_spl, 
    ", B=", n_bg, 
    "\n",
    # "KS D:\n",
    # "Retention vs Background = ", signif(D_ret_bg, 2), "\n",
    # "Spliced vs Background = ", signif(D_spl_bg, 2), "\n",
    # "Retnetion vs Spliced =", signif(D_ret_spl,2),
    "KS p-value:\n",
    "R vs B= ", signif(p_ret_bg, 2), "\n",
    "S vs B= ", signif(p_spl_bg, 2), "\n",
    "R vs S =", signif(p_ret_spl,2)
  ))

ecdf_plot_len_splicing_stat <- ecdf_plot_len_splicing +
  geom_text(
    data = df_pvals,
    aes(
      x = 6 , y = 0.05,   
      label = label
    ),
    color = "black",
    size = 2.5,
    inherit.aes = FALSE,
    hjust=1,
    vjust=0
  )

ecdf_plot_len_splicing_stat

ggsave(filename = "ecdf_plot_len_splicing_stat.pdf",
       ecdf_plot_len_splicing_stat,width = 9.3,height = 4.5)
#### intron retention (delta_PSI > 0) X-KD, Y-log10(intron len) ####
my_comp <- list(
  c("U1", "U2"),
  c("U2", "U4"),
  c("U4", "U6"),
  c("U1", "U4"),
  c("U2", "U6"),
  c("U1", "U6")
)

#### plot fig5E ####
RI_boxplot <- all_df %>%
  filter(delta_PSI > 0 &
           intron_length > 30) %>%
  ggplot(aes(x = KD, y = log10(intron_length), fill = KD)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA,alpha=0.3) +
  facet_grid(~cell_type) +
  stat_compare_means(
    comparisons = my_comp,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(5.3, 5.3, 5.3, 
                5.6, 6.0, 6.3)  
  ) +
  theme_bw(base_size = 14) +
  labs(x = "KD", 
       y = expression(log[10]*"(Length of retained introns)"))+
  theme(legend.position = "none",
        #strip.background = element_rect(fill = "grey90"),
        #strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
  )

RI_boxplot

ggsave(filename = "RI_boxplot.pdf",
       RI_boxplot,width = 4.3,height = 4.2)