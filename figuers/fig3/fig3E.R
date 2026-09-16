library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)
library(ggpubr)
library(rstatix)

selected_df <- read.csv("human_geneID_RNAexp_selected_df.csv")
selected_tissues <- c("Muscle","Heart","Liver",
                      "Kidney","Testis","Brain","Nerve","Lung")
# color define
color_main <- wes_palette("Darjeeling2", 5, type = "discrete")  # for 0~4
color_extra <- wes_palette("Darjeeling1", 2, type = "discrete") # for -2 (Human), -1 (Mammal)

ageScore_colors <- c(
  "-2" = color_extra[1],
  "-1" = color_extra[2],
  "0"  = color_main[1],
  "1"  = color_main[2],
  "2"  = color_main[3],
  "3"  = color_main[4],
  "4"  = color_main[5],
  "5" = "#CD3278", # for  age 4 GC-rich
  "6" = "cornsilk1" # for age 4 mixed AT-/GC-rich
)

top_short_plot <- ggplot(selected_df, aes(x = tissue, 
                                          y = log_rowMean, 
                                          fill = as.factor(ageScore))) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) +
  scale_fill_manual(values = ageScore_colors,
                    name = "Intron age",
                    labels = c(
                      "-2" = "Human-specific genes",
                      "-1" = "Mammal-specific genes",
                      "0" = "0","1" = "1","2" = "2", "3" = "3",
                      "4" = "4 AT-rich","5" = "4 GC-rich",
                      "6" = "4 AT/GC-mixed"
                    ),
                    guide = guide_legend(title.position = "top", nrow = 9)) +
  theme_linedraw() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank()
  ) +
  labs(x="",y = expression(log[2]*"(TPM + 1)"))

top_short_plot

# statistical analysis
kruskal_short_df <- selected_df %>%
  mutate(ageScore = as.factor(ageScore)) %>%
  group_by(tissue) %>%
  kruskal_test(log_rowMean ~ ageScore) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    p.adj.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    ),
    x = factor(tissue, levels = selected_tissues),
    xend = x,
    group1 = x, 
    group2 = x,
    y.position = 10.5
  )

top_short_plot_sig <- top_short_plot +
  stat_pvalue_manual(
    data = kruskal_short_df,
    label = "p.adj.signif",
    y.position = "y.position",
    tip.length = 0,
    inherit.aes = FALSE
  )


top_short_plot_sig

ggsave(filename = "top_short_plot_sig_geneLevel.pdf", plot = top_short_plot_sig, 
       width = 15, height = 6, units = "in")
