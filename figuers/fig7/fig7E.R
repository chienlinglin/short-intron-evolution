library(dplyr)
library(tidyr)
library(ggplot2)

###### MUS81 #######
group_sample <- read.csv("alt_exon_ratio.csv")

summary_df <- group_sample %>%
  group_by(groups, altMidExon_isoform) %>%
  summarise(
    mean_ratio = mean(percentage),
    sd_ratio   = sd(percentage),
    n          = n(),
    sem_ratio  = sd_ratio / sqrt(n),
    .groups = "drop"
  ) 

mean_nanopore_ratio <- ggplot(summary_df, aes(x = groups,
                                              y = mean_ratio,
                                              fill = groups)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_ratio - sem_ratio,
                    ymax = mean_ratio + sem_ratio),
                width = 0.2) +
  geom_jitter(data = group_sample,
              aes(x = groups, y = percentage),
              width = 0.15, size = 2, alpha = 0.6,
              inherit.aes = FALSE) +
  facet_wrap(~ altMidExon_isoform) +
  xlab("") +
  ylab("Isofom ratio") +
  theme_linedraw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, color = "black", face = "bold")
  )

ggsave(filename = "mean_nanopore_ratio.pdf",
       plot = mean_nanopore_ratio, width = 5,height = 4
)
