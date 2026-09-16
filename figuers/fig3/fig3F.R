library(dplyr)
library(gridExtra)
library(tibble)
library(tidytext)
library(forcats)
library(patchwork)

go_top15_all_group <- read.csv("go_top15_all_group.csv")
all_group <- go_top15_all_group %>%
  mutate(Description_unique = paste(type, Description, sep = "___")) %>%
  group_by(type) %>%
  mutate(
    Description_unique = fct_reorder(Description_unique, FoldEnrichment)
  ) %>%
  ungroup()


final_plot_older <- all_group  %>%
  filter(!type %in% c("Intron age = 0 or 1", "Intron age = 2 or 3")) %>%
  mutate(type = factor(type,
                       levels = c("Intron age = 4 and AT-rich",
                                  "Intron age = 4 mixed AT/GC",
                                  "Intron age = 4 and GC-rich"))
  ) %>%
  ggplot(aes(x = `Description_unique`, y = `FoldEnrichment`)) + 
  geom_point(stat = 'identity', 
             aes(color = `p.adjust`, 
                 shape = `ONTOLOGY`,
                 size = `Count`)) +
  scale_color_gradient(low = "red", high = "black", limits = c(0, 0.2),
                       guide = guide_colorbar(
                         title.position = "top",
                         label.theme = element_text(size = 8)   
                       )) + 
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_size_continuous(
    range = c(1, 6),
    limits = c(min(all_group$Count), max(all_group$Count)),
    breaks = unique(c(
      min(all_group$Count),
      pretty(all_group$Count, n = 4),
      max(all_group$Count)))
  )+
  scale_shape_manual(values = c("BP" = 16, "CC" = 17, "MF" = 18),drop = FALSE) + 
  facet_wrap( ~ type, scales="free_y", ncol = 3) +
  coord_flip() +
  guides(shape = guide_legend(override.aes = list(size = 9))) +
  theme_linedraw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, color = "black", face = "bold")) +
  labs(x="")

final_plot_older

ggsave(filename = "goTermOrthoShortIntron_ageScore_younger_group_2.pdf",
       plot = final_plot_younger, width = 14, height = 6)

ggsave(filename = "goTermOrthoShortIntron_ageScore_older_group_2.pdf",
       plot = final_plot_older, width = 19, height = 6.5)