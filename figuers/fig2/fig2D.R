library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(ggplot2)

#### import intron PSI ####
real_intron_result_dataframe <- read_csv("hg19_IR_psi.csv")

real_intron_result_dataframe_pivot <- real_intron_result_dataframe %>% 
  tidyr::pivot_longer(c("adipose","adrenal_gland","artery","bladder",
                        "brain","breast","cells","cervix","colon","esophagus",
                        "fallopian_tube","heart","kidney","liver","lung",
                        "minor_salivary_gland","muscle","nerve","ovary","pancreas",
                        "pituitary","prostate","skin","small_intestine","spleen","stomach","testis",
                        "thyroid","uterus","vagina","whole_blood")) %>% na.omit()
ordered_names <- real_intron_result_dataframe_pivot %>%
  filter(width <= 100) %>%
  group_by(name) %>%
  summarize(median_PSI = median(value, na.rm = TRUE)) %>%
  arrange(median_PSI) %>%
  pull(name)

real_intron_result_dataframe_pivot$name <- factor(real_intron_result_dataframe_pivot$name,
                                                  levels = ordered_names)

selected_tissues <- c("adrenal_gland","muscle","brain","nerve","heart","minor_salivary_gland",
                      "lung","testis","liver","kidney","bladder","skin")

my_comparisons <- list(
  c("Short GC-rich intron", "Short AT-rich intron"),
  c("Short AT-rich intron", "Long intron"),
  c("Short GC-rich intron", "Long intron")
)
#### plot fig4D ####
selecetedTissue_shortINT_psi <- real_intron_result_dataframe_pivot %>%
  filter(name %in% selected_tissues) %>%
  mutate(name = str_to_sentence(name)) %>%
  ggplot(aes(x = intron, y = value, fill = intron)) + 
  geom_boxplot(alpha = 0.7, outliers = FALSE) +
  facet_wrap(~ name,ncol=6) +
  stat_compare_means(
    comparisons = my_comparisons,
    label       = "p.signif",
    method      = "wilcox.test",
    label.y     = c(0.9, 0.9, 1.03) 
  ) +
  ylab("PSI (retained intron)") + 
  xlab("") + 
  labs(fill = "") +
  theme_linedraw() +
  coord_cartesian(ylim = c(NA, 1.2)) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, color = "black", face = "bold"))

selecetedTissue_shortINT_psi

ggsave(selecetedTissue_shortINT_psi,
       filename = "selecetedTissue_shortINT_psi.pdf",
       width = 7,height = 6,units = "in")