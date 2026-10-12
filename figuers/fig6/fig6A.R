library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
# sample 1 = WT
# sample 2 = SON KD/ SRRM2 KD
# IncLevelDifference: average(IncLevel1) - average(IncLevel2)
### SON
son_summary <- read.delim("GSE242898/rMATs/rMATsOutput_WTvsSON/ASE_son_summary.txt")

son_sig_only_long <- son_summary %>%
  transmute(
    EventType,
    `WT Higher Inclusion` = SigEventsJCSample1HigherInclusion,
    `KD Higher Inclusion` = SigEventsJCSample2HigherInclusion
  ) %>%
  pivot_longer(cols = -EventType, names_to = "Direction", values_to = "Count") %>%
  filter(Count > 0)
pdf("GSE242898/rMATs/rScripts/plot/son_summary.pdf",width = 6,height = 7)
ggplot(son_sig_only_long, aes(x = EventType, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = comma(Count)), position = position_stack(vjust = 0.5), size = 6, color = "black") +
  scale_fill_manual(values = c(
    "WT Higher Inclusion" = "#86A8C6",  
    "KD Higher Inclusion" = "#D97C7C"   
  )) +
  labs(
    x = "Alternative Splicing Event Type (KD SON)",
    y = "Number of Significant Events",
    fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.25, 0.88))
dev.off()