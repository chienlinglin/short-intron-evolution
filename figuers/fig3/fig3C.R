library(dplyr)
library(ggrepel)
library(ggplot2)
RBP_clip_result <- read.csv("RBP_clip_result.csv")

#### plot fig5D ####
selected_RBPs <- c("SFPQ","PUS1","METTL14","METTL3","ALKBH5",
                   "FUS","ELAVL1", "U2AF2","FUBP3","HNRNPC",
                   "SRSF1","HNRNPH","HNRNPL","HNRNPH1","AQR"
)
res_plot <- res %>%
  as.data.frame() %>%
  mutate(
    RBP = as.character(RBP),
    log2OR = log2(odds_ratio),
    log10FDR = -log10(padj),
    RBP_label = case_when(preference == "short" ~ RBP,
                          RBP %in% selected_RBPs ~ RBP,
                          TRUE ~ NA
    )
  )

RBP_enrichment <- res_plot %>%
  ggplot(
    aes(x = log2OR,
        y = log10FDR,
        color = preference,
        label = RBP_label)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(
    name = "",   # legend title
    values = c(
      short = "#D7191C",
      none  = "grey70",
      long  = "#2C7BB6"
    ),
    labels = c(
      short = "Short intron\nenriched",
      none  = "No preference",
      long  = "Long intron\nenriched"
    )
  ) +
  xlab(expression(log[2]*"(Odds Ratio)")) +
  ylab(expression(-log[10]*"(FDR)"))+
  geom_text_repel(na.rm = TRUE, max.overlaps = 25,
                  box.padding = 1,  
                  size = 3,
                  point.padding = 0.1,
                  segment.color = "grey50",
                  color = "black"
  ) +
  theme_linedraw(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.13, 0.85),
        legend.text = element_text(face = "bold",size=10),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        strip.background = element_blank(),
        legend.key = element_blank()) 

RBP_enrichment

ggsave(RBP_enrichment,
       filename = "RBP_enrichment.pdf",
       width = 5.5,height = 4.4)
