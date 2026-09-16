library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)

human_short_counts <- read.csv("human_uniqShortIntron_counts.csv")
ageScore_colors <- setNames(wes_palette("Darjeeling2", length(unique(human_short_counts$ageScore)),
                                        type = "discrete"), 
                            levels(human_short_counts$ageScore))

human_ortho_short_intron_ratio <- ggplot(human_short_counts,
                                         aes(x = "Human short\northologous introns", y = percent, fill = factor(ageScore))) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  scale_fill_manual(values = ageScore_colors) +
  coord_flip() +
  labs(x = "",
       y = "",
       fill = "Age Score") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 17),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.5, 0),             
    legend.justification = c("center", "bottom"),
    legend.direction = "horizontal"         
  )

human_ortho_short_intron_ratio

ggsave(human_ortho_short_intron_ratio,filename = "human_ortho_short_intron_ratio.pdf",
       width=8.5,height = 1.7)