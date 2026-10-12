library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(wesanderson)
library(circlize)
library(Matrix)

length_ageScore <- read.csv("short_intron_evolution_heatmap.csv") %>%
  mutate(ageScore = factor(ageScore, 
                           levels = sort(unique(ageScore))
                           )
         )

#### Length ####
heatmap_data_length <- length_ageScore %>%
  select(Human, Mouse, Chicken, Turtle, Frog, Zebrafish) %>%
  as.matrix()

ageScore_colors <- setNames(wes_palette("Darjeeling2", 
                                        length(unique(length_ageScore$ageScore)),
                                        type = "discrete"), 
                            levels(length_ageScore$ageScore))

row_annotation <- rowAnnotation(
  AgeGroup = length_ageScore$ageScore, 
  col = list(AgeGroup = ageScore_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "Age Score")
)

col_fun_length <- colorRamp2(
  c(min(heatmap_data_length), median(heatmap_data_length), max(heatmap_data_length)),  # Map min, median, max
  wes_palette("Zissou1", 3, type = "continuous")  # Use Zissou1 for gradient
)


#### heatmapGC group by intron age ####
heatmap_data_GC <- length_ageScore %>%
  mutate(ageScore = factor(ageScore, levels = sort(unique(ageScore)))) %>%
  select(human_intronGC, mouse_intronGC, chicken_intronGC, 
         turtle_intronGC, frog_intronGC, zebrafish_intronGC) %>% 
  rename("Human" = "human_intronGC",
         "Mouse" = "mouse_intronGC",
         "Chicken" = "chicken_intronGC",
         "Turtle" = "turtle_intronGC",
         "Frog" = "frog_intronGC",
         "Zebrafish" = "zebrafish_intronGC") %>%
  as.matrix()  

col_fun_GC <- colorRamp2(
  c(min(heatmap_data_GC), median(heatmap_data_GC), max(heatmap_data_GC)),  # Map min, median, max
  wes_palette("Zissou1", 3, type = "continuous")  # Use Zissou1 for gradient
)

#### final combined plots ####
shortIntron_Length <- Heatmap(
  heatmap_data_length,
  col = col_fun_length,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 45,
  left_annotation = row_annotation,
  name = "log10(Intron length)",
  column_title = "Intron length"
)

shortIntron_GC <- Heatmap(
  heatmap_data_GC,
  col = col_fun_GC,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 45,
  name = "GC-content",
  column_title = "Intron GC"
)


png("shortIntron_Length_GC_heatmap.png", 
    width = 7, height = 5, units = "in", res = 300)

draw(
  shortIntron_Length + shortIntron_GC,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)

dev.off()
