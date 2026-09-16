library(dplyr)
library(tidyr)
library(ggplot2)
length_ageScore <- read.csv("/data2/Data/ang_chu_data/short_intron/orthoFinder/ScriptsForOrthoIntrons/short_intron_evolution_heatmap.csv")

### 1.  ageScore == 4
age4_data <- length_ageScore %>% 
  filter(ageScore == 4) %>%
  mutate(oriIntronID_1 = oriIntronID) %>%
  separate(col = oriIntronID_1,
           into = c("humanID","mouseID",'chickenID',"turtleID","frogID","zebrafishID"),sep = ",",
  )
## GC cutoff 
age4_data %>% dplyr::select(humanID,human_intronGC) %>% distinct() %>% pull(human_intronGC) %>% hist()
abline(v=0.45) 
### 2. IQR
age4_data <- age4_data %>%
  rowwise() %>%
  mutate(length_IQR = IQR(c_across(c(Human, Mouse, Chicken, Turtle, Frog, Zebrafish))),
         huamnGC_group = ifelse(human_intronGC<0.45,"AT-rich","GC-rich"))
### ks test 
ks_res <- ks.test(age4_data$length_IQR[age4_data$huamnGC_group == "AT-rich"],
                  age4_data$length_IQR[age4_data$huamnGC_group == "GC-rich"])
ks_stat <- round(ks_res$statistic, 3)
ks_pval <- format.pval(ks_res$p.value, digits = 3, eps = 0.001)

### gg_hitogram
lenIQR_histo <- ggplot(age4_data, aes(x = length_IQR, fill = huamnGC_group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  scale_fill_manual(values = c("AT-rich" = "blue", "GC-rich" = "red")) +
  labs(x = "Cross-species intron length variation (IQR)", 
       y = "Count",
       fill = "Human intron\n(Age = 4)") +
  theme_linedraw(base_size = 18)+
  annotate("text", 
           x = 0.5, y = 2000,   
           label = paste0("Kolmogorov-Smirnov test",
                          "\nD = ", ks_stat,
                          "\nP ", ks_pval),
           hjust = 0, vjust = 1, size = 4)+
  theme(legend.position = c(0.85, 0.5),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3))

lenIQR_histo

ggsave(lenIQR_histo,filename = "lenIQR_histo.pdf",width = 5,height = 6)