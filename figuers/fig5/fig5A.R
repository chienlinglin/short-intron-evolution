library(dplyr)
library(tidyr)
library(Biostrings)
library(ggseqlogo)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

unique_intronPOS <- read.csv("unique_intronPOS.csv") # derived from short-intron-evolution/data_processing/human_introns_examined.R

## intronic - 5'SS:15 nt & 3'SS:20
## exonic - 10 nt
unique_SSregions <- unique_intronPOS %>%
  mutate(region5_start = ifelse(strand=="+", start-10, end-14),
         region5_end = ifelse(strand=="+", start+14, end+10),
         region3_start = ifelse(strand=="+", end-19, start-10),
         region3_end = ifelse(strand=="+", end+10, start+19))

## get sequences
unique_SSregions$region5_seq <- 
  as.character(getSeq(genome, unique_SSregions$seqnames,
                      start=unique_SSregions$region5_start, 
                      end=unique_SSregions$region5_end,
                      strand=unique_SSregions$strand))

unique_SSregions$region3_seq <- 
  as.character(getSeq(genome, unique_SSregions$seqnames,
                      start=unique_SSregions$region3_start, 
                      end=unique_SSregions$region3_end,
                      strand=unique_SSregions$strand))


## ggseqlogo
# all introns
all_5 <- ggseqlogo(unique_SSregions$region5_seq)+ ggtitle("all intron 5' region")
all_3 <- ggseqlogo(unique_SSregions$region3_seq,method = "prob")+ ggtitle("all intron 3' region")

## function for 5ss
fun_5 <- function(df,titles){
  df %>%
    pull(region5_seq) %>%
    ggseqlogo(method = 'prob' ) +
    annotate('rect', xmin = 10.5, xmax = 12.5, ymin = -0.05, ymax = 1.05, 
             alpha = .1, col='black', fill='burlywood2') + 
    annotate('text', x=11.5, y=1.1, label="5'SS") +  
    xlab(titles) +
    coord_cartesian(xlim = c(1,30)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) 
}
## function for 3ss
fun_3 <- function(df,titles){
  df %>%
    pull(region3_seq) %>%
    ggseqlogo(method = 'prob' )+
    annotate('rect', xmin = 18.5, xmax = 20.5, ymin = -0.05, ymax = 1.05, 
             alpha = .1, col='black', fill='burlywood2') + 
    annotate('text', x=19.5, y=1.1, label="3'SS") +
    xlab(titles)+ 
    coord_cartesian(xlim = c(1,30)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}
# short introns
short_only <- unique_SSregions %>% filter(intron_length<=100) %>%
  mutate(GC_type = ifelse(
    GC_content >= 0.45,"GC-rich","AT-rich"
  ))
shortGC_5 <- fun_5(short_only %>%
                     filter(GC_type=="GC-rich"),"GC-rich short introns")
shortAT_5 <- fun_5(short_only %>%
                     filter(GC_type=="AT-rich"),"AT-rich short introns")
shortGC_3 <- fun_3(short_only %>%
                     filter(GC_type=="GC-rich"),"GC-rich short introns")
shortAT_3 <- fun_3(short_only %>%
                     filter(GC_type=="AT-rich"),"AT-rich short introns")

# long introns
long_only <- unique_SSregions %>% filter(intron_length>100)
long_5 <- fun_5(long_only,"Long introns")
long_3 <- fun_3(long_only,"Long introns")

# plot fig5A
pdf("introns_ggseq.pdf", width = 10, height = 5.5)

gridExtra::grid.arrange(shortGC_5, shortGC_3,
                        shortAT_5, shortAT_3,
                        long_5,long_3, ncol = 2)
dev.off()