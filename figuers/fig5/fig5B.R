library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(rtracklayer)
library(data.table)
library(scales)
library(ggplot2)
#### 1. Process reference genome ####
txdb <- txdbmaker::makeTxDbFromGFF("gencode.v38.annotation.gtf")

# Retrieve exons and CDS
exons_gr <- exonsBy(txdb, by = "tx", use.names = TRUE)
cds_gr <- cdsBy(txdb, by = "tx", use.names = TRUE)

# Convert to data.table
exons_dt <- as.data.table(as.data.frame(unlist(exons_gr), row.names = NULL))
exons_dt$transcript_id <- names(unlist(exons_gr))

cds_dt <- as.data.table(as.data.frame(unlist(cds_gr), row.names = NULL))
cds_dt$transcript_id <- names(unlist(cds_gr))

#### 2. Calculate exonic lengths for each region ####
# Get transcript boundaries
transcript_bounds <- exons_dt[, .(
  tx_start = min(start),
  tx_end = max(end),
  strand = unique(as.character(strand)),
  seqnames = unique(as.character(seqnames)),
  total_exon_length = sum(width)
), by = transcript_id]

# Get CDS boundaries and exonic length
cds_bounds <- cds_dt[, .(
  cds_start = min(start),
  cds_end = max(end),
  cds_exon_length = sum(width)
), by = transcript_id]

# Merge
transcript_regions <- merge(transcript_bounds, cds_bounds, 
                            by = "transcript_id", all.x = TRUE)

# Calculate exonic lengths for 5'UTR and 3'UTR
# IMPORTANT: Handle cases where UTR-CDS junction is within a single exon
calculate_region_exon_lengths <- function(exons_dt, cds_dt) {
  result_list <- list()
  tx_ids <- unique(exons_dt$transcript_id)
  total_tx <- length(tx_ids)
  
  for (i in seq_along(tx_ids)) {
    tx_id <- tx_ids[i]
    tx_exons <- exons_dt[transcript_id == tx_id]
    tx_cds <- cds_dt[transcript_id == tx_id]
    
    if (i %% 500 == 0 || i == total_tx) {
      cat(sprintf("Processing %d / %d (%.1f%%)\r", i, total_tx, i / total_tx * 100))
      flush.console()
    }
    
    if (nrow(tx_cds) == 0) next
    
    strand <- unique(tx_exons$strand)
    cds_start <- min(tx_cds$start)
    cds_end <- max(tx_cds$end)
    
    utr5_length <- 0
    utr3_length <- 0
    
    for (j in 1:nrow(tx_exons)) {
      exon_start <- tx_exons$start[j]
      exon_end <- tx_exons$end[j]
      
      if (strand == "+") {
        # 5'UTR
        if (exon_end < cds_start) {
          utr5_length <- utr5_length + (exon_end - exon_start + 1)
        } else if (exon_start < cds_start && exon_end >= cds_start) {
          utr5_length <- utr5_length + (cds_start - exon_start)
        }
        # 3'UTR
        if (exon_start > cds_end) {
          utr3_length <- utr3_length + (exon_end - exon_start + 1)
        } else if (exon_start <= cds_end && exon_end > cds_end) {
          utr3_length <- utr3_length + (exon_end - cds_end)
        }
        
      } else {
        # Negative strand
        if (exon_start > cds_end) {
          utr5_length <- utr5_length + (exon_end - exon_start + 1)
        } else if (exon_start <= cds_end && exon_end > cds_end) {
          utr5_length <- utr5_length + (exon_end - cds_end)
        }
        
        if (exon_end < cds_start) {
          utr3_length <- utr3_length + (exon_end - exon_start + 1)
        } else if (exon_start < cds_start && exon_end >= cds_start) {
          utr3_length <- utr3_length + (cds_start - exon_start)
        }
      }
    }
    
    result_list[[tx_id]] <- data.frame(
      transcript_id = tx_id,
      utr5_exon_length = utr5_length,
      utr3_exon_length = utr3_length
    )
  }
  
  cat("\nDone!\n")
  return(rbindlist(result_list))
}

# Calculate UTR exonic lengths
utr_lengths <- calculate_region_exon_lengths(exons_dt, cds_dt)
saveRDS(utr_lengths,"utr_lengths.RDS")
# Merge with transcript_regions
transcript_regions <- transcript_regions %>%
  left_join(utr_lengths, by = "transcript_id") %>%
  mutate(
    utr5_exon_length = ifelse(is.na(utr5_exon_length), 0, utr5_exon_length),
    utr3_exon_length = ifelse(is.na(utr3_exon_length), 0, utr3_exon_length),
    cds_exon_length = ifelse(is.na(cds_exon_length), 0, cds_exon_length)
  )
saveRDS(transcript_regions,"mRNA_transcript_regions.RDS")

#### 3. Import intron data ####
gtf <- import("gencode.v38.annotation.gtf")
tx_type_df <- as.data.frame(mcols(gtf)[, c("transcript_id", "transcript_type")])
tx_type_df <- tx_type_df %>%
  mutate(transcript_id = str_replace(transcript_id, "\\..*", ""))

protein_coding_transcript <- unique(tx_type_df$transcript_id[tx_type_df$transcript_type == "protein_coding"])

intron_df_full_GRCh38 <- readRDS("intron_df_GC_GRCh38.rds") # derived from short-intron-evolution/data_processing/getIntrons_GRCh38_gtf.R
intron_df_protein_ori <- intron_df_full_GRCh38 %>%
  filter(transcript_id %in% protein_coding_transcript) 

saveRDS(intron_df_protein_ori,"intron_df_protein.RDS")

intron_df_protein <- intron_df_protein_ori %>%
  na.omit()

# For each intron, upstream exon end position
intron_with_exon_boundary <- intron_df_protein %>%
  mutate(
    upstream_exon_end = ifelse(strand == "+", 
                               start - 1, 
                               end + 1)
  )

#### 4. Map introns to regions using upstream exon end ####
intron_with_region <- intron_with_exon_boundary %>%
  left_join(transcript_regions %>%
              mutate(transcript_id = sub("\\..*", "", transcript_id)), 
            by = c("transcript_id", "seqnames", "strand")) %>%
  filter(!is.na(cds_start)) %>%
  mutate(
    region = case_when(
      strand == "+" & upstream_exon_end < cds_start ~ "5UTR",
      strand == "+" & upstream_exon_end >= cds_start & upstream_exon_end <= cds_end ~ "CDS",
      strand == "+" & upstream_exon_end > cds_end ~ "3UTR",
      
      strand == "-" & upstream_exon_end > cds_end ~ "5UTR",
      strand == "-" & upstream_exon_end >= cds_start & upstream_exon_end <= cds_end ~ "CDS",
      strand == "-" & upstream_exon_end < cds_start ~ "3UTR",
      
      TRUE ~ NA_character_
    )
  )

#### 5. Calculate cumulative exonic position for upstream exon end ####
# This function calculates how many EXONIC bases are transcribed
# before reaching the upstream_exon_end position
# IMPORTANT: Properly handles UTR-CDS junctions within single exons
exons_dt$transcript_id <- str_replace(exons_dt$transcript_id, "\\..*", "")
transcript_regions$transcript_id <- str_replace(transcript_regions$transcript_id, "\\..*", "")

# Convert to data.table for speed
setDT(intron_with_region)
setDT(exons_dt)
setDT(transcript_regions)

# Add unique intron ID
intron_with_region[, intron_id := .I]

# Prepare transcript info
transcript_info <- transcript_regions[, .(transcript_id, cds_start, cds_end)]

# Step 1: Join introns with their transcript's exons
intron_exon_dt <- intron_with_region[, .(intron_id, transcript_id, strand, 
                                         upstream_exon_end, region)][
                                           exons_dt[, .(transcript_id, exon_start = start, exon_end = end)],
                                           on = "transcript_id",
                                           allow.cartesian = TRUE
                                         ]

# Join CDS boundaries
intron_exon_dt <- intron_exon_dt[transcript_info, 
                                 on = "transcript_id", 
                                 nomatch = NULL]

# Step 2: Calculate contributions (vectorized)

# Positive strand
intron_exon_dt[strand == "+", `:=`(
  is_before = exon_end < upstream_exon_end,
  contains_up = exon_start <= upstream_exon_end & exon_end >= upstream_exon_end
)]

intron_exon_dt[strand == "+" & is_before == TRUE & region == "5UTR", 
               contribution := pmax(0, pmin(exon_end, cds_start - 1) - exon_start + 1)]

intron_exon_dt[strand == "+" & is_before == TRUE & region == "CDS", 
               contribution := pmax(0, pmin(exon_end, cds_start - 1) - exon_start + 1) +
                 pmax(0, pmin(exon_end, cds_end) - pmax(exon_start, cds_start) + 1)]

intron_exon_dt[strand == "+" & is_before == TRUE & region == "3UTR", 
               contribution := pmax(0, pmin(exon_end, cds_start - 1) - exon_start + 1) +
                 pmax(0, pmin(exon_end, cds_end) - pmax(exon_start, cds_start) + 1) +
                 pmax(0, exon_end - pmax(exon_start, cds_end + 1) + 1)]

intron_exon_dt[strand == "+" & contains_up == TRUE & region == "5UTR",
               contribution := pmax(0, pmin(upstream_exon_end, cds_start - 1) - exon_start + 1)]

intron_exon_dt[strand == "+" & contains_up == TRUE & region == "CDS",
               contribution := fifelse(exon_start < cds_start, cds_start - exon_start, 0) +
                 pmax(0, upstream_exon_end - pmax(exon_start, cds_start) + 1)]

intron_exon_dt[strand == "+" & contains_up == TRUE & region == "3UTR",
               contribution := fifelse(exon_start < cds_start, cds_start - exon_start, 0) +
                 pmax(0, pmin(upstream_exon_end, cds_end) - pmax(exon_start, cds_start) + 1) +
                 fifelse(upstream_exon_end > cds_end, upstream_exon_end - cds_end, 0)]

# Negative strand
intron_exon_dt[strand == "-", `:=`(
  is_before = exon_start > upstream_exon_end,
  contains_up = exon_start <= upstream_exon_end & exon_end >= upstream_exon_end
)]

intron_exon_dt[strand == "-" & is_before == TRUE & region == "5UTR",
               contribution := pmax(0, exon_end - pmax(exon_start, cds_end + 1) + 1)]

intron_exon_dt[strand == "-" & is_before == TRUE & region == "CDS",
               contribution := pmax(0, exon_end - pmax(exon_start, cds_end + 1) + 1) +
                 pmax(0, pmin(exon_end, cds_end) - pmax(exon_start, cds_start) + 1)]

intron_exon_dt[strand == "-" & is_before == TRUE & region == "3UTR",
               contribution := pmax(0, exon_end - pmax(exon_start, cds_end + 1) + 1) +
                 pmax(0, pmin(exon_end, cds_end) - pmax(exon_start, cds_start) + 1) +
                 pmax(0, pmin(exon_end, cds_start - 1) - exon_start + 1)]

intron_exon_dt[strand == "-" & contains_up == TRUE & region == "5UTR",
               contribution := pmax(0, exon_end - pmax(upstream_exon_end, cds_end + 1) + 1)]

intron_exon_dt[strand == "-" & contains_up == TRUE & region == "CDS",
               contribution := fifelse(exon_end > cds_end, exon_end - cds_end, 0) +
                 pmax(0, pmin(exon_end, cds_end) - upstream_exon_end + 1)]

intron_exon_dt[strand == "-" & contains_up == TRUE & region == "3UTR",
               contribution := fifelse(exon_end > cds_end, exon_end - cds_end, 0) +
                 pmax(0, pmin(exon_end, cds_end) - pmax(upstream_exon_end, cds_start) + 1) +
                 fifelse(upstream_exon_end < cds_start, cds_start - exon_start, 0)]

# Replace NA with 0
intron_exon_dt[is.na(contribution), contribution := 0]

# Step 3: Sum contributions by intron
pre_cumulative_pos <- intron_exon_dt[, .(pre_cumulative_exon_pos = sum(contribution, na.rm = TRUE)), 
                                     by = intron_id]


# Join back to original data
intron_with_region <- intron_with_region[pre_cumulative_pos, 
                                         on = "intron_id"]

intron_with_region[, cumulative_exon_pos :=
                     fifelse(
                       region == "5UTR",
                       pre_cumulative_exon_pos,
                       fifelse(
                         region == "CDS",
                         pre_cumulative_exon_pos - utr5_exon_length,
                         fifelse(
                           region == "3UTR",
                           pre_cumulative_exon_pos - utr5_exon_length - cds_exon_length,
                           NA_real_
                         )
                       )
                     )
]
# Convert back to data.frame if needed
intron_with_region <- as.data.frame(intron_with_region)

saveRDS(intron_with_region,file = "intron_with_region.RDS")
#### 6. Normalize positions within intron-less coordinate system ####
intron_normalized <- intron_with_region %>%
  filter(!is.na(region)) %>%
  mutate(
    region_rel_pos = case_when(
      region == "5UTR" ~ cumulative_exon_pos / utr5_exon_length,
      region == "CDS" ~ cumulative_exon_pos / cds_exon_length,
      region == "3UTR" ~ cumulative_exon_pos / utr3_exon_length,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(region_rel_pos >= 0 & region_rel_pos <= 1)

#### 7. Compute scaling factors ####
avg_exon_lengths <- intron_normalized %>%
  group_by(region) %>%
  summarise(
    avg_exon_length = mean(case_when(
      region == "5UTR" ~ utr5_exon_length,
      region == "CDS" ~ cds_exon_length,
      region == "3UTR" ~ utr3_exon_length
    ), na.rm = TRUE)
  )

cds_avg_exon_length <- avg_exon_lengths$avg_exon_length[avg_exon_lengths$region == "CDS"]

scale_factors <- avg_exon_lengths %>%
  mutate(scale_factor = avg_exon_length / cds_avg_exon_length)

intron_normalized <- intron_normalized %>%
  left_join(scale_factors, by = "region") %>%
  mutate(
    global_pos = case_when(
      region == "5UTR" ~ region_rel_pos * scale_factor,
      region == "CDS" ~ scale_factors$scale_factor[scale_factors$region == "5UTR"] + 
        region_rel_pos * scale_factor,
      region == "3UTR" ~ scale_factors$scale_factor[scale_factors$region == "5UTR"] + 
        scale_factors$scale_factor[scale_factors$region == "CDS"] +
        region_rel_pos * scale_factor
    )
  )

#### 8. Plot ####
boundary_5utr_cds <- scale_factors$scale_factor[scale_factors$region == "5UTR"]
boundary_cds_3utr <- boundary_5utr_cds + scale_factors$scale_factor[scale_factors$region == "CDS"]

mRNA_intronGCAT_position <- intron_normalized %>%
  mutate(Introns = ifelse(intron_length>100,
                          "Long",
                          ifelse(
                            intron_GC >= 0.45,"GC-rich\nshort","AT-rich\nshort"
                          )),
         Introns = factor(Introns, levels = c("GC-rich\nshort","AT-rich\nshort", "Long"))) %>%
  ggplot(aes(x = global_pos, fill = Introns, color = Introns)) +
  geom_density(alpha = 0.2, linewidth = 0.5) +
  geom_vline(xintercept = c(boundary_5utr_cds, boundary_cds_3utr), 
             linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(
    breaks = c(boundary_5utr_cds / 2, 
               boundary_5utr_cds + (boundary_cds_3utr - boundary_5utr_cds) / 2,
               boundary_cds_3utr + (max(intron_normalized$global_pos) - boundary_cds_3utr) / 2),
    labels = c("5'UTR", "CDS", "3'UTR")
  ) +
  theme_linedraw(base_size = 14) +
  labs(
    y = "Intron density",
    x = "Normalized position (mRNA)",
    fill = "Introns",
    color = "Introns"
  ) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.87,0.8),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank()
  )
mRNA_intronGCAT_position
ggsave(filename = "mRNA_intronGCAT_position.pdf",plot = mRNA_intronGCAT_position,
       width = 5.4,height = 4,units = "in"
)