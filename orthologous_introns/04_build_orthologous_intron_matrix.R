library(dplyr)
library(tidyr)

# ============================================================
# Configuration
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

ROOT <- if (length(args) >= 1) {
  normalizePath(args[1])
} else {
  normalizePath(".")
}

MMSEQS_DIR <- file.path(
  ROOT,
  "MMseqs2_pairwise_output"
)

RDS_DIR <- file.path(
  ROOT,
  "ScriptsForOrthoIntrons",
  "RDS_data"
)

RESULT_DIR <- file.path(
  MMSEQS_DIR,
  "Scripts",
  "ResultRDS"
)

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. Read MMseqs2 RBH results
# ============================================================

load_mmseq_data <- function(dir_path) {
  
  tsv_files <- list.files(
    path = dir_path,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
  
  column_names <- c(
    "Query_ID",
    "Target_ID",
    "Identity",
    "Alignment_Length",
    "Mismatches",
    "Gap_Opens",
    "Query_Start",
    "Query_End",
    "Target_Start",
    "Target_End",
    "E_value",
    "Bit_Score"
  )
  
  result <- lapply(
    tsv_files,
    function(file) {
      
      if (file.info(file)$size == 0) {
        return(NULL)
      }
      
      x <- read.delim(
        file,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
      
      if (nrow(x) == 0) {
        return(NULL)
      }
      
      colnames(x) <- column_names
      
      x
    }
  )
  
  bind_rows(result)
}


mmseq <- list(
  Mouse = load_mmseq_data(
    file.path(MMSEQS_DIR, "Human_Mouse")
  ),
  Chicken = load_mmseq_data(
    file.path(MMSEQS_DIR, "Human_Chicken")
  ),
  Turtle = load_mmseq_data(
    file.path(MMSEQS_DIR, "Human_Turtle")
  ),
  Frog = load_mmseq_data(
    file.path(MMSEQS_DIR, "Human_Frog")
  ),
  Zebrafish = load_mmseq_data(
    file.path(MMSEQS_DIR, "Human_Zebrafish")
  )
)


# ============================================================
# 2. Convert junction IDs to genomic intron IDs
# ============================================================

make_intron_lookup <- function(rds_file) {
  
  readRDS(rds_file) %>%
    mutate(
      intronID = paste0(
        "chr",
        chrom,
        ":",
        intron_start,
        "-",
        intron_end,
        ":",
        ifelse(
          strand == "+",
          "pos",
          "neg"
        )
      ),
      intronLength =
        intron_end - intron_start + 1
    ) %>%
    select(
      junctionID,
      intronID,
      intronLength
    ) %>%
    distinct()
}


human_lookup <- make_intron_lookup(
  file.path(
    RDS_DIR,
    "humanExtractedJunctionProteins.rds"
  )
)

species_lookup <- list(
  Mouse = make_intron_lookup(
    file.path(
      RDS_DIR,
      "mouseExtractedJunctionProteins.rds"
    )
  ),
  Chicken = make_intron_lookup(
    file.path(
      RDS_DIR,
      "chickenExtractedJunctionProteins.rds"
    )
  ),
  Turtle = make_intron_lookup(
    file.path(
      RDS_DIR,
      "turtleExtractedJunctionProteins.rds"
    )
  ),
  Frog = make_intron_lookup(
    file.path(
      RDS_DIR,
      "frogExtractedJunctionProteins.rds"
    )
  ),
  Zebrafish = make_intron_lookup(
    file.path(
      RDS_DIR,
      "zebrafishExtractedJunctionProteins.rds"
    )
  )
)


map_rbh_to_introns <- function(
    mmseq_data,
    target_lookup,
    species_name
) {
  
  human_table <- human_lookup %>%
    rename(
      Query_ID = junctionID,
      intronID_Human = intronID,
      intronLength_Human = intronLength
    )
  
  target_table <- target_lookup %>%
    rename(
      Target_ID = junctionID
    )
  
  names(target_table)[
    names(target_table) == "intronID"
  ] <- paste0(
    "intronID_",
    species_name
  )
  
  names(target_table)[
    names(target_table) == "intronLength"
  ] <- paste0(
    "intronLength_",
    species_name
  )
  
  mmseq_data %>%
    select(
      Query_ID,
      Target_ID
    ) %>%
    distinct() %>%
    left_join(
      human_table,
      by = "Query_ID"
    ) %>%
    left_join(
      target_table,
      by = "Target_ID"
    ) %>%
    select(
      -Query_ID,
      -Target_ID
    ) %>%
    drop_na() %>%
    distinct()
}


intron_matches <- lapply(
  names(mmseq),
  function(species_name) {
    
    map_rbh_to_introns(
      mmseq_data =
        mmseq[[species_name]],
      target_lookup =
        species_lookup[[species_name]],
      species_name =
        species_name
    )
  }
)

names(intron_matches) <- names(mmseq)


# ============================================================
# 3. Construct human-centered orthologous intron matrix
# ============================================================

allHumanOrthoCDSIntrons <- human_lookup %>%
  select(
    intronID,
    intronLength
  ) %>%
  distinct() %>%
  rename(
    intronID_Human = intronID,
    intronLength_Human = intronLength
  )


allJoined <- allHumanOrthoCDSIntrons

for (species_name in names(intron_matches)) {
  
  allJoined <- full_join(
    allJoined,
    intron_matches[[species_name]],
    by = c(
      "intronID_Human",
      "intronLength_Human"
    )
  )
}


saveRDS(
  allJoined,
  file.path(
    RESULT_DIR,
    "allJoined.rds"
  )
)


# ============================================================
# 4. Filter human introns
# ============================================================

valid_hg38 <-
  "^chr[0-9XYMT]+:[0-9]+-[0-9]+:(pos|neg)$"

filtered_allJoined <- allJoined %>%
  filter(
    grepl(
      valid_hg38,
      intronID_Human
    ),
    intronLength_Human >= 10
  )


# ============================================================
# 5. Prepare cross-species intron-length matrix
# ============================================================

species_names <- c(
  "Mouse",
  "Chicken",
  "Turtle",
  "Frog",
  "Zebrafish"
)

filtered_allJoined_reframe <-
  filtered_allJoined


# Replace missing intron IDs
for (species_name in species_names) {
  
  id_col <- paste0(
    "intronID_",
    species_name
  )
  
  filtered_allJoined_reframe[[id_col]][
    is.na(
      filtered_allJoined_reframe[[id_col]]
    )
  ] <- paste0(
    "NA_intronID_",
    species_name
  )
}


# Combined orthologous intron identifier
filtered_allJoined_reframe <-
  filtered_allJoined_reframe %>%
  mutate(
    intronID_ALL = paste0(
      intronID_Human,
      ",",
      intronID_Mouse,
      ",",
      intronID_Chicken,
      ",",
      intronID_Turtle,
      ",",
      intronID_Frog,
      ",",
      intronID_Zebrafish,
      ","
    )
  )


# Remove individual intron ID columns
filtered_allJoined_reframe <-
  filtered_allJoined_reframe %>%
  select(
    -starts_with("intronID_"),
    intronID_ALL
  )


# Missing orthologous introns are represented by zero
# after log10 transformation:
# NA -> 1 -> log10(1) = 0
filtered_allJoined_reframe <-
  filtered_allJoined_reframe %>%
  mutate(
    across(
      starts_with("intronLength"),
      ~ replace_na(.x, 1)
    ),
    across(
      starts_with("intronLength"),
      log10
    )
  )


filtered_allJoined_heatMap <-
  filtered_allJoined_reframe %>%
  select(
    -intronID_ALL
  )


rownames(filtered_allJoined_heatMap) <-
  filtered_allJoined_reframe$intronID_ALL


colnames(filtered_allJoined_heatMap) <-
  c(
    "Human",
    "Mouse",
    "Chicken",
    "Turtle",
    "Frog",
    "Zebrafish"
  )


saveRDS(
  filtered_allJoined_heatMap,
  file.path(
    RESULT_DIR,
    "filtered_allJoined_heatMap.rds"
  )
)