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

INPUT_FILE <- file.path(
  ROOT,
  "MMseqs2_pairwise_output",
  "Scripts",
  "ResultRDS",
  "filtered_allJoined_heatMap.rds"
)

OUTPUT_FILE <- file.path(
  ROOT,
  "ScriptsForOrthoIntrons",
  "RDS_data",
  "intronBySpecies_2_ageScore.rds"
)


# ============================================================
# 1. Load orthologous intron matrix
# ============================================================

filtered_allJoined_heatMap <- readRDS(INPUT_FILE)

intronBySpecies <- data.frame(
  oriIntronID = rownames(filtered_allJoined_heatMap)
)


# ============================================================
# 2. Parse species-specific intron IDs
# ============================================================

intronBySpecies <- intronBySpecies %>%
  separate(
    oriIntronID,
    into = c(
      "human",
      "mouse",
      "chicken",
      "turtle",
      "frog",
      "zebrafish"
    ),
    sep = ",",
    extra = "drop",
    remove = FALSE
  ) %>%
  mutate(
    across(
      human:zebrafish,
      ~ ifelse(
        grepl("^NA_intronID_", .x),
        NA,
        .x
      )
    )
  )


# ============================================================
# 3. Assign intron evolutionary age
# ============================================================

intronBySpecies_age <- intronBySpecies %>%
  mutate(
    ageScore = case_when(
      
      # Conserved to zebrafish
      !is.na(zebrafish) ~ 4,
      
      # Conserved to frog, but not zebrafish
      is.na(zebrafish) &
        !is.na(frog) ~ 3,
      
      # Conserved to an amniote species
      is.na(zebrafish) &
        is.na(frog) &
        (!is.na(turtle) | !is.na(chicken)) ~ 2,
      
      # Conserved to mouse only among non-human species
      is.na(zebrafish) &
        is.na(frog) &
        is.na(turtle) &
        is.na(chicken) &
        !is.na(mouse) ~ 1,
      
      # Human-specific
      is.na(zebrafish) &
        is.na(frog) &
        is.na(turtle) &
        is.na(chicken) &
        is.na(mouse) &
        !is.na(human) ~ 0,
      
      TRUE ~ NA_real_
    )
  )


# ============================================================
# 4. Export intron age scores
# ============================================================

intronBySpecies_2_ageScore <- intronBySpecies_age %>%
  select(
    oriIntronID,
    ageScore
  ) %>%
  distinct()


saveRDS(
  intronBySpecies_2_ageScore,
  OUTPUT_FILE
)