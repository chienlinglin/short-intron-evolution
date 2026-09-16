library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)
library(biomaRt)

# ============================================================
# Configuration
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

ROOT <- if (length(args) >= 1) {
  normalizePath(args[1])
} else {
  normalizePath(".")
}

ORTHOFINDER_DIR <- file.path(
  ROOT,
  "ensemblProteinSequence",
  "OrthoFinder",
  "Results_Nov21_1"
)

PROTEIN_DIR <- file.path(
  ROOT,
  "ensemblProteinSequence"
)

GTF_DIR <- file.path(
  ROOT,
  "ensemblGTFfiles"
)

RDS_DIR <- file.path(
  ROOT,
  "ScriptsForOrthoIntrons",
  "RDS_data"
)

FUNCTION_DIR <- file.path(
  ROOT,
  "function"
)

PAIRED_ID_DIR <- file.path(
  ROOT,
  "extractedJunctionPeptides",
  "PairedIDTable"
)

dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PAIRED_ID_DIR, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. Import OrthoFinder results
# ============================================================

hg38 <- read.delim(
  file.path(
    ORTHOFINDER_DIR,
    "Orthologues",
    "Homo_sapiens.GRCh38.pep.all.tsv"
  ),
  header = TRUE,
  sep = "\t"
)

# Retain ortholog groups represented in all six species
hg38_expanded <- hg38 %>%
  separate_rows(
    Homo_sapiens.GRCh38.pep.all,
    sep = ", "
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = Orthologs
  ) %>%
  na.omit() %>%
  mutate(
    Homo_sapiens.GRCh38.pep.all =
      gsub("\\..*", "", Homo_sapiens.GRCh38.pep.all)
  )


# ============================================================
# 2. Convert protein IDs to transcript IDs
# ============================================================

hg38_ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  version = 113
)

human_protein_ids <- unique(
  hg38_expanded$Homo_sapiens.GRCh38.pep.all
)

human_annotations <- getBM(
  attributes = c(
    "ensembl_peptide_id",
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "external_gene_name"
  ),
  filters = "ensembl_peptide_id",
  values = human_protein_ids,
  mart = hg38_ensembl
)

human_annotations$external_gene_name[
  human_annotations$external_gene_name == ""
] <- NA


hg38_expanded_with_transcripts <- hg38_expanded %>%
  left_join(
    human_annotations,
    by = c(
      "Homo_sapiens.GRCh38.pep.all" =
        "ensembl_peptide_id"
    )
  ) %>%
  rename(
    Human = Homo_sapiens.GRCh38.pep.all,
    Zebrafish = Danio_rerio.GRCz11.pep.all,
    `Western painted turtle` =
      Chrysemys_picta_bellii.Chrysemys_picta_bellii.3.0.3.pep.all,
    Mouse = Mus_musculus.GRCm39.pep.all,
    `Western clawed frog` =
      Xenopus_tropicalis.UCB_Xtro_10.0.pep.all,
    Chicken =
      Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all,
    hg38_transcript_id = ensembl_transcript_id,
    hg38_ensembl_gene_id = ensembl_gene_id,
    hg38_gene_name = external_gene_name
  ) %>%
  select(
    Orthogroup,
    Human,
    hg38_transcript_id,
    hg38_ensembl_gene_id,
    hg38_gene_name,
    Mouse,
    Chicken,
    `Western painted turtle`,
    `Western clawed frog`,
    Zebrafish
  )

saveRDS(
  hg38_expanded_with_transcripts,
  file.path(
    RDS_DIR,
    "hg38_expanded_with_transcripts.rds"
  )
)


# ============================================================
# 3. Annotate orthologous proteins from other species
# ============================================================

species <- c(
  "Zebrafish",
  "Western painted turtle",
  "Mouse",
  "Western clawed frog",
  "Chicken"
)

species_list_cleaned <- lapply(
  species,
  function(species_name) {
    
    hg38_expanded_with_transcripts %>%
      select(
        Human,
        hg38_transcript_id,
        hg38_ensembl_gene_id,
        hg38_gene_name,
        all_of(species_name)
      ) %>%
      rename(
        Species_Data = all_of(species_name)
      ) %>%
      separate_rows(
        Species_Data,
        sep = ", "
      ) %>%
      mutate(
        Species_Data =
          gsub("\\..*", "", Species_Data)
      ) %>%
      distinct()
  }
)

names(species_list_cleaned) <- species


# Ensembl connections
# Version settings below follow the original analysis.

mm39_ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl"
)

chicken_ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "ggallus_gene_ensembl"
)

turtle_ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "cpbellii_gene_ensembl"
)

xenTro10_ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "xtropicalis_gene_ensembl",
  version = 111
)

zebrafish_ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "drerio_gene_ensembl"
)

ensembl_connections <- list(
  "Mouse" = mm39_ensembl,
  "Chicken" = chicken_ensembl,
  "Western painted turtle" = turtle_ensembl,
  "Western clawed frog" = xenTro10_ensembl,
  "Zebrafish" = zebrafish_ensembl
)


get_annotations <- function(protein_ids, ensembl) {
  
  if (length(protein_ids) == 0) {
    return(data.frame())
  }
  
  getBM(
    attributes = c(
      "ensembl_peptide_id",
      "ensembl_transcript_id",
      "ensembl_gene_id",
      "external_gene_name"
    ),
    filters = "ensembl_peptide_id",
    values = protein_ids,
    mart = ensembl
  )
}


annotated_species_list <- lapply(
  names(species_list_cleaned),
  function(species_name) {
    
    species_data <-
      species_list_cleaned[[species_name]]
    
    protein_ids <-
      unique(species_data$Species_Data)
    
    annotations <- get_annotations(
      protein_ids,
      ensembl_connections[[species_name]]
    )
    
    species_data %>%
      left_join(
        annotations,
        by = c(
          "Species_Data" = "ensembl_peptide_id"
        )
      )
  }
)

names(annotated_species_list) <-
  names(species_list_cleaned)

saveRDS(
  annotated_species_list,
  file.path(
    RDS_DIR,
    "annotated_species_list.rds"
  )
)


# ============================================================
# 4. Retrieve CDS intron coordinates from GTF
# ============================================================

source(
  file.path(
    FUNCTION_DIR,
    "process_gtf_introns_with_cds.R"
  )
)


# Human
human_gtf <- process_gtf_introns_with_cds(
  file.path(
    GTF_DIR,
    "Homo_sapiens.GRCh38.113.chr_patch_hapl_scaff.gtf"
  )
)

hg38_OrthoTransCDS <-
  hg38_expanded_with_transcripts %>%
  select(
    Human,
    hg38_transcript_id,
    hg38_ensembl_gene_id,
    hg38_gene_name
  ) %>%
  merge(
    human_gtf,
    by.x = "hg38_transcript_id",
    by.y = "transcript_id",
    all.x = TRUE
  )

saveRDS(
  hg38_OrthoTransCDS,
  file.path(RDS_DIR, "hg38_OrthoTransCDS.rds")
)


# Helper for non-human species
build_species_cds <- function(
    species_name,
    gtf_file
) {
  
  gtf_data <- process_gtf_introns_with_cds(
    file.path(GTF_DIR, gtf_file)
  )
  
  annotated_species_list[[species_name]] %>%
    select(
      Human,
      hg38_transcript_id,
      Species_Data,
      ensembl_transcript_id
    ) %>%
    merge(
      gtf_data,
      by.x = "ensembl_transcript_id",
      by.y = "transcript_id",
      all.x = TRUE
    )
}


mm39_OrthoTransCDS <- build_species_cds(
  "Mouse",
  "Mus_musculus.GRCm39.113.gtf"
)

galGal6_OrthoTransCDS <- build_species_cds(
  "Chicken",
  "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.113.gtf"
)

turtle_OrthoTransCDS <- build_species_cds(
  "Western painted turtle",
  "Chrysemys_picta_bellii.Chrysemys_picta_bellii-3.0.3.113.gtf"
)

xenTro10_OrthoTransCDS <- build_species_cds(
  "Western clawed frog",
  "Xenopus_tropicalis.UCB_Xtro_10.0.113.gtf"
)

danRer11_OrthoTransCDS <- build_species_cds(
  "Zebrafish",
  "Danio_rerio.GRCz11.113.chr_patch_hapl_scaff.gtf"
)


saveRDS(
  mm39_OrthoTransCDS,
  file.path(RDS_DIR, "mm39_OrthoTransCDS.rds")
)

saveRDS(
  galGal6_OrthoTransCDS,
  file.path(RDS_DIR, "galGal6_OrthoTransCDS.rds")
)

saveRDS(
  turtle_OrthoTransCDS,
  file.path(RDS_DIR, "turtle_OrthoTransCDS.rds")
)

saveRDS(
  xenTro10_OrthoTransCDS,
  file.path(RDS_DIR, "xenTro10_OrthoTransCDS.rds")
)

saveRDS(
  danRer11_OrthoTransCDS,
  file.path(RDS_DIR, "danRer11_OrthoTransCDS.rds")
)


# ============================================================
# 5. Extract exon-junction peptides
# ============================================================

source(
  file.path(
    FUNCTION_DIR,
    "JunctionProteinExtraction.R"
  )
)


humanExtractedJunctionProteins <-
  JunctionProteinExtraction(
    protein_path = file.path(
      PROTEIN_DIR,
      "Homo_sapiens.GRCh38.pep.all.fa"
    ),
    Species_OrthoTransCDS = hg38_OrthoTransCDS,
    SpeciesProteinID = "Human",
    SpeciesTranscriptID = "hg38_transcript_id"
  )


mouseExtractedJunctionProteins <-
  JunctionProteinExtraction(
    protein_path = file.path(
      PROTEIN_DIR,
      "Mus_musculus.GRCm39.pep.all.fa"
    ),
    Species_OrthoTransCDS = mm39_OrthoTransCDS,
    SpeciesProteinID = "Species_Data",
    SpeciesTranscriptID = "ensembl_transcript_id"
  ) %>%
  distinct()


chickenExtractedJunctionProteins <-
  JunctionProteinExtraction(
    protein_path = file.path(
      PROTEIN_DIR,
      "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa"
    ),
    Species_OrthoTransCDS = galGal6_OrthoTransCDS,
    SpeciesProteinID = "Species_Data",
    SpeciesTranscriptID = "ensembl_transcript_id"
  ) %>%
  distinct()


turtleExtractedJunctionProteins <-
  JunctionProteinExtraction(
    protein_path = file.path(
      PROTEIN_DIR,
      "Chrysemys_picta_bellii.Chrysemys_picta_bellii-3.0.3.pep.all.fa"
    ),
    Species_OrthoTransCDS = turtle_OrthoTransCDS,
    SpeciesProteinID = "Species_Data",
    SpeciesTranscriptID = "ensembl_transcript_id"
  ) %>%
  distinct()


frogExtractedJunctionProteins <-
  JunctionProteinExtraction(
    protein_path = file.path(
      PROTEIN_DIR,
      "Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa"
    ),
    Species_OrthoTransCDS = xenTro10_OrthoTransCDS,
    SpeciesProteinID = "Species_Data",
    SpeciesTranscriptID = "ensembl_transcript_id"
  ) %>%
  distinct()


zebrafishExtractedJunctionProteins <-
  JunctionProteinExtraction(
    protein_path = file.path(
      PROTEIN_DIR,
      "Danio_rerio.GRCz11.pep.all.fa"
    ),
    Species_OrthoTransCDS = danRer11_OrthoTransCDS,
    SpeciesProteinID = "Species_Data",
    SpeciesTranscriptID = "ensembl_transcript_id"
  ) %>%
  distinct()


junction_objects <- list(
  human = humanExtractedJunctionProteins,
  mouse = mouseExtractedJunctionProteins,
  chicken = chickenExtractedJunctionProteins,
  turtle = turtleExtractedJunctionProteins,
  frog = frogExtractedJunctionProteins,
  zebrafish = zebrafishExtractedJunctionProteins
)

for (nm in names(junction_objects)) {
  saveRDS(
    junction_objects[[nm]],
    file.path(
      RDS_DIR,
      paste0(nm, "ExtractedJunctionProteins.rds")
    )
  )
}


# ============================================================
# 6. Export human-species protein pairing tables
# ============================================================

for (species_name in names(annotated_species_list)) {
  
  sub_data <-
    annotated_species_list[[species_name]]
  
  file_name <- paste0(
    gsub(" ", "_", species_name),
    "_vs_human.csv"
  )
  
  write.csv(
    sub_data,
    file = file.path(
      PAIRED_ID_DIR,
      file_name
    ),
    row.names = FALSE
  )
}