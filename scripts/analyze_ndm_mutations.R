# ============================================================
# NDM MUTATION ANALYSIS
# ============================================================
# Project:
# Development of an automated pipeline for mutation analysis
# of NDM variants associated with beta-lactam antibiotic resistance
#
# Purpose:
# Reproduce the descriptive mutation-frequency tables used
# for downstream analysis and visualisation.
#
# Authoritative input:
#   analysis/mutations/NDM_mutations.csv
#
# Outputs:
#   overall_mutation_frequency.csv
#   top5_mutations_per_variant.csv
#   variant_counts.csv
#   variant_mutation_frequency.csv
#   variant_mutation_matrix_counts.csv
#   variant_mutation_matrix_percentages.csv
#
# R version used during development:
#   4.3.3
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_root <- normalizePath(".", mustWork = TRUE)

input_file <- file.path(
  project_root,
  "analysis",
  "mutations",
  "NDM_mutations.csv"
)

output_dir <- file.path(
  project_root,
  "audit",
  "r_reconstruction_test"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read authoritative mutation data
# ------------------------------------------------------------

mutations <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "accession",
  "variant",
  "reference_position",
  "reference_residue",
  "observed_residue",
  "mutation",
  "msa_column"
)

missing_columns <- setdiff(required_columns, names(mutations))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. Basic validation
# ------------------------------------------------------------

if (any(is.na(mutations$accession))) {
  stop("Missing accession values detected.")
}

if (any(is.na(mutations$mutation))) {
  stop("Missing mutation values detected.")
}

if (any(is.na(mutations$variant))) {
  stop("Missing variant values detected.")
}

if (any(mutations$variant == "Unknown")) {
  stop("Unknown variant labels detected in authoritative input.")
}

# Number of primary sequences in the final dataset.
# This is deliberately fixed from the curated final dataset.
total_sequences <- 30

# ------------------------------------------------------------
# 4. Overall mutation frequency
# ------------------------------------------------------------

overall_mutation_frequency <- mutations %>%
  count(mutation, name = "sequence_count") %>%
  mutate(
    percentage_of_all_30_sequences =
      round(100 * sequence_count / total_sequences, 2)
  ) %>%
  arrange(desc(sequence_count), mutation)

write_csv(
  overall_mutation_frequency,
  file.path(output_dir, "overall_mutation_frequency.csv")
)

# ------------------------------------------------------------
# 5. Mutation-bearing sequence counts by variant
# ------------------------------------------------------------

variant_counts <- mutations %>%
  distinct(accession, variant) %>%
  count(variant, name = "sequence_count") %>%
  arrange(
    as.numeric(str_extract(variant, "\\d+"))
  )

write_csv(
  variant_counts,
  file.path(output_dir, "variant_counts.csv")
)

# ------------------------------------------------------------
# 6. Variant-specific mutation frequencies
# ------------------------------------------------------------

variant_mutation_frequency <- mutations %>%
  distinct(accession, variant, mutation) %>%
  count(variant, mutation, name = "mutation_count") %>%
  left_join(
    variant_counts,
    by = "variant"
  ) %>%
  mutate(
    percentage_within_variant =
      round(100 * mutation_count / sequence_count, 2)
  ) %>%
  arrange(
    as.numeric(str_extract(variant, "\\d+")),
    desc(mutation_count),
    mutation
  )

write_csv(
  variant_mutation_frequency,
  file.path(output_dir, "variant_mutation_frequency.csv")
)

# ------------------------------------------------------------
# 7. Top five mutations per variant
# ------------------------------------------------------------

top5_mutations_per_variant <- variant_mutation_frequency %>%
  group_by(variant) %>%
  arrange(desc(mutation_count), mutation, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(
    variant,
    mutation,
    mutation_count,
    sequence_count,
    percentage_within_variant
  )

write_csv(
  top5_mutations_per_variant,
  file.path(output_dir, "top5_mutations_per_variant.csv")
)

# ------------------------------------------------------------
# 8. Variant × mutation count matrix
# ------------------------------------------------------------

variant_mutation_matrix_counts <- mutations %>%
  distinct(accession, variant, mutation) %>%
  count(variant, mutation, name = "count") %>%
  pivot_wider(
    names_from = mutation,
    values_from = count,
    values_fill = 0
  ) %>%
  arrange(
    as.numeric(str_extract(variant, "\\d+"))
  )

write_csv(
  variant_mutation_matrix_counts,
  file.path(output_dir, "variant_mutation_matrix_counts.csv")
)

# ------------------------------------------------------------
# 9. Variant × mutation percentage matrix
# ------------------------------------------------------------

variant_mutation_matrix_percentages <- variant_mutation_matrix_counts

mutation_columns <- setdiff(
  names(variant_mutation_matrix_percentages),
  "variant"
)

variant_n <- variant_counts %>%
  select(variant, sequence_count)

variant_mutation_matrix_percentages <- variant_mutation_matrix_percentages %>%
  left_join(
    variant_n,
    by = "variant"
  )

for (col in mutation_columns) {
  variant_mutation_matrix_percentages[[col]] <-
    round(
      100 *
        variant_mutation_matrix_percentages[[col]] /
        variant_mutation_matrix_percentages$sequence_count,
      2
    )
}

variant_mutation_matrix_percentages <- variant_mutation_matrix_percentages %>%
  select(-sequence_count)

write_csv(
  variant_mutation_matrix_percentages,
  file.path(output_dir, "variant_mutation_matrix_percentages.csv")
)

# ------------------------------------------------------------
# 10. Final checks
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("NDM MUTATION R ANALYSIS — TEST RUN\n")
cat("============================================================\n\n")

cat("Input file:\n", input_file, "\n\n")

cat("Total mutation observations:", nrow(mutations), "\n")
cat("Unique accessions:", n_distinct(mutations$accession), "\n")
cat("Unique mutation types:", n_distinct(mutations$mutation), "\n")
cat("Mutation-bearing variants:", n_distinct(mutations$variant), "\n")
cat("Unknown variants:", sum(mutations$variant == "Unknown"), "\n\n")

cat("NDM-16 R264H check:\n")

ndm16_check <- mutations %>%
  filter(
    accession == "A0A0K1Z5K4",
    variant == "NDM-16",
    mutation == "R264H"
  )

if (nrow(ndm16_check) != 1) {
  stop("NDM-16 R264H validation failed.")
}

print(ndm16_check)

cat("\nOutput files created:\n")

output_files <- list.files(
  output_dir,
  pattern = "\\.csv$",
  full.names = FALSE
)

print(sort(output_files))

if (length(output_files) != 6) {
  stop("Expected exactly six R output files.")
}

cat("\nFINAL TEST VERDICT: PASS\n")
cat("============================================================\n")
