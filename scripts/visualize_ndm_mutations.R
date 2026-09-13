#!/usr/bin/env Rscript

# ============================================================
# NDM MUTATION VISUALISATION PIPELINE
# ============================================================
#
# Purpose:
#   Generate reproducible, publication-quality visualisations
#   from the validated NDM mutation dataset.
#
# Authoritative mutation input:
#   analysis/mutations/NDM_mutations.csv
#
# Primary sequence metadata:
#   analysis/reference/primary_sequences.tsv
#
# MSA:
#   analysis/msa/NDM_MSA.fasta
#
# Output:
#   audit/visualisation_test/
#
# Important:
#   This script does NOT regenerate mutation calls.
#   It uses the validated mutation CSV as the authoritative
#   mutation dataset.
#
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(viridis)
  library(patchwork)
})

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("--file=", args_all)]

if (length(script_arg) == 0) {
  stop("FAIL: Could not determine the visualization script path.")
}

script_path <- sub("^--file=", "", script_arg[1])

BASE <- normalizePath(
  file.path(
    dirname(script_path),
    ".."
  ),
  mustWork = FALSE
)

MUTATION_FILE <- file.path(
  BASE, "analysis", "mutations", "NDM_mutations.csv"
)

PRIMARY_FILE <- file.path(
  BASE, "analysis", "reference", "primary_sequences.tsv"
)

MSA_FILE <- file.path(BASE, "analysis", "msa", "NDM_MSA.fasta")

OUT_DIR <- file.path(BASE, "results", "figures")

DIRS <- character(0)



# ------------------------------------------------------------
# BASIC INPUT VALIDATION
# ------------------------------------------------------------

required_mutation_columns <- c(
  "accession",
  "variant",
  "reference_position",
  "reference_residue",
  "observed_residue",
  "mutation",
  "msa_column"
)

mutations <- read_csv(
  MUTATION_FILE,
  show_col_types = FALSE
)

primary <- read_tsv(
  PRIMARY_FILE,
  col_names = FALSE,
  show_col_types = FALSE
)

# The primary metadata file stores the FASTA-style header in the
# first column rather than as a named accession column.
if (ncol(primary) < 1) {
  stop("FAIL: Primary sequence metadata is empty.")
}

names(primary)[1] <- "header"

primary <- primary |>
  mutate(
    accession = case_when(
      str_detect(header, "^(tr|sp)\\|") ~ str_split_fixed(header, "\\|", 3)[,2],
      TRUE ~ word(header, 1)
    ),
    variant = str_extract(
      header,
      regex("NDM-\\d+", ignore_case = TRUE)
    ) |>
      str_to_upper()
  )

if (any(is.na(primary$accession))) {
  stop("FAIL: Could not extract accession from primary sequence metadata.")
}

if (any(is.na(primary$variant))) {
  stop("FAIL: Could not extract NDM variant from primary sequence metadata.")
}

if (!all(required_mutation_columns %in% names(mutations))) {
  stop("FAIL: Required mutation columns are missing.")
}

if (nrow(mutations) != 106) {
  stop(
    paste(
      "FAIL: Expected 106 mutation observations;",
      "found", nrow(mutations)
    )
  )
}

if (n_distinct(mutations$accession) != 27) {
  stop("FAIL: Expected 27 mutation-bearing accessions.")
}

if (n_distinct(mutations$mutation) != 82) {
  stop("FAIL: Expected 82 unique mutation types.")
}

if (any(mutations$variant == "Unknown")) {
  stop("FAIL: Unknown variant labels detected.")
}

if (!"accession" %in% names(primary)) {
  stop("FAIL: Primary sequence metadata lacks accession column.")
}

# ------------------------------------------------------------
# NDM-16 CRITICAL CHECK
# ------------------------------------------------------------

ndm16 <- mutations |>
  filter(accession == "A0A0K1Z5K4")

if (
  nrow(ndm16) != 1 ||
  ndm16$variant != "NDM-16" ||
  ndm16$mutation != "R264H"
) {
  stop("FAIL: NDM-16/R264H validation check failed.")
}

# ------------------------------------------------------------
# DATASET COMPOSITION
# ------------------------------------------------------------

extract_variant <- function(x) {
  m <- str_extract(x, regex("NDM-\\d+", ignore_case = TRUE))
  str_to_upper(m)
}

# The primary metadata file contains the accession and sequence
# information. Variant names are extracted from the sequence
# headers when required.
if (!"variant" %in% names(primary)) {

  header_column <- names(primary)[
    which(
      sapply(
        primary,
        function(x) any(str_detect(as.character(x), "NDM-\\d+"))
      )
    )[1]
  ]

  if (is.na(header_column)) {
    stop("FAIL: Could not identify a sequence-header column.")
  }

  primary <- primary |>
    mutate(
      variant = extract_variant(.data[[header_column]])
    )
}

if (any(is.na(primary$variant))) {
  stop("FAIL: Some primary sequences lack NDM variant labels.")
}

variant_counts_all <- primary |>
  count(variant, name = "sequence_count") |>
  arrange(
    as.numeric(str_extract(variant, "\\d+"))
  )

# ------------------------------------------------------------
# MUTATION BURDEN
# ------------------------------------------------------------

mutation_burden <- primary |>
  select(accession, variant) |>
  left_join(
    mutations |>
      count(accession, name = "mutation_count"),
    by = "accession"
  ) |>
  mutate(
    mutation_count = replace_na(mutation_count, 0L)
  ) |>
  arrange(
    as.numeric(str_extract(variant, "\\d+")),
    desc(mutation_count),
    accession
  )

# ------------------------------------------------------------
# OVERALL MUTATION FREQUENCY
# ------------------------------------------------------------

overall_frequency <- mutations |>
  count(mutation, name = "sequence_count") |>
  mutate(
    percentage = 100 * sequence_count / nrow(primary)
  ) |>
  arrange(desc(sequence_count), mutation)

# ------------------------------------------------------------
# MUTATION POSITION LANDSCAPE
# ------------------------------------------------------------

position_frequency <- mutations |>
  count(
    reference_position,
    name = "sequence_count"
  ) |>
  arrange(reference_position)

position_labels <- mutations |>
  count(
    reference_position,
    mutation,
    name = "sequence_count"
  ) |>
  group_by(reference_position) |>
  arrange(desc(sequence_count), mutation) |>
  slice_head(n = 1) |>
  ungroup()

# ------------------------------------------------------------
# VARIANT × MUTATION MATRIX
# ------------------------------------------------------------

variant_sizes <- primary |>
  count(variant, name = "sequence_count")

heatmap_data <- mutations |>
  count(variant, mutation, name = "mutation_count") |>
  complete(
    variant,
    mutation,
    fill = list(mutation_count = 0)
  ) |>
  left_join(
    variant_sizes,
    by = "variant"
  ) |>
  mutate(
    percentage = 100 * mutation_count / sequence_count
  )

# ------------------------------------------------------------
# CO-OCCURRENCE
# ------------------------------------------------------------

sequence_mutations <- mutations |>
  group_by(accession, variant) |>
  summarise(
    mutation_set = paste(
      sort(unique(mutation)),
      collapse = " + "
    ),
    mutation_count = n_distinct(mutation),
    .groups = "drop"
  )

cooccurrence <- sequence_mutations |>
  count(
    mutation_set,
    mutation_count,
    name = "sequence_count"
  ) |>
  arrange(
    desc(sequence_count),
    desc(mutation_count),
    mutation_set
  )

# ------------------------------------------------------------
# THEME
# ------------------------------------------------------------

theme_ndm <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    panel.grid.minor = element_blank(),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      10, 15, 10, 15
    )
  )

# ------------------------------------------------------------
# SAVE HELPER
# ------------------------------------------------------------

save_plot <- function(plot, filename, width = 8, height = 6) {

  png_file <- file.path(
    dirname(filename),
    paste0(basename(filename), ".png")
  )

  pdf_file <- file.path(
    dirname(filename),
    paste0(basename(filename), ".pdf")
  )

  ggsave(
    png_file,
    plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
}

# ------------------------------------------------------------
# FIGURE 1 — DATASET COMPOSITION
# ------------------------------------------------------------

p1_data <- variant_counts_all |>
  mutate(
    variant = factor(
      variant,
      levels = paste0(
        "NDM-",
        sort(unique(as.numeric(str_extract(variant, "\\d+"))))
      )
    ),
    label = paste0("n = ", sequence_count)
  )

p1 <- ggplot(
  p1_data,
  aes(
    x = variant,
    y = sequence_count
  )
) +
  geom_col(
    width = 0.72,
    fill = "#3B6EA8"
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.4,
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Composition of the 30-sequence NDM dataset",
    subtitle = "Curated primary sequences represented by each NDM variant; n denotes the number of sequences",
    x = "NDM variant",
    y = "Number of curated sequences"
  ) +
  theme_ndm +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

save_plot(
  p1,
  file.path(OUT_DIR, "figure1_variant_composition"),
  width = 10,
  height = 6
)

# ------------------------------------------------------------
# FIGURE 2 — MUTATION FREQUENCY

p2_data <- overall_frequency |>
  slice_head(n = 15) |>
  mutate(
    mutation = fct_reorder(mutation, sequence_count),
    label = paste0(
      "n = ",
      sequence_count,
      " (",
      number(
        100 * sequence_count / nrow(primary),
        accuracy = 0.1,
        trim = TRUE
      ),
      "%)"
    )
  )

p2 <- ggplot(
  p2_data,
  aes(
    x = mutation,
    y = sequence_count
  )
) +
  geom_col(
    width = 0.72,
    fill = "#2F6F9F"
  ) +
  geom_text(
    aes(
      label = label,
      y = sequence_count +
        case_when(
          mutation %in% c("D130G", "A233V") ~ 0.35,
          mutation %in% c("P28A", "D95N", "D130N") ~ 0.85,
          TRUE ~ 0.45
        )
    ),
    size = 3.0,
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Most frequently observed NDM mutations",
    subtitle = "Top 15 mutation types detected among the 30 curated sequences",
    x = "Amino-acid substitution",
    y = "Number of sequences carrying mutation"
  ) +
  theme_ndm +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

save_plot(
  p2,
  file.path(OUT_DIR, "figure2_top_mutations"),
  width = 10,
  height = 6
)


# FIGURE 3 — MUTATION POSITION LANDSCAPE

key_mutations <- c(
  "V88L",
  "D130G",
  "D130N",
  "M154L",
  "E170K",
  "A233V",
  "A263P",
  "R264H"
)

position_labels <- overall_frequency |>
  filter(mutation %in% key_mutations) |>
  left_join(
    mutations |>
      distinct(
        mutation,
        reference_position
      ),
    by = "mutation"
  ) |>
  mutate(
    label_x = case_when(
      mutation == "A263P" ~ reference_position - 3.5,
      mutation == "R264H" ~ reference_position + 3.5,
      TRUE ~ reference_position
    ),
    label_y = sequence_count + case_when(
      mutation == "A263P" ~ 0.9,
      mutation == "R264H" ~ 1.8,
      TRUE ~ 0.45
    ),
    label_hjust = case_when(
      mutation == "A263P" ~ 1,
      mutation == "R264H" ~ 0,
      TRUE ~ 0.5
    )
  )

p3 <- ggplot(
  position_frequency,
  aes(
    x = reference_position,
    y = sequence_count
  )
) +
  geom_segment(
    aes(
      xend = reference_position,
      y = 0,
      yend = sequence_count
    ),
    linewidth = 0.7,
    colour = "grey65"
  ) +
  geom_point(
    aes(
      fill = sequence_count
    ),
    shape = 21,
    size = 3.2,
    colour = "white",
    stroke = 0.5
  ) +
  geom_text(
    data = position_labels,
    aes(
      x = label_x,
      y = label_y,
      label = mutation,
      hjust = label_hjust
    ),
    size = 3.1,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis(
    option = "C",
    direction = -1,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 12),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Mutation landscape across the NDM protein",
    subtitle = "All observed substitutions are shown; positions are mapped to the CBZ39493.1 reference coordinate system",
    x = "Reference protein position",
    y = "Number of sequences with mutation"
  ) +
  theme_ndm

save_plot(
  p3,
  file.path(OUT_DIR, "figure3_mutation_landscape"),
  width = 11,
  height = 6
)


# FIGURE 4 — MUTATION BURDEN

mutation_burden_plot <- mutation_burden |>
  arrange(mutation_count, accession) |>
  mutate(
    sequence_label = paste0(
      accession,
      " (",
      variant,
      ")"
    ),
    sequence_label = factor(
      sequence_label,
      levels = sequence_label
    )
  )

median_burden <- median(
  mutation_burden_plot$mutation_count
)

p4 <- ggplot(
  mutation_burden_plot,
  aes(
    x = mutation_count,
    y = sequence_label
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = mutation_count,
      yend = sequence_label
    ),
    colour = "grey75",
    linewidth = 0.7
  ) +
  geom_point(
    size = 3.5,
    colour = "#2F6F9F"
  ) +
  geom_vline(
    xintercept = median_burden,
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#B44C3A"
  ) +
  annotate(
    "text",
    x = median_burden,
    y = Inf,
    label = paste0(
      "Median = ",
      median_burden
    ),
    vjust = 1.5,
    hjust = -0.05,
    fontface = "bold",
    colour = "#B44C3A"
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 10),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Mutation burden per NDM sequence",
    subtitle = "Each point represents one of the 30 curated primary sequences; zero-mutation sequences are retained",
    x = "Number of detected substitutions per sequence",
    y = "Sequence"
  ) +
  theme_ndm +
  theme(
    axis.text.y = element_text(
      size = 7.5
    )
  )

save_plot(
  p4,
  file.path(OUT_DIR, "figure4_mutation_burden"),
  width = 10,
  height = 10
)


# FIGURE 5 — VARIANT × MUTATION HEATMAP

top_heatmap_mutations <- overall_frequency |>
  slice_head(n = 20) |>
  pull(mutation)

p5_data <- heatmap_data |>
  filter(
    mutation %in% top_heatmap_mutations
  ) |>
  mutate(
    variant_number = as.numeric(
      str_extract(
        variant,
        "\\d+"
      )
    ),
    variant = factor(
      variant,
      levels = paste0(
        "NDM-",
        sort(
          unique(variant_number)
        )
      )
    )
  )

variant_labels <- primary |>
  count(
    variant,
    name = "sequence_count"
  ) |>
  mutate(
    label = paste0(
      variant,
      " (n = ",
      sequence_count,
      ")"
    )
  )

y_labels <- setNames(
  variant_labels$label,
  variant_labels$variant
)

p5 <- ggplot(
  p5_data,
  aes(
    x = mutation,
    y = variant,
    fill = percentage
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = ifelse(
        mutation_count > 0,
        paste0(
          mutation_count,
          "\n",
          number(
            percentage,
            accuracy = 0.1,
            trim = TRUE
          ),
          "%"
        ),
        ""
      ),
      colour = percentage >= 50
    ),
    size = 2.5
  ) +
  scale_colour_manual(
    values = c(
      "TRUE" = "white",
      "FALSE" = "black"
    ),
    guide = "none"
  ) +
  scale_fill_viridis(
    option = "D",
    direction = -1,
    limits = c(0, 100),
    oob = squish,
    name = "Within-variant\nfrequency (%)"
  ) +
  scale_y_discrete(
    labels = y_labels
  ) +
  labs(
    title = "Mutation profile across NDM variants",
    subtitle = "Counts and within-variant percentages for the 20 most frequent mutation types; zero indicates no mutation was observed in the selected dataset and does not imply biological absence",
    x = "Mutation",
    y = "NDM variant"
  ) +
  theme_ndm +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 8.5
    )
  )

save_plot(
  p5,
  file.path(OUT_DIR, "figure5_variant_mutation_heatmap"),
  width = 13,
  height = 8
)


# FIGURE 6 — PAIRWISE CO-OCCURRENCE

pairwise_cooccurrence <- mutations |>
  distinct(
    accession,
    mutation
  ) |>
  group_by(accession) |>
  summarise(
    mutation_list = list(sort(unique(mutation))),
    .groups = "drop"
  )

pairwise_rows <- lapply(
  seq_len(nrow(pairwise_cooccurrence)),
  function(i) {

    muts <- pairwise_cooccurrence$mutation_list[[i]]

    if (length(muts) < 2) {
      return(NULL)
    }

    pairs <- combn(
      muts,
      2,
      simplify = FALSE
    )

    data.frame(
      accession = pairwise_cooccurrence$accession[i],
      mutation_1 = vapply(
        pairs,
        function(x) x[1],
        character(1)
      ),
      mutation_2 = vapply(
        pairs,
        function(x) x[2],
        character(1)
      ),
      stringsAsFactors = FALSE
    )
  }
)

pairwise_rows <- Filter(
  Negate(is.null),
  pairwise_rows
)

if (length(pairwise_rows) == 0) {

  pairwise_cooccurrence <- tibble(
    accession = character(),
    mutation_1 = character(),
    mutation_2 = character(),
    sequence_count = integer(),
    mutation_pair = character()
  )

} else {

  pairwise_cooccurrence <- bind_rows(
    pairwise_rows
  ) |>
    count(
      mutation_1,
      mutation_2,
      name = "sequence_count"
    ) |>
    mutate(
      mutation_pair = paste(
        mutation_1,
        mutation_2,
        sep = " + "
      )
    ) |>
    arrange(
      desc(sequence_count),
      mutation_pair
    )
}

if (nrow(pairwise_cooccurrence) > 0) {

  top_pair_count <- max(
    pairwise_cooccurrence$sequence_count
  )

  if (top_pair_count >= 2) {

    p6_data <- pairwise_cooccurrence |>
      filter(
        sequence_count >= 2
      ) |>
      slice_head(n = 15)

  } else {

    p6_data <- pairwise_cooccurrence |>
      slice_head(n = 15)

  }

  p6_data <- p6_data |>
    mutate(
      mutation_pair = fct_reorder(
        mutation_pair,
        sequence_count
      )
    )

  p6 <- ggplot(
    p6_data,
    aes(
      x = sequence_count,
      y = mutation_pair
    )
  ) +
    geom_col(
      width = 0.7,
      fill = "#5B4B8A"
    ) +
    geom_text(
      aes(
        label = paste0(
          "n = ",
          sequence_count
        )
      ),
      hjust = -0.2,
      size = 3.2,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      title = "Recurring mutation pairs",
      subtitle = "Pairwise co-occurrence detected within the same NDM sequence; recurring pairs were observed in at least two sequences",
      x = "Number of sequences containing both mutations",
      y = "Mutation pair"
    ) +
    theme_ndm

} else {

  p6 <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No mutation pairs were detected",
      size = 5
    ) +
    labs(
      title = "Recurring mutation pairs",
      subtitle = "No sequences contained two or more detected mutations"
    ) +
    theme_void()

}

save_plot(
  p6,
  file.path(OUT_DIR, "figure6_mutation_cooccurrence"),
  width = 10,
  height = 7
)

# FIGURE 7 — INTEGRATED SUMMARY

metric_data <- tibble(
  metric = c(
    "Curated sequences",
    "NDM variants",
    "Mutation observations",
    "Mutation-bearing sequences",
    "Zero-mutation sequences",
    "Unique mutation types"
  ),
  value = c(
    nrow(primary),
    n_distinct(primary$variant),
    nrow(mutations),
    n_distinct(mutations$accession),
    sum(mutation_burden$mutation_count == 0),
    n_distinct(mutations$mutation)
  )
)

metric_data$metric <- factor(
  metric_data$metric,
  levels = rev(metric_data$metric)
)

p7_left <- ggplot(
  metric_data,
  aes(
    x = value,
    y = metric
  )
) +
  geom_col(
    fill = "#3B6EA8",
    width = 0.68
  ) +
  geom_text(
    aes(
      label = value
    ),
    hjust = -0.2,
    fontface = "bold",
    size = 4
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Dataset and mutation metrics",
    x = NULL,
    y = NULL
  ) +
  theme_ndm

p7_right_data <- overall_frequency |>
  slice_head(n = 8) |>
  mutate(
    mutation = fct_reorder(
      mutation,
      sequence_count
    )
  )

p7_right <- ggplot(
  p7_right_data,
  aes(
    x = sequence_count,
    y = mutation
  )
) +
  geom_col(
    fill = "#5B4B8A",
    width = 0.68
  ) +
  geom_text(
    aes(
      label = sequence_count
    ),
    hjust = -0.2,
    fontface = "bold",
    size = 3.5
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Most frequently observed substitutions",
    x = "Number of sequences",
    y = NULL
  ) +
  theme_ndm

p7 <- (
  p7_left +
  p7_right
) +
  plot_annotation(
    title = "Integrated overview of NDM sequence variation",
    subtitle = "Summary of the 30 curated sequences and their detected amino-acid substitutions"
  )

save_plot(
  p7,
  file.path(OUT_DIR, "figure7_ndm_mutation_summary"),
  width = 14,
  height = 8
)

report <- c(
  "NDM MUTATION VISUALISATION REPORT",
  paste(rep("=", 60), collapse = ""),
  "",
  paste("Mutation input:", MUTATION_FILE),
  paste("Primary metadata:", PRIMARY_FILE),
  paste("MSA:", MSA_FILE),
  "",
  paste("Primary sequences:", nrow(primary)),
  paste("Mutation observations:", nrow(mutations)),
  paste("Mutation-bearing accessions:", n_distinct(mutations$accession)),
  paste("Unique mutation types:", n_distinct(mutations$mutation)),
  paste("NDM variants:", n_distinct(primary$variant)),
  "",
  "NDM-16 check: PASS — A0A0K1Z5K4 = NDM-16 / R264H",
  "",
  "Figures generated:",
  "figure1_variant_composition",
  "figure2_top_mutations",
  "figure3_mutation_landscape",
  "figure4_mutation_burden",
  "figure5_variant_mutation_heatmap",
  "figure6_mutation_cooccurrence",
  "figure7_ndm_mutation_summary",
  "",
  "REPORT STATUS: PASS — Visualisation completed successfully."
)

REPORT_DIR <- file.path(BASE, "results", "reports")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

writeLines(
  report,
  file.path(REPORT_DIR, "visualisation_report.txt")
)

cat(paste(report, collapse = "\n"))
cat("\n")
