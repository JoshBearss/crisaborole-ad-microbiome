install.packages("qiime2R")

library(tidyverse)
library(phyloseq)
library(qiime2R)
library(vegan)
library(ggh4x)
library(ggbreak)
library(svglite)

# --- theme (same vibe as your script) ---
apatheme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(),
    text             = element_text(family = "Arial")
  )

taxa_colors <- c(
  "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770", "#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861", "red",
  "blue", "#483D8B", "gray", "darkcyan", "darkgoldenrod4", "green", "#8B0000",
  "darkgreen", "yellow"
)

results_dir <- "/Users/joshbearss/Documents/MesinkovskaLab/drive-download-20241215T185334Z-001/results"

# --- read the mixed table (counts + embedded metadata rows) ---
otu_raw <- readr::read_tsv(
  file.path(results_dir, "feature-table.tsv"),
  comment = "#",
  show_col_types = FALSE,
  progress = FALSE
)
names(otu_raw)[1] <- "Feature ID"

# --- extract embedded metadata (Treatment + subject) ---
treatment_row <- otu_raw %>%
  filter(`Feature ID` == "class") %>%
  select(-`Feature ID`) %>%
  pivot_longer(everything(), names_to = "SampleID", values_to = "Treatment") %>%
  mutate(Treatment = case_when(
    str_to_lower(Treatment) == "treated" ~ "T",
    str_to_lower(Treatment) == "untreated" ~ "U",
    TRUE ~ as.character(Treatment)
  ))

subject_row <- otu_raw %>%
  filter(`Feature ID` == "Subject ID") %>%
  select(-`Feature ID`) %>%
  pivot_longer(everything(), names_to = "SampleID", values_to = "subject")

# --- parse Timepoint from sample name (…-V2-…, …-V3-…, etc.) ---
sample_meta <- left_join(subject_row, treatment_row, by = "SampleID") %>%
  mutate(
    Timepoint = str_extract(SampleID, "V[0-9]+"),
    Date = str_extract(SampleID, "\\d{2}[A-Za-z]{3}\\d{4}")  # e.g., 06Apr2021
  )

# --- build counts matrix: drop the non-count rows and coerce to numeric ---
otu_counts <- otu_raw %>%
  filter(!`Feature ID` %in% c("class", "V", "Subject ID")) %>%
  mutate(across(-`Feature ID`, ~ readr::parse_number(.)))

otu_mat <- otu_counts %>%
  column_to_rownames("Feature ID") %>%
  as.matrix()

otu_mat[is.na(otu_mat)] <- 0

OTUtab <- otu_table(otu_mat, taxa_are_rows = TRUE)

# --- taxonomy from taxonomy.qza ---
tax_qza <- read_qza(file.path(results_dir, "taxonomy.qza"))
tax_mat <- parse_taxonomy(tax_qza$data) %>% as.matrix()
tax <- tax_table(tax_mat)

# --- sample_data: make rownames = sample IDs ---
md <- sample_meta %>%
  column_to_rownames("SampleID")

sam <- sample_data(md)

# --- build phyloseq object ---
dat <- phyloseq(OTUtab, tax, sam)
dat

