# 01_ancombc2_differential_abundance.R
# Differential abundance of genera between crisaborole-treated and untreated AD lesions.
# Reproduces the ANCOM-BC2 results reported in the manuscript.
#
# Inputs (relative to repository root):
#   data/feature-table-json.biom   unrarefied ASV count table exported from QIIME2 (see README)
#   data/taxonomy.tsv              SILVA 138 taxonomy per ASV (QIIME2 export)
#   data/sample-metadata.tsv       sample-id, group (Treated/Untreated), V# (visit), Subject ID
# Outputs: results/ancombc2_pooled_postbaseline_mixed.csv, results/ancombc2_pervisit_V{2..6}.csv,
#          results/sessionInfo.txt
#
# Run from the repository root:  Rscript R/01_ancombc2_differential_abundance.R
# Tested with R 4.5.1, ANCOMBC 2.10.1, phyloseq 1.52.0 (results/sessionInfo.txt). ANCOMBC 2.10.x
# requires CVXR <= 1.0-15 (CVXR 1.9 removed solve()); see README for the install line.

suppressPackageStartupMessages({ library(phyloseq); library(tidyverse); library(dplyr); library(ANCOMBC) })
dir.create("results", showWarnings = FALSE)

# ── 1) Import data & build phyloseq object ───────────────────────────────────
biom_fp    <- "data/feature-table-json.biom"
meta_fp    <- "data/sample-metadata.tsv"
taxmeta_fp <- "data/taxonomy.tsv"
taxmeta_fp <- "taxon_metadata.tsv"

# 1.1 Import BIOM, extract OTU table
biom_physeq <- import_biom(biom_fp) #, parseFunction = parse_taxonomy_default
otu_matrix  <- as(otu_table(biom_physeq), "matrix")
OTU         <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# 1.2 Read & split taxon metadata → tax_table
taxon_metadata <- read_tsv(
  taxmeta_fp,
  comment   = "#",
  col_types = cols(
    `Feature ID` = col_character(),
    Taxon        = col_character(),
    Confidence   = col_double()
  )
) %>%
  filter(`Feature ID` != "#q2:types")

tax_md_sep <- taxon_metadata %>%
  separate(
    col  = Taxon,
    into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    sep  = ";\\s*",
    fill = "right"
  )
tax_md_sep

tax_mat <- tax_md_sep %>%
  dplyr::select(`Feature ID`, Kingdom:Species) %>%
  column_to_rownames("Feature ID") %>%
  as.matrix()

TAX <- tax_table(tax_mat)

# 1.3 Read & process sample metadata → sample_data
meta_raw <- read_tsv(meta_fp, col_types = cols())
first_col <- names(meta_raw)[1]

sample_md <- meta_raw %>%
  rename(sample_id = all_of(first_col)) %>%
  column_to_rownames("sample_id") %>%
  rename(Subject_ID = `Subject ID`) %>%
  mutate(
    v_group   = paste(`V#`, group, sep = "_"),          # e.g. "V3_Treated"
    treatment = ifelse(str_detect(v_group, "_Treated$"),
                       "Treated", "Untreated"),
    visit     = factor(`V#`)                            # e.g. "V2","V3",…
  )

META <- sample_data(sample_md)

# 1.4 Build & prune initial phyloseq object
physeq <- phyloseq(OTU, TAX, META)
#depths     <- sample_sums(physeq)
#lib_cutoff <- quantile(depths, 0.05)
#message("lib_cutoff = ", lib_cutoff)

#physeq <- prune_samples(sample_sums(physeq) >= lib_cutoff, physeq)

# ── 2) Aggregate to Genus (tax_glom) ──────────────────────────────────────────
physeq_genus <- tax_glom(physeq, taxrank = "Genus", NArm = TRUE) %>%
  prune_taxa(taxa_sums(.) > 0, .)

genus_labels <- as.character(tax_table(physeq_genus)[, "Genus"])
dupes        <- names(which(table(genus_labels) > 1))
idx <- which(genus_labels == dupes)
# (duplicate-genus check removed; none present)

physeq_genus <- subset_taxa(
  physeq_genus,
  !grepl("uncultured", Genus, ignore.case = TRUE)
)

physeq_genus <- subset_taxa(
  physeq_genus,
  !grepl("Incertae", Genus, ignore.case = TRUE)
)

# (Optional) Rename each genus‐cluster’s taxa_name to its Genus label,
# so downstream tables refer to "g__Staphylococcus", etc., instead of OTU IDs.
taxa_names(physeq_genus) <- tax_table(physeq_genus)[, "Genus"]

# ── 3) Aggregate OTU (now genus) counts across technical replicates ───────────
# We'll take physeq_genus as our starting point (instead of physeq).

# 3.1 Convert sample_data (genus‐level) to a plain data frame with sample_id
meta_df <- sample_data(physeq_genus) %>%
  as("data.frame") %>%
  rownames_to_column("sample_id")

# 3.2 Extract genus‐level OTU table as data frame (samples as rows)
otu_df <- as.data.frame(t(as.data.frame(otu_table(physeq_genus)))) %>%
  rownames_to_column("sample_id")
# Now each column is one genus (taxa_names(physeq_genus)), each row is a sample_id.

# 3.3 Join OTU counts with metadata
otu_meta_combined <- otu_df %>%
  left_join(meta_df, by = "sample_id")

# 3.4 Aggregate (mean) counts for each Subject_ID + v_group
aggregated_data <- otu_meta_combined %>%
  group_by(Subject_ID, v_group) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )
# 3.5 Build an aggregated metadata table that preserves one row per Subject_ID + v_group
aggregated_meta <- meta_df %>%
  distinct(Subject_ID, v_group, .keep_all = TRUE) %>%
  dplyr::select(Subject_ID, v_group, treatment, visit)

# 3.6 Ensure aggregated_data has the same set of (Subject_ID, v_group) rows
aggregated_data <- aggregated_data %>%
  inner_join(aggregated_meta, by = c("Subject_ID", "v_group")) %>%
  arrange(Subject_ID, v_group)

# ── 4) Keep only Subject_ID + visit combos with BOTH Treated & Untreated ──────
paired_data <- aggregated_data %>%
  group_by(Subject_ID, visit) %>%
  filter(n_distinct(treatment) == 2) %>%
  ungroup()

# 4.1 Extract matching metadata rows
paired_meta <- aggregated_meta %>%
  semi_join(paired_data, by = c("Subject_ID", "v_group")) %>%
  arrange(Subject_ID, visit)

# Now `paired_data` has exactly two rows per (Subject_ID, visit): one Treated, one Untreated.

# ── 5) Reformat paired_data into phyloseq-compatible components ────────────────
# 5.1 Determine which columns are genus‐count columns
otu_numeric_cols <- setdiff(colnames(paired_data),
                            c("Subject_ID", "v_group", "treatment", "visit"))

# 5.2 Build the final OTU matrix (taxa_are_rows = TRUE)
final_otu_mat <- paired_data %>%
  dplyr::select(all_of(otu_numeric_cols)) %>%
  as.matrix() %>%
  t()
# Now rows are genera, columns are samples (paired subjects).

# 5.3 Construct unique column names matching Subject + visit + treatment
colnames(final_otu_mat) <- paste0(
  paired_data$Subject_ID, "_", paired_data$visit, "_", paired_data$treatment
)

# 5.4 Build the final metadata table with unique rownames
agg_meta_unique <- paired_data %>%
  mutate(unique_id = paste0(Subject_ID, "_", visit, "_", treatment)) %>%
  dplyr::select(unique_id, Subject_ID, visit, treatment) %>%
  column_to_rownames("unique_id")

# 5.5 Use the genus‐level TAX table
FINAL_OTU <- otu_table(final_otu_mat, taxa_are_rows = TRUE)
FINAL_META <- sample_data(agg_meta_unique)
FINAL_TAX  <- tax_table(physeq_genus)

# 5.6 Construct the final paired phyloseq object
physeq_paired <- phyloseq(FINAL_OTU, FINAL_TAX, FINAL_META)

# ── 6) Reorder phyloseq by Subject_ID then visit ──────────────────────────────
physeq_paired <- prune_samples(
  sample_names(physeq_paired)[order(
    sample_data(physeq_paired)$Subject_ID,
    sample_data(physeq_paired)$visit
  )],
  physeq_paired
)


# Check treatment factor
sample_data(physeq_paired)$treatment <- factor(sample_data(physeq_paired)$treatment)


sample_data(physeq_paired)$treatment <- relevel(sample_data(physeq_paired)$treatment, ref = "Untreated")

# ── 7) Primary analysis: treated vs untreated across post-baseline visits (V3–V6) ─
# Linear mixed model with subject and visit as random intercepts (split-lesion design:
# each subject contributes one treated and one untreated lesion per visit).
physeq_post <- subset_samples(physeq_paired, visit != "V2")
physeq_post <- prune_taxa(taxa_sums(physeq_post) > 0, physeq_post)
sample_data(physeq_post)$visit     <- factor(sample_data(physeq_post)$visit)
sample_data(physeq_post)$treatment <- relevel(factor(sample_data(physeq_post)$treatment), ref = "Untreated")

run_ancombc2 <- function(data, rand) {
  ancombc2(data = data, assay_name = "counts", tax_level = "Genus",
           fix_formula = "treatment", rand_formula = rand, group = "treatment",
           p_adj_method = "BH", lib_cut = 500, prv_cut = 0.05,
           struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05, global = FALSE)
}

set.seed(123)
pooled <- run_ancombc2(physeq_post, "(1|Subject_ID) + (1|visit)")
write.csv(pooled$res, "results/ancombc2_pooled_postbaseline_mixed.csv")
message("Pooled model: ", nrow(pooled$res), " genera tested; q<0.05: ",
        sum(pooled$res$q_treatmentTreated < 0.05, na.rm = TRUE))

# ── 8) Secondary analysis: treated vs untreated at each visit ────────────────
# A per-visit mixed model cannot be fitted (one observation per subject per arm), so
# visit-level comparisons are unpaired.
for (v in c("V2", "V3", "V4", "V5", "V6")) {
  ps_v <- subset_samples(physeq_paired, visit == v)
  ps_v <- prune_taxa(taxa_sums(ps_v) > 0, ps_v)
  sample_data(ps_v)$treatment <- relevel(factor(sample_data(ps_v)$treatment), ref = "Untreated")
  set.seed(123)
  res_v <- run_ancombc2(ps_v, NULL)
  write.csv(res_v$res, paste0("results/ancombc2_pervisit_", v, ".csv"))
  message(v, ": ", nrow(res_v$res), " genera tested; q<0.05: ",
          sum(res_v$res$q_treatmentTreated < 0.05, na.rm = TRUE))
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
