# =============================================================================
# Mesinkovska Skin Microbiome - Figure Recreation (NO qiime2R, NO EcolUtils)
# Uses:
#   - results/feature-table.tsv  (counts + embedded metadata rows)
#   - results/taxonomy.qza       (unzipped to taxonomy.tsv)
# Outputs recreated plots with *_rebuilt.svg names.
# =============================================================================

# ----------------------------
# 0) Packages
# ----------------------------
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggh4x)
library(ggbreak)
library(svglite)

# ----------------------------
# 1) Config
# ----------------------------
results_dir <- "/Users/joshbearss/Documents/MesinkovskaLab/drive-download-20241215T185334Z-001/results"
setwd(results_dir)

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

# Optional: a clean two-color palette for treatment plots
treat_colors <- c(T = "#DA5724", U = "#5F7FC7")

# ----------------------------
# 2) Load feature table TSV (contains embedded metadata rows)
# ----------------------------
otu_raw <- readr::read_tsv(
  file.path(results_dir, "feature-table.tsv"),
  comment = "#",
  show_col_types = FALSE,
  progress = FALSE
)
names(otu_raw)[1] <- "Feature ID"

# Extract embedded Treatment + subject rows
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

# Parse Timepoint and Date from sample IDs like: 001-06Apr2021-V3-T
sample_meta <- left_join(subject_row, treatment_row, by = "SampleID") %>%
  mutate(
    # Extract visit timepoint ONLY from the "-V#-" segment (avoids NOV2021 -> V2021)
    Timepoint = stringr::str_match(SampleID, "-(V\\d+)(?:-|$)")[, 2],
    
    # Extract date token like 01NOV2021 (optional; keeps your original intent)
    Date = stringr::str_extract(SampleID, "\\d{2}[A-Za-z]{3}\\d{4}"),
    
    # Optional: base ID without technical replicate suffix "-2"
    SampleID_base = stringr::str_replace(SampleID, "-2$", "")
  )


# Build counts matrix: drop the metadata rows and convert to numeric counts
otu_counts <- otu_raw %>%
  filter(!`Feature ID` %in% c("class", "V", "Subject ID")) %>%
  mutate(across(-`Feature ID`, ~ readr::parse_number(.)))

otu_mat <- otu_counts %>%
  column_to_rownames("Feature ID") %>%
  as.matrix()

otu_mat[is.na(otu_mat)] <- 0

OTUtab <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Sample data
md <- sample_meta %>% column_to_rownames("SampleID")
sam <- sample_data(md)

# ----------------------------
# 3) Read taxonomy.qza WITHOUT qiime2R
#    taxonomy.qza is a ZIP; we unzip and read the taxonomy.tsv inside.
# ----------------------------
if (!dir.exists("taxonomy_unzipped")) {
  unzip("taxonomy.qza", exdir = "taxonomy_unzipped")
}

tax_path <- file.path("taxonomy_unzipped", "data", "taxonomy.tsv")
tax_raw <- readr::read_tsv("taxonomy.tsv")

# Split QIIME2 Taxon string into ranks
# Typically "d__Bacteria; p__Firmicutes; c__...; ...; s__..."
tax_parsed <- tax_raw %>%
  select(`Feature ID`, Taxon) %>%
  separate(
    Taxon,
    into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    sep = ";\\s*",
    fill = "right",
    remove = TRUE
  ) %>%
  mutate(across(-`Feature ID`, ~ gsub("^.*__", "", .))) %>%  # remove rank prefixes
  column_to_rownames("Feature ID")

tax_mat <- as.matrix(tax_parsed)
tax <- tax_table(tax_mat)

# ----------------------------
# 4) Build phyloseq object
# ----------------------------
dat <- phyloseq(OTUtab, tax, sam)
dat

# =============================================================================
# FIGURES
# =============================================================================

# ----------------------------
# 5) Rarefaction curves
# ----------------------------
tab <- otu_table(dat)
class(tab) <- "matrix"
tab <- t(tab)

rar_curve_df <- vegan::rarecurve(
  tab,
  step  = 1000,
  label = FALSE,
  sample = 3514,
  xlab  = "Read Depth",
  col   = "orange",
  tidy  = TRUE
)

rar_curves <- rar_curve_df %>%
  filter(!Site %in% c("lysis-buffer", "lysis-buffer-2", "Mock", "Mock2", "negative", "negative2")) %>%
  ggplot(aes(x = Sample, y = Species, group = Site)) +
  geom_line() +
  geom_vline(xintercept = 3251, color = "red", linetype = "dashed") +
  labs(x = "Number of sequences", y = "Number of ASVs") +
  apatheme +
  theme(text = element_text(size = 12))

print(rar_curves)

ggsave("Rarefaction_curve_rebuilt.svg",
       rar_curves, width = 25, height = 20, units = "cm", dpi = 600, device = svglite)

rar_curves_zoom <- rar_curves +
  coord_cartesian(xlim = c(0, 10000)) #+
  #labs(subtitle = "Zoomed view (0–10,000 reads)")

print(rar_curves_zoom)

ggsave(
  "new_figs/Rarefaction_curve_zoomed.svg",
  rar_curves_zoom,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 600
)


# ----------------------------
# 6) Filter unwanted taxa (Chloroplast / Mito / Euk / Unassigned)
# ----------------------------
dat2 <- dat %>%
  subset_taxa(Genus != "Chloroplast" | is.na(Genus)) %>%
  subset_taxa(Genus != "Mitochondria" | is.na(Genus)) %>%
  subset_taxa(Kingdom != "Eukaryota" | is.na(Kingdom)) %>%       # after prefix strip above
  subset_taxa(Kingdom != "Unassigned" | is.na(Kingdom))

# ----------------------------
# 7) Sequence count by sample (pre-rarefaction, like your script)
# ----------------------------
sampsums <- tibble(
  sample_id = sample_names(dat2),
  count = as.numeric(sample_sums(dat2))
) %>%
  arrange(count) %>%
  mutate(sample_id = fct_reorder(sample_id, count))

seq_count_plot <- ggplot(sampsums, aes(x = sample_id, y = count)) +
  geom_col() +
  geom_hline(yintercept = 3251, color = "red", linetype = "dashed") +
  coord_cartesian(ylim = c(0, 20000)) +   # zoom instead of break
  labs(x = "Sample ID", y = "Sequence count") +
  apatheme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

print(seq_count_plot)


ggsave("seq_count_by_sample_rebuilt.svg",
       seq_count_plot, width = 50, height = 25, units = "cm", dpi = 600, device = svglite)

# ----------------------------
# 8) Rarefy (replacement for EcolUtils::rrarefy.perm)
# ----------------------------
set.seed(1)
rar_data <- rarefy_even_depth(
  dat2,
  sample.size = 3251,
  rngseed = 1,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = TRUE
)

rar_data

# Remove any samples without Timepoint (if present)
rar_data_nona <- subset_samples(rar_data, !is.na(Timepoint))

# ----------------------------
# 9) Alpha diversity (Shannon + Observed)
# ----------------------------
divdf <- estimate_richness(rar_data, measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson"))
divdf <- cbind(divdf, as.data.frame(sample_data(rar_data)))

# Shannon by Treatment
p_shannon_treat <- ggplot(divdf, aes(Treatment, Shannon, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.8) +
  scale_fill_manual(values = treat_colors) +
  apatheme

print(p_shannon_treat)

ggsave("shannon_by_treatment_rebuilt.svg",
       p_shannon_treat, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Shannon by Treatment + subject
p_shannon_sub <- ggplot(divdf, aes(subject, Shannon, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), alpha = 0.8) +
  scale_fill_manual(values = treat_colors) +
  apatheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_shannon_sub)

ggsave("shannon_by_treatmentANDsubject_rebuilt.svg",
       p_shannon_sub, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Shannon by Treatment + Timepoint faceted by subject
p_shannon_time_sub <- divdf %>%
  filter(!is.na(Timepoint)) %>%
  ggplot(aes(Timepoint, Shannon, color = Treatment, group = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~subject) +
  scale_color_manual(values = treat_colors) +
  apatheme

print(p_shannon_time_sub)

ggsave("shannon_by_treatmentANDsubjectANDtimepoint_rebuilt.svg",
       p_shannon_time_sub, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Shannon means +/- SD by Treatment + Timepoint
divsum <- divdf %>%
  filter(!is.na(Timepoint)) %>%
  group_by(Timepoint, Treatment) %>%
  summarize(avg = mean(Shannon), sd = sd(Shannon), n = n(), .groups = "drop")

p_shannon_means <- ggplot(divsum, aes(Timepoint, avg, color = Treatment, group = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 0.2) +
  ylab("Shannon") +
  scale_color_manual(values = treat_colors) +
  apatheme

print(p_shannon_means)

ggsave("shannon_by_treatmentANDtimepoint_means_rebuilt.svg",
       p_shannon_means, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Observed richness by Treatment
p_obs_treat <- ggplot(divdf, aes(Treatment, Observed, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.8) +
  scale_fill_manual(values = treat_colors) +
  apatheme

print(p_obs_treat)

ggsave("richness_by_treatment_rebuilt.svg",
       p_obs_treat, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Observed by Treatment + subject
p_obs_sub <- ggplot(divdf, aes(subject, Observed, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), alpha = 0.8) +
  scale_fill_manual(values = treat_colors) +
  apatheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_obs_sub)

ggsave("richness_by_treatmentANDsubject_rebuilt.svg",
       p_obs_sub, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Observed by Treatment + Timepoint faceted by subject
p_obs_time_sub <- divdf %>%
  filter(!is.na(Timepoint)) %>%
  ggplot(aes(Timepoint, Observed, color = Treatment, group = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~subject) +
  scale_color_manual(values = treat_colors) +
  apatheme

print(p_obs_time_sub)

ggsave("richness_by_treatmentANDsubjectANDtimepoint_rebuilt.svg",
       p_obs_time_sub, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# Observed means +/- SD by Treatment + Timepoint
divsumO <- divdf %>%
  filter(!is.na(Timepoint)) %>%
  group_by(Timepoint, Treatment) %>%
  summarize(avg = mean(Observed), sd = sd(Observed), n = n(), .groups = "drop")

p_obs_means <- ggplot(divsumO, aes(Timepoint, avg, color = Treatment, group = Treatment, shape = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 0.2) +
  ylab("Observed") +
  scale_color_manual(values = treat_colors) +
  apatheme

print(p_obs_means)

ggsave("richness_by_treatmentANDtimepoint_means_rebuilt.svg",
       p_obs_means, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# ----------------------------
# 10) Taxa barplot (Genus)
#     (same filtering concept: keep genera with mean rel abundance >= 0.0075)
# ----------------------------
genus_data <- tax_glom(rar_data_nona, taxrank = "Genus")
phy_g <- transform_sample_counts(genus_data, \(x) x / sum(x))

# keep genera with mean rel abundance >= 0.0075
keep_taxa <- taxa_names(phy_g)[
  taxa_names(phy_g) %in% taxa_names(filter_taxa(phy_g, \(x) mean(x) >= 0.0075, TRUE))
]

phy_keep <- prune_taxa(keep_taxa, phy_g)
df_genus <- psmelt(phy_keep) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("V2","V3","V4","V5","V6"))) %>%
  arrange(Genus)

# order samples by subject, Treatment, Timepoint (like your original facet ordering)
sample_order <- df_genus %>%
  distinct(Sample, subject, Treatment, Timepoint) %>%
  arrange(subject, Treatment, Timepoint) %>%
  pull(Sample)

sample_timepoint_labels <- df_genus %>%
  distinct(Sample, Timepoint) %>%
  deframe()

genus_names <- sort(unique(df_genus$Genus))
if (length(genus_names) >= 2) genus_names[2] <- "Rhizobium"

p_genus <- ggplot(df_genus, aes(x = factor(Sample, levels = sample_order), y = Abundance, fill = Genus)) +
  geom_col() +
  scale_fill_manual(values = taxa_colors, labels = genus_names) +
  facet_nested(~ subject + Treatment, scale = "free", space = "free_x") +
  scale_x_discrete(labels = sample_timepoint_labels) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Timepoints") +
  ylab("Abundance (%)") +
  theme_bw() +
  theme(
    panel.spacing = unit(0.125, "line"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

print(p_genus)

ggsave("new_figs/taxa_barplot_genus_rebuilt.svg",
       p_genus, width = 100, height = 25, units = "cm", dpi = 600, device = svglite)

sample_meta %>%
  filter(is.na(Timepoint)) %>%
  select(SampleID, subject, Treatment) %>%
  distinct()
sample_meta %>%
  filter(subject == "016") %>%
  select(SampleID, Timepoint)

# =========================
# Readable Genus Barplot
# - Side-by-side alignment restored
# - Nested facets: subject + Treatment
# - Top 10 genera + "Other"
# =========================

library(tidyverse)
library(ggplot2)
library(ggh4x)

# --- 1) Keep top N genera and lump the rest as "Other"
top_n <- 8

top_genera <- df_genus %>%
  group_by(Genus) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = top_n) %>%
  pull(Genus)

df_genus_readable <- df_genus %>%
  mutate(
    Genus = if_else(Genus %in% top_genera, as.character(Genus), "Other"),
    Genus = factor(Genus, levels = c(sort(top_genera), "Other"))
  )

# --- 2) Preserve sample ordering (subject → Treatment → Timepoint)
df_genus_readable <- df_genus_readable %>%
  mutate(Sample = factor(Sample, levels = sample_order))

# --- 3) Build a safe fill palette (includes "Other")
genus_levels <- levels(df_genus_readable$Genus)
fill_pal <- setNames(rep(taxa_colors, length.out = length(genus_levels)), genus_levels)
fill_pal["Other"] <- "grey80"

# --- 4) Plot (SIDE-BY-SIDE, nested facets)
p_genus_readable <- ggplot(
  df_genus_readable,
  aes(x = Sample, y = Abundance, fill = Genus)
) +
  geom_col(width = 0.95) +
  facet_nested(
    ~ subject + Treatment,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values = fill_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Sequential visits per subject",
    y = "Relative abundance (%)",
    fill = "Genus"
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0.25, "lines"),
    panel.border  = element_rect(color = "gray70", fill = NA, linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(size = 9),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm")
  )

print(p_genus_readable)

ggsave(
  "taxa_barplot_genus_side_by_side_readable.svg",
  plot  = p_genus_readable,
  width = 50,
  height = 25,
  units  = "cm",
  dpi    = 600
)



# ----------------------------
# 11) NMDS plots (Bray-Curtis)
# ----------------------------
nmds <- ordinate(rar_data_nona, method = "NMDS", distance = "bray")
nmds_df <- plot_ordination(rar_data_nona, nmds, justDF = TRUE) %>%
  bind_cols(as.data.frame(sample_data(rar_data_nona)))

p_nmds_treat <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = Treatment)) +
  geom_point(size = 4, alpha = 0.75) +
  scale_color_manual(values = treat_colors) +
  annotate("text", x = -0.4, y = -0.2,
           label = paste("Stress =", round(nmds$stress, 3)), size = 3.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(p_nmds_treat)

ggsave("nmds_by_treatment_rebuilt.svg",
       p_nmds_treat, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# NMDS by subject + Treatment as shape
subjects <- sort(unique(nmds_df$subject))
subj_cols <- setNames(colorRampPalette(taxa_colors)(length(subjects)), subjects)

p_nmds_sub <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = subject, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.75) +
  scale_color_manual(values = subj_cols) +
  annotate("text", x = -0.4, y = -0.2,
           label = paste("Stress =", round(nmds$stress, 3)), size = 3.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(p_nmds_sub)

ggsave("nmds_by_subject_rebuilt.svg",
       p_nmds_sub, width = 20, height = 20, units = "cm", dpi = 600, device = svglite)

# NMDS faceted by subject with Timepoint labels
p_nmds_facet <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_text(aes(label = Timepoint), size = 2.5, vjust = 1.2, show.legend = FALSE) +
  scale_color_manual(values = treat_colors) +
  facet_wrap(~subject) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(p_nmds_facet)

ggsave("nmds_by_treatment_facetby_subject_rebuilt.svg",
       p_nmds_facet, width = 35, height = 30, units = "cm", dpi = 600, device = svglite)

# =============================================================================
# End
# =============================================================================





library(tidyverse)
library(ggplot2)

# df_genus is your melted phyloseq table (psmelt output) with columns:
# Sample, Abundance, Genus, subject, Treatment, Timepoint

top_n <- 8

# Ensure Timepoint is ordered
df_genus2 <- df_genus %>%
  filter(!is.na(Timepoint)) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("V2","V3","V4","V5","V6")))

# Top genera (by mean abundance across all samples)
top_genera <- df_genus2 %>%
  group_by(Genus) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = top_n) %>%
  pull(Genus)

# Lump the rest into "Other"
df_genus2 <- df_genus2 %>%
  mutate(Genus = if_else(Genus %in% top_genera, as.character(Genus), "Other"))

# Compute mean composition per Treatment x Timepoint x Genus
df_mean <- df_genus2 %>%
  group_by(Treatment, Timepoint, Genus) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Treatment, Timepoint) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%   # renormalize to 100%
  ungroup()

# Build a safe palette for the genera shown
genus_levels <- df_mean %>%
  distinct(Genus) %>%
  arrange(Genus) %>%
  pull(Genus)

fill_pal <- setNames(rep(taxa_colors, length.out = length(genus_levels)), genus_levels)
if ("Other" %in% names(fill_pal)) fill_pal["Other"] <- "grey80"

# Plot: grouped by Treatment, x = Timepoint (labels visible!)
p_genus_grouped <- ggplot(df_mean, aes(x = Timepoint, y = Abundance, fill = Genus)) +
  geom_col(width = 0.85) +
  facet_wrap(~ Treatment, nrow = 1) +
  scale_fill_manual(values = fill_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Timepoint",
    y = "Mean relative abundance (%)",
    fill = "Genus",
    title = "Genus-level composition over time",
    subtitle = "Stacked bars show mean composition by Treatment at each visit (Top genera + Other)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 11)
  )

print(p_genus_grouped)

ggsave(
  "new_figs/taxa_barplot_genus_grouped_Treatment_Timepoint.svg",
  plot = p_genus_grouped,
  width = 26,
  height = 14,
  units = "cm",
  dpi = 600
)


# =============================================================================
# OPTION 2 (FIXED): Per-subject detailed stacked bars
# - x = Timepoint (V2–V6) so Treated + Untreated align vertically
# - facet_grid(Treatment ~ subject)
# - Top N genera + Other
# - Aggregates technical replicates (-2) by averaging within each Timepoint
# =============================================================================

library(tidyverse)
library(ggplot2)

top_n <- 5

# 1) Clean/ensure Timepoint is correct and ordered
df2 <- df_genus %>%
  filter(!is.na(Timepoint)) %>%
  mutate(
    Timepoint = factor(Timepoint),
    Timepoint = factor(
      Timepoint,
      levels = sort(unique(Timepoint)),
      labels = seq_along(sort(unique(Timepoint)))
    )
  )

# 2) Keep top genera + lump rest as Other (improves readability)
top_genera <- df2 %>%
  group_by(Genus) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = top_n) %>%
  pull(Genus)

df2 <- df2 %>%
  mutate(
    Genus = if_else(Genus %in% top_genera, as.character(Genus), "Other"),
    Genus = factor(Genus, levels = c(sort(top_genera), "Other"))
  )

# 3) Aggregate to ONE bar per subject x Treatment x Timepoint (averages replicates like "-2")
df_agg <- df2 %>%
  group_by(subject, Treatment, Timepoint, Genus) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(subject, Treatment, Timepoint) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%  # renormalize to 100%
  ungroup()

# 4) Safe fill palette (includes Other)
genus_levels <- levels(df_agg$Genus)
fill_pal <- setNames(rep(taxa_colors, length.out = length(genus_levels)), genus_levels)
fill_pal["Other"] <- "grey80"

# 5) Plot: Treatment rows, subject columns, Timepoint x-axis (aligned!)
p_option2_aligned <- ggplot(df_agg, aes(x = Timepoint, y = Abundance, fill = Genus)) +
  geom_col(width = 0.85) +
  facet_grid(Treatment ~ subject) +
  scale_fill_manual(values = fill_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Timepoint",
    y = "Relative abundance (%)",
    fill = "Genus",
    title = "Genus-level composition by subject over time",
    subtitle = "Treated and Untreated aligned by timepoint within each subject (Top genera + Other)"
  ) +
  theme_bw(base_size = 20) +   # <-- BIG base size
  theme(
    # Panels
    panel.spacing = unit(0.5, "lines"),
    panel.border  = element_rect(color = "gray70", fill = NA, linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Axis text
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 12)),
    
    # Facet strips
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    strip.text = element_text(size = 18, face = "bold"),
    
    # Titles
    plot.title = element_text(size = 26, face = "bold"),
    plot.subtitle = element_text(size = 20),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.key.height = unit(0.7, "cm"),
    legend.key.width = unit(0.6, "cm")
  )

print(p_option2_aligned)

ggsave(
  "new_figs/taxa_barplot_genus_option2_aligned_AAD.svg",
  plot  = p_option2_aligned,
  width = 55,
  height = 22,
  units  = "cm",
  dpi    = 600
)

vivid_palette <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#984EA3",  # purple
  "#FF7F00",  # orange
  "#FFFF33",  # yellow
  "#A65628",  # brown
  "#F781BF"   # pink
)

genus_levels <- levels(df_agg$Genus)

fill_pal <- setNames(
  vivid_palette[seq_along(genus_levels)],
  genus_levels
)

fill_pal["Other"] <- "grey85"
