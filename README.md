# Topical crisaborole and the atopic dermatitis skin microbiome — analysis code and data

Code and processed data for: *[Title]*, Bearss JD, Greene RK, Phong C, Nguyen C, Babadjouni A,
Mesinkovska NA, Juhasz M. Submitted to the *Journal of Investigative Dermatology* (Letter to the Editor).

Twenty patients with atopic dermatitis applied 2% crisaborole ointment to one lesion and left a
second lesion untreated. Both lesions were swabbed at baseline (V2), week 2 (V3), week 4 (V4),
week 12 (V5, end of treatment) and week 16 (V6, 4 weeks after cessation). 16S rRNA amplicon
sequencing was performed by the UCI Microbiome Center.

## Data availability
Raw sequence reads: NCBI SRA BioProject **PRJNA______** (accession to be inserted before publication).
Processed data in this repository:

| File | Description |
|---|---|
| `data/feature-table-json.biom` | Unrarefied ASV × sample count table exported from QIIME2 (171 samples incl. 6 controls, 5,545 ASVs). Input to ANCOM-BC2. |
| `data/taxonomy.tsv` | SILVA 138 taxonomy and confidence per ASV. |
| `data/sample-metadata.tsv` | De-identified sample metadata: sample ID, arm (Treated/Untreated), visit, subject. |
| `qiime2/table1.qza`, `filtered-table.qza`, `taxonomy.qza` | QIIME2 artifacts from the core with full provenance (DADA2 parameters, classifier). |
| `core_report/core_feature_table.tsv` | Core's exported feature table used for the rarefied diversity analyses. |

Sample IDs follow `subject-date-visit-arm` (e.g. `015-26OCT2021-V2-T`); a `-2` suffix marks a
re-amplified technical replicate. Controls: `Mock`, `Mock2` (ZymoBIOMICS D6305), `negative`,
`negative2` (PCR no-template), `lysis-buffer`, `lysis-buffer-2` (extraction blanks).

## Sequencing and upstream processing (UCI Microbiome Center)
DNA: ZymoBIOMICS 96 MagBead kit. PCR: Earth Microbiome Project 515F-Y/926R (V4–V5, ~411 bp),
AccuStart II ToughMix, 30 cycles. Sequencing: Illumina MiSeq 2×300 v3. Only forward reads were
usable. QIIME2 2022.2: `demux emp-single` → `dada2 denoise-single` (trim-left 5, trunc-len 211,
max-ee 2, chimera consensus) → `feature-classifier classify-sklearn` with a naive Bayes classifier
trained on SILVA 138 SSURef NR99 extracted with 515F/926R (RESCRIPt), confidence 0.7.
All parameters are recoverable from the provenance inside the `.qza` files (`qiime tools peek`
or view.qiime2.org).

## Analyses
### Diversity, ordination and PERMANOVA (`core_report/`)
Performed by the UCI Microbiome Center in R 4.1: rarefaction to 3,251 sequences (EcolUtils,
100 iterations), Shannon index and richness, Bray–Curtis NMDS, PERMANOVA in PRIMER-e with
treatment nested in subject. Scripts: `skin_alpha_beta_permanova.Rmd`, `MicrobiomeFigures.R`,
`Figures2.R` (the Rmd also merges an extended clinical metadata sheet that is not distributed, so the
confounder columns in the PERMANOVA are not re-runnable from this repository); PERMANOVA output in `PERMANOVA_subject_treatment_time.xls`. Final table after QC:
4,426 ASVs, 174 experimental samples + 6 controls.

### Differential abundance (`R/01_ancombc2_differential_abundance.R`)
ANCOM-BC2 on unrarefied counts aggregated to genus (uncultured / *Incertae sedis* removed),
technical replicates averaged, restricted to subject–visit pairs with both arms (142 samples).
`lib_cut = 500`, `prv_cut = 0.05`, BH-adjusted q-values.

* **Primary:** treated vs untreated across post-baseline visits V3–V6 (118 samples), linear mixed
  model with subject and visit random intercepts → `results/ancombc2_pooled_postbaseline_mixed.csv`
* **Secondary:** treated vs untreated at each visit, unpaired (a per-visit mixed model has one
  observation per subject per arm and cannot be fitted) → `results/ancombc2_pervisit_V2..V6.csv`
* **Sensitivity:** subject-only random intercept and no random effect for the pooled comparison
  → `results/sensitivity/`. `R/99_sensitivity_all_variants.R` reruns every variant and writes a log.

`results/run_log_2026-09-17.txt` and `results/sessionInfo.txt` record the run that produced the
reported numbers.

## Reproducing
```r
# R >= 4.5. ANCOMBC 2.10.x needs CVXR <= 1.0-15 (CVXR 1.9 removed solve()):
install.packages("CVXR", repos = "https://packagemanager.posit.co/cran/2025-06-02")
BiocManager::install(c("ANCOMBC", "phyloseq")); install.packages("tidyverse")
```
```sh
Rscript R/01_ancombc2_differential_abundance.R     # ~15 min; mixed models dominate
```

## Citation
Please cite the manuscript and this repository (Zenodo DOI: 10.5281/zenodo.______).
16S sequencing and initial bioinformatic analysis were performed by the UCI Microbiome Center.
Funded by Pfizer Inc.; the sponsor had no role in study design, analysis, or publication.

## License
Code: MIT (see `LICENSE`). Data and results: CC BY 4.0.
