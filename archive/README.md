# Archive

Earlier code kept for provenance. Nothing here is run by `make`. The paths inside these files
point to their original locations and will not resolve in this repository.

## `pilot_report_2026/`

The first version of the pilot analysis (April-October 2026), an R Markdown report.

| File | Replaced by |
|---|---|
| `methylation_report.Rmd`, `render.R` | `scripts/02_pilot_methylation/01_preprocess_sesame.R` (preprocessing, QC, PCA) and `03_exploratory_dmp_dmr.R` (paired DMP/DMR candidates) |
| `install_deps.R` | `env/renv.lock` |
| `17_story_figures.R` | first-draft figures, superseded by `scripts/07_figures/` |

The report applied sesame `prepSesame("QCDPB")` twice (see `docs/CHANGELOG.md`).

## `gse204943_original_scripts/`

The scripts from the earlier GSE204943 project that made the frozen objects
(`crucial_objects/beta.rds`, `pd.rds`, `ann_auto.rds` and the `results/` tables) used before
stage 01 existed. `scripts/01_gse204943/` reimplements the steps the manuscript depends on.

| Original script | Output used here | Reimplemented in |
|---|---|---|
| `01_data_processing_QC.R` | `beta.rds`, `pd.rds`, annotation | `01_minfi_preprocess.R` |
| `25_cell_type_deconvolution.R` | `cell_type_proportions.csv` | `02_cell_composition.R` |
| `10_PMD_erosion_metrics.R`, `16_survival_analysis.R` | `hypo_PMD_blocks.csv`, `PMD_erosion_score_per_sample.csv` | `03_pmd_blocks.R` |
| `27_CNV_inference.R` | `CNV_segments_all.csv` | `04_copy_number_conumee2.R` |
| `07_locus_specific_TFAP2A_ZNF471.R` | `ZNF471_CpG_details.csv` | `05_locus_cpg_sets.R` |

Known problems in these scripts, fixed in the reimplementation:

- `01_data_processing_QC.R` applies the probe filter to the RGChannelSet (see `docs/CHANGELOG.md`).
- `27_CNV_inference.R` uses `tot_int_f` without defining it. It was presumably the total
  (methylated + unmethylated) intensity matrix from an interactive session.
- `10_PMD_erosion_metrics.R` was run interactively in sections, and `16_survival_analysis.R`
  reads objects from that session. The per-sample erosion score is defined in `16`.
