# Pipeline: scripts, inputs and outputs

Paths are relative to the configured locations: `<out>` = `outputs` (config entry `outputs`),
`<ext>` = `<out>/extended`, `<gse>` = GSE204943 objects (`gse204943_objects`), `<rna>` = pilot
RNA-seq folder (`rnaseq`), `<geo>` = public downloads (`geo`). Tables go to `<ext>/tables/`.
The sentinel column gives the file that Make uses to decide whether a script must rerun.

## 00 Download

| Script | Inputs | Outputs |
|---|---|---|
| `00_download/download_geo.sh` | internet (GEO FTP) | `<geo>/`: GSE85514, GSE85195, GSE26549, GSE234379 series matrices; GSE234379 normalised betas; GPL22311 and GPL6480 tables; GSE227919 files; GSE204943 sample sheet and IDATs |

## 01 GSE204943 (EPIC v1, Mumbai)

| Script | Inputs | Main outputs (sentinel in bold) |
|---|---|---|
| `01_minfi_preprocess.R` | IDATs, sample sheet | **`<out>/gse204943/crucial_objects/beta.rds`**, `pd.rds`, `ann_auto.rds`; `intermediate/RGset_qc.rds`, `Mset_noob.rds` |
| `02_cell_composition.R` | beta, pd | **`results/cell_type_proportions.csv`** |
| `03_pmd_blocks.R` | `Mset_noob.rds`, pd | **`results/hypo_PMD_blocks.csv`**, `PMD_erosion_score_per_sample.csv` |
| `04_copy_number_conumee2.R` | `RGset_qc.rds`, pd | **`results/CNV_segments_all.csv`** |
| `05_locus_cpg_sets.R` | beta, pd, annotation | **`results/ZNF471_CpG_details.csv`** |

When `gse204943_objects` points to a frozen copy outside `<out>`, stage 01 is not part of
`make all`, and downstream scripts read the frozen files.

## 02 Pilot methylation (EPIC v2)

| Script | Inputs | Main outputs |
|---|---|---|
| `01_preprocess_sesame.R` | pilot IDATs, sample sheet | **`<out>/processed.rds`** (betas, M values, uncollapsed betas, sample sheet, PCA, QC); `qc_summary.csv` |
| `02_identity_sex_copy_number.R` | pilot IDATs | **`<ext>/cnv_objects.rds`**; genotype concordance, sex, arm-level and 9p21 copy-number tables |
| `03_exploratory_dmp_dmr.R` (`make exploratory`) | `processed.rds` | **`dmp_lesion_vs_normal.csv`**, `dmrs_lesion_vs_normal.csv` (not used in the manuscript) |

## 03 Methylation cohorts

| Script | Inputs | Main outputs |
|---|---|---|
| `01_composition_signature_scores.R` | `processed.rds`; `<gse>` objects and results | **`<ext>/external_methylation.rds`**, `gse_signature.rds`; EpiDISH fractions, signature CpGs, tumour-likeness (GSE204943 and pilot), PCA projection, replication, PMD/ZNF471 tables |
| `02_opl_groups.R` | `external_methylation.rds`, `<gse>/results/CNV_segments_all.csv`, GSE227919 | **`<ext>/story_analyses.rds`**; OPL methylome groups vs copy number; GSE227919 arm dosage |
| `03_mumbai_crosscheck.R` | above, GSE85514, GSE85195, GPL tables | **`<ext>/mumbai_crosscheck.rds`**; conumee vs array CGH, aCGH arms by OPL group, methylation-expression correlations |
| `04_gse234379_kzfp.R` | above, GSE234379 betas, `processed.rds` | **`<ext>/gse234379_kzfp.rds`**; signature transfer, paired tests, KZFP by OPL group, pilot KZFP deltas |

## 04 Expression cohorts and outcome

| Script | Inputs | Main outputs |
|---|---|---|
| `01_gse227919_programmes.R` | GSE227919, `external_methylation.rds`, `processed.rds` | **`<ext>/external_rnaseq.rds`**; programme and gene tests, methylation-linked gene sets |
| `02_gse26549_transformation.R` | GSE26549, `external_rnaseq.rds`, `mumbai_crosscheck.rds` | **`<ext>/gse26549_progression.rds`**; methylome-state proxy genes, Cox models, random-set and label-permutation nulls (~10 min) |

## 05 Pilot RNA-seq

| Script | Inputs | Main outputs |
|---|---|---|
| `00_strandedness_check.sh` (`make kallisto`) | first 300,000 read pairs of S1 | `<rna>/strandedness_check.tsv` |
| `01_kallisto_sequential.sh` (`make kallisto`) | FASTQs, kallisto index | `<rna>/kallisto/<library>/abundance.tsv.gz`, `run_info.json` |
| `02_gene_matrices.R` | kallisto outputs | **`<rna>/matrices/gene_tpm.tsv.gz`**, gene counts, `samples.tsv`, `gene_annotation.tsv` |
| `03_rna_methylation_validation.R` | matrices, `cnv_objects.rds`, `external_methylation.rds`, GSE26549 and KZFP outputs | **`rnaseq_proxy_programmes.csv`**; QC/sex, composition vs EpiDISH, arm dosage vs copy number, 9p21, KZFP |
| `04_rna_analyses.R` | matrices, `processed.rds`, stage 03-04 outputs | **`<ext>/rnaseq_story.rds`**; programme scores, paired DE and hallmarks, promoter methylation vs expression, chr19 KZFP enrichment, external test of the lesion programme |

## 06 Audit

| Script | Inputs | Main outputs |
|---|---|---|
| `01_results_audit.R` | stage 02-04 outputs, `<gse>` | **`<ext>/audit.rds`**; `audit_checks.csv` (key claims: reported vs recomputed, status `confirmed` or `CHECK`), `audit_gse26549_forest_exact_ci.csv`; OPL score gap, mixture BIC and dip test |
| `02_cohort_verification.R` | stage 02-04 outputs, `<gse>` | **`verification_kzfp_cpg_accounting.csv`**; Mumbai OPL outcome by group (Supplementary Table 2), sex-adjusted signature, HPV-negative subset, tumour-only PMD blocks, GSE227919 subsets |
| `03_manuscript_numbers.R` | audit, verification and stage 02-05 outputs | **`manuscript_numbers.csv`**: every number quoted in the text (section, statement, value), formatted as in the text |

## 07 Figures

| Script | Inputs | Outputs (`<out>/figures_nature/`, PDF + PNG) |
|---|---|---|
| `01_figures_1_to_4.R` | stage 02-04 and audit outputs | **`Fig1`**, `Fig2`, `Fig3`, `Fig4`, `Supplementary_Fig1`, `Supplementary_Fig2` |
| `02_figures_5_6.R` | stage 05 outputs, GSE26549, audit, `resources/cnv_events_curated.csv` | **`Fig5`**, `Fig6`, `Supplementary_Fig3` |

## Shared code

| File | Content |
|---|---|
| `R/lib/config.R` | repository root, paths from `config/paths*.yml` and `OPMD_*` overrides, `TAB()`, `FIG()`, `require_files()` |
| `R/lib/sample_sheet.R` | parses the GenomeStudio sample sheet of the pilot chip |
| `R/lib/utils.R` | small helpers shared by the pilot scripts |
| `R/lib/make_paths.R` | writes the configured output paths for the Makefile |

## Tools

| File | Use |
|---|---|
| `tools/compare_tables.R` | lists the tables (and columns) that differ between two output folders |
| `tools/build_manuscript_docx.py` | formats a pandoc Markdown manuscript as a Nature Communications-style .docx |
