# Where each figure and table comes from

The figure scripts only draw. Every value they plot is computed by the analysis script(s)
listed here and saved in `outputs/extended/`.

| Item | Content | Computed in | Drawn in |
|---|---|---|---|
| Fig. 1a-c | pilot PCA, PC1 vs epithelial fraction, EpiDISH fractions | `02_pilot_methylation/01_preprocess_sesame.R`, `03_methylation_cohorts/01_composition_signature_scores.R` | `07_figures/01_figures_1_to_4.R` |
| Fig. 1d-f | tumour-likeness (GSE204943 and pilot), composition independence, pilot ranking | `03_methylation_cohorts/01_composition_signature_scores.R` | same |
| Fig. 2a-d | KZFP promoter methylation vs PMD index, composition, groups, pilot deltas | `03_methylation_cohorts/01_composition_signature_scores.R`, `04_gse234379_kzfp.R` | same |
| Fig. 2e | GSE234379 adjacent-OSCC pairs | `03_methylation_cohorts/04_gse234379_kzfp.R` | same |
| Fig. 2f | KZFP methylation-expression correlation (GSE85195) | `03_methylation_cohorts/03_mumbai_crosscheck.R`, `06_audit/01_results_audit.R` (within-group) | same |
| Fig. 3a-b | conumee vs array CGH; aCGH arms by OPL group | `03_methylation_cohorts/02_opl_groups.R`, `03_mumbai_crosscheck.R` | same |
| Fig. 3c-e | pilot arm copy number, 9p21 bins, CDKN2A methylation | `02_pilot_methylation/02_identity_sex_copy_number.R` | same |
| Fig. 4b-d | KM curves, Cox HRs, random-set and permutation nulls | `04_expression_outcome/02_gse26549_transformation.R`, `06_audit/01_results_audit.R` (exact CIs) | same |
| Fig. 4e-f | integrated pilot profile; model | stage 02-04 outputs (panel f is a schematic) | same |
| Fig. 5a-c | RNA composition and programme scores vs EpiDISH | `05_pilot_rnaseq/03_rna_methylation_validation.R`, `04_rna_analyses.R` | `07_figures/02_figures_5_6.R` |
| Fig. 5d-e | arm dosage vs array copy number; 9p21 expression | `05_pilot_rnaseq/03_rna_methylation_validation.R` | same |
| Fig. 5f | cross-assay detection of copy-number events | RNA calls computed in the figure script from `rnaseq_arm_dosage_vs_cn_lpo.csv`; array, exome and 9p21 calls from `resources/cnv_events_curated.csv` | same |
| Fig. 6a-b | methylome-state proxy and programme shifts in the pilot | `05_pilot_rnaseq/03_rna_methylation_validation.R`, `04_rna_analyses.R` | same |
| Fig. 6c-e | hallmarks, silencing of hypermethylated promoters, KZFP enrichment | `05_pilot_rnaseq/04_rna_analyses.R` | same |
| Fig. 6f | outcome specificity: proxy vs lesion programme | `04_expression_outcome/02_gse26549_transformation.R`, `05_pilot_rnaseq/04_rna_analyses.R`, `06_audit/01_results_audit.R` | same |
| Supplementary Fig. 1 | score and proxy robustness | `06_audit/01_results_audit.R`, `04_expression_outcome/0*.R`, `03_methylation_cohorts/02_opl_groups.R` | `07_figures/01_figures_1_to_4.R` |
| Supplementary Fig. 2 | pilot genome-wide copy number, genotype concordance, 9p21 probe failure | `02_pilot_methylation/02_identity_sex_copy_number.R` | same |
| Supplementary Fig. 3 | RNA QC, strandedness, sex, paired DE, muscle content, RNA-only dosage | `05_pilot_rnaseq/00_strandedness_check.sh`, `01`-`04` | `07_figures/02_figures_5_6.R` |
| Supplementary Table 1 | public datasets | text (checked in `06_audit/02_cohort_verification.R`) | - |
| Supplementary Table 2 | GSE204943 leukoplakias: histology, outcome, methylome group | `06_audit/02_cohort_verification.R` (`verification_mumbai_opl_outcome.csv`) | - |
| Supplementary Table 3 | gene sets for RNA programme and cell-type scores | defined in `05_pilot_rnaseq/04_rna_analyses.R` (`cells`) and `03_rna_methylation_validation.R` | - |

Every number quoted in the text is listed, in reading order and formatted as in the text, in
`outputs/extended/tables/manuscript_numbers.csv` (`06_audit/03_manuscript_numbers.R`). The key
claims are also checked against their reported values in `audit_checks.csv`
(`06_audit/01_results_audit.R`).
