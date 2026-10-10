# Changelog

## 2026-10 Repository assembled from the analysis folder

### Script mapping

| Earlier file (`analysis/`) | Now |
|---|---|
| `methylation_report.Rmd` + `render.R` | `scripts/02_pilot_methylation/01_preprocess_sesame.R`, `03_exploratory_dmp_dmr.R` |
| `R/01_sample_sheet.R`, `R/utils.R` | `R/lib/sample_sheet.R`, `R/lib/utils.R` |
| `R/10_sigdf_qc_cnv.R` | `scripts/02_pilot_methylation/02_identity_sex_copy_number.R` |
| `R/11_external_methylation.R` | `scripts/03_methylation_cohorts/01_composition_signature_scores.R` |
| `R/12_external_rnaseq.R` | `scripts/04_expression_outcome/01_gse227919_programmes.R` |
| `R/13_story_analyses.R` | `scripts/03_methylation_cohorts/02_opl_groups.R` |
| `R/14_mumbai_cohort_crosscheck.R` | `scripts/03_methylation_cohorts/03_mumbai_crosscheck.R` |
| `R/15_gse26549_progression.R` | `scripts/04_expression_outcome/02_gse26549_transformation.R` |
| `R/16_gse234379_and_kzfp.R` | `scripts/03_methylation_cohorts/04_gse234379_kzfp.R` |
| `R/17_story_figures.R` | archived (superseded) |
| `R/18_results_audit.R` | `scripts/06_audit/01_results_audit.R` |
| `R/19_nature_figures.R` | `scripts/07_figures/01_figures_1_to_4.R` |
| `R/20_cohort_verification.R` | `scripts/06_audit/02_cohort_verification.R` |
| `R/21_rnaseq_matrices.R` | `scripts/05_pilot_rnaseq/02_gene_matrices.R` |
| `R/22_rnaseq_methylation_integration.R` | `scripts/05_pilot_rnaseq/03_rna_methylation_validation.R` |
| `R/23_rnaseq_analyses.R` | `scripts/05_pilot_rnaseq/04_rna_analyses.R` |
| `R/24_rnaseq_figures.R` | `scripts/07_figures/02_figures_5_6.R` |
| `rnaseq/run_kallisto_sequential.sh` | `scripts/05_pilot_rnaseq/01_kallisto_sequential.sh` |
| earlier GSE204943 project, `scripts/{01,07,10,16,25,27}` | `scripts/01_gse204943/01-05` (new; originals in `archive/`) |

The ported scripts differ from the originals only in how they find files (`R/lib/config.R`)
and in the changes listed below.

### Changes that can affect numbers

1. **Pilot sesame preprocessing is applied once.** The report ran
   `openSesame(prefixes, func = NULL)` (which already applies the default `prep = "QCDPB"`),
   then `openSesame(sdfs, func = getBetas)` on the resulting SigDFs, which applies QCDPB a
   second time. So quality masking, channel inference, non-linear dye-bias correction,
   pOOBAH and noob were all applied twice. `01_preprocess_sesame.R` applies them once.
   `SESAME_PREP_PASSES=2` restores the earlier behaviour for comparison.
   Effect: pilot β values shift by a median of 0.020 (single-pass values are 0.012 lower on
   average; per-sample r = 0.998), and 1,743 more CpGs pass filtering (762,517). Numbers that
   depend on pilot β values change, without changing any conclusion. For example:
   - PC1 vs epithelial fraction: r = 0.90 becomes 0.84;
   - P4 / P5 shifts vs own adjacent: +0.96 / +0.61 become +0.89 / +0.56;
   - CDKN2A promoter β in the three 9p21-deleted samples: up to 0.04 becomes up to 0.08,
     against 0.05-0.06 in adjacent tissue.

2. **GSE204943 probe filtering.** In the earlier pipeline
   (`archive/gse204943_original_scripts/01_data_processing_QC.R`), the per-probe keep vector
   (866,091 probes) was used to subset the RGChannelSet, whose rows are bead addresses
   (1,051,815). R recycles a logical index that is too short, so 11,593 address rows were
   dropped in a pattern unrelated to detection P. As a result, the frozen `beta.rds`
   (855,005 CpGs) retains 9,975 probes that failed the detection filter and lacks 10,971
   probes that passed. `01_minfi_preprocess.R` filters the noob MethylSet, whose rows are
   probes. Noob is also estimated on the unfiltered array, as intended.
   Effect: small. Examples:
   - PMD group difference adjusted for composition: P = 3.1e-4 becomes 3.6e-4;
   - conumee vs array CGH median per-lesion r: 0.92 / 0.67 becomes 0.91 / 0.66;
   - delta BIC: 17 becomes 18.
   The OPL groups (12 of 22 tumour-like), the transformation analyses and every outcome
   statistic are unchanged.

3. **Pilot copy-number segmentation is seeded.** DNAcopy `segment()` assesses change points
   by permutation and the earlier script set no seed, so repeated runs gave slightly different
   segments (e.g. 406 vs 405 leave-patient-out segments). `02_identity_sex_copy_number.R` now
   calls `set.seed(1)` before segmenting.

4. **The earlier pilot betas came from an April 2026 knitr cache.** `methylation_report.Rmd`
   cached its normalisation chunk (`cache = TRUE`), so the October `processed.rds` reused β
   values computed on 28 April. Those β values are identical to a fresh computation today
   wherever both are measured, but 46-48 probes per sample (replicate `_TC11`/`_BC11` probes)
   were masked in April and are not masked today: they pass today's quality mask and pOOBAH
   (P <= 0.015). After filtering and collapsing, today's matrix has 760,774 CpGs, against
   760,650 in April. This is a change in the sesame environment, not in the code.

5. **GSE204943 objects can be rebuilt from IDATs** (stage 01). The earlier objects came from
   interactive sessions. One script used an undefined object (`tot_int_f` in the copy-number
   script), and the per-sample PMD erosion score was defined in a different script from the
   blocks. Stage 01 reimplements each step with the original parameters.

### Changes that do not affect numbers

- Paths come from `config/paths.yml`; there are no hard-coded locations.
- Values that were typed into figure scripts are now read from the analysis outputs:
  - the Fig. 6f confidence interval (from the audit table);
  - the Fig. 2f group sizes (from the audit);
  - the Supplementary Fig. 3b strandedness percentages (from `00_strandedness_check.sh`,
    whose first run is recorded in `outputs/rnaseq/strandedness_check.tsv`).
- Curated copy-number calls for Fig. 5f moved to `resources/cnv_events_curated.csv`.
- Figure device: Quartz PDF on macOS, `cairo_pdf` elsewhere.

## Validation of the repository

All runs used the same software (R 4.5.2, Bioconductor 3.22; `env/renv.lock`).
`tools/compare_tables.R` compared the tables (relative tolerance 1e-6), and PNG figures were
compared pixel by pixel.

| Run | Pilot preprocessing | GSE204943 objects | Purpose | Result |
|---|---|---|---|---|
| A0 | earlier `processed.rds` and pilot copy-number objects, as inputs | frozen | Do the ported scripts reproduce the earlier analysis? | All 64 tables equal. 8 of 9 figures pixel-identical. In Fig. 6f the error bar ends moved by less than a pixel, because the exact CI (1.253-2.582) is now read from the audit table instead of a rounded value typed into the script. |
| A | legacy (`SESAME_PREP_PASSES=2`), today's environment | frozen | earlier method, recomputed today | 760,774 CpGs (April: 760,650; see change 4). The signature shares 3,999 of 4,000 CpGs with the earlier one. |
| B | single preprocessing (current) | frozen | effect of change 1 | 2 of 18 audited claims change beyond rounding (PC1 vs epithelial r; pilot shifts) |
| C | single preprocessing (current) | rebuilt from IDATs (stage 01) | effect of change 2; the configuration of record (`outputs/`) | 21 of 21 audited claims confirmed once the text was updated. Of the numbers quoted in the manuscript, 33 changed: 10 because of change 1 alone, 6 because of change 2 alone and 17 because of both. None reverses a conclusion. |

Stage 01 rebuilt against the frozen objects:

| Output | Frozen | Rebuilt | Agreement |
|---|---|---|---|
| CpGs after filtering | 855,005 | 856,001 | 9,975 failing probes removed, 10,971 passing probes restored |
| EpiDISH fractions (118 samples) | | | r >= 0.9998, max abs. difference <= 0.009 |
| Hypomethylated PMD blocks | 781 (92.1 Mb) | 804 (100.3 Mb) | 96.1% of frozen bp recovered |
| Per-sample PMD erosion score | | | r = 0.9997 |
| ZNF471 CpGs | 13 (all promoter) | 13 (all promoter) | same CpGs and regions; group means within 0.003 |
| conumee2 segments (96 lesions) | 3,489 | 3,526 | chromosome-level mean log ratio r = 0.9985, median abs. difference 0.001 |

### Other corrections found while validating

- **Unscripted numbers now come from the pipeline.** The delta BIC (mixture models) and the
  Hartigan's dip test P value for the OPL scores had been computed outside the scripts. They
  are now computed in `06_audit/01_results_audit.R`, with a fixed seed for the simulated P value.
- **Numbers typed into figure scripts are now read from outputs:**
  - the strandedness percentages: `05_pilot_rnaseq/00_strandedness_check.sh` records the test;
  - the Fig. 6f confidence interval;
  - the Fig. 2f group sizes.
- **Rounding.** The multivariable proxy P value was printed as 4.5 x 10^-3; the value is
  0.00445, so it now reads 4.4 x 10^-3.
- **Figures:**
  - overlapping labels fixed in Figs 1b, 1e, 2f, 5a, 5b, 6a and 6e;
  - true minus signs, and no negative zero, in Figs 1f, 2d and 4e;
  - Fig. 1f distances shown with one decimal (P3 8.0 vs P4 8.3).
- **Legend.** Supplementary Fig. 1e now states that no interval was computed for the proxy.
- `06_audit/03_manuscript_numbers.R` writes every number quoted in the text to
  `manuscript_numbers.csv`, so the effect of a future change on the text can be read off
  with `tools/compare_tables.R`.

### Manuscript files updated (9 Oct 2026)

All changed numbers were replaced in `Methylation_RNAseq_sections_NatCommun.md` and in the
older `Methylation_sections_NatCommun.md`. Both .docx files were rebuilt, and the figures were
replaced with the repository versions. The previous text and figures are in
`analysis/output/figures_nature/archive_pre_repo_rerun_2026-10-09/`.

### Open points for the authors

- `ZNF667-AS1` is counted among the "19q13 KZFPs" by the symbol rule (prefix ZNF), but it is an
  antisense lncRNA, not a zinc-finger protein. It is one of the 12 silenced "KZFPs" in Fig. 6e.
- The P2 adjacent tissue carries a 12p loss call (log2 ratio -0.08) that the text does not
  mention. It was present in the earlier run too.
- In the unadjusted ranking of Fig. 1f, P3 (8.0) and P4 (8.3) are nearly tied for closest to
  the tumour. The text's statement still holds, but only narrowly.
