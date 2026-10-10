# Methods

This document describes, script by script, what the analysis code does: procedure, parameters,
software, statistics and outputs. It is the reference for the manuscript Methods and must be
kept in step with the code (see [Keeping this document current](#keeping-this-document-current)).
Every script in `scripts/` and the shared code in `R/lib/` is described here, and
`make methods-check` verifies that this stays true.

**Conventions**
- Tests are two-sided unless stated otherwise.
- "BH" is the Benjamini–Hochberg false discovery rate.
- "Spearman" and "Pearson" refer to correlation coefficients.
- Stochastic steps use the seeds given here.
- "Values that come from the data", such as a group boundary or a noise floor, are listed in the
  table at the end with their current values. A rerun that changes one of them fails
  `make methods-check` until this document is updated.

## Software

R 4.5.2 with Bioconductor 3.22. Every R package and its version is pinned in `env/renv.lock`.
The command-line tools are in `env/environment.yml`.

| Package / tool | Version | Used for |
|---|---|---|
| minfi | 1.56.0 | GSE204943 preprocessing, PMD blocks |
| IlluminaHumanMethylationEPICanno.ilm10b4.hg19 | 0.6.0 | EPIC v1 annotation (hg19) |
| conumee2 | 2.1.2 (GitHub) | GSE204943 copy number |
| sesame / sesameData | 1.28.1 / 1.28.0 | pilot EPIC v2 preprocessing, manifests, cytobands, transcripts (hg38) |
| IlluminaHumanMethylationEPICv2anno.20a1.hg38 | 1.0.0 | EPIC v2 annotation |
| EpiDISH | 2.26.0 | cell composition |
| limma | 3.66.0 | differential methylation, per-gene expression tests |
| DMRcate | 3.6.0 | exploratory pilot DMRs |
| DNAcopy | 1.84.0 | pilot copy-number segmentation |
| survival | 3.8.6 | Cox models |
| mclust, diptest | 6.1.2, 0.77.2 | OPL score distribution |
| kallisto | 0.52.0 | pilot RNA-seq quantification |
| tximport, edgeR, msigdbr | 1.38.2, 4.8.2, 26.1.0 | gene matrices, paired DE, hallmark sets |

## Data

- **Pilot DNA methylation.** Eight samples on one Infinium MethylationEPIC v2.0 BeadChip:
  - tumour and adjacent tissue from patient 2;
  - lesion and adjacent tissue from patients 3–5.
- **Pilot RNA-seq.** Eleven libraries, 2 × 101 bp over two lanes:
  - lesion or tumour and adjacent tissue from patients 2 and 4–7;
  - lesion only from patient 3.
- **Public cohorts** (Supplementary Table 1; `data/README.md`):
  - GSE204943: EPIC v1, Mumbai; 22 healthy mucosa, 22 OPL, 74 OSCC.
  - GSE85514: array CGH; same Mumbai patients.
  - GSE85195: expression arrays; same Mumbai patients.
  - GSE234379: 450K, São Paulo; 46 OSCC with adjacent tissue.
  - GSE26549: HuGene 1.0 ST; 86 leukoplakias with oral cancer outcome.
  - GSE227919: RNA-seq; 66 oral biopsies.

## Stage 00 – Download

### 00_download/download_geo.sh

- **What it downloads.** The series matrices, platform tables and supplementary files of the
  public cohorts, from GEO.
- **Checksums.** Each file is checked against the MD5 of the copy used for the manuscript. A
  mismatch is reported (GEO sometimes regenerates files) but does not stop the run.
- **GSE204943.** The RAW archive is extracted only if no IDATs are present. The GSE204943
  series matrix and file list are also fetched, for the integrity audit.

## Stage 01 – GSE204943 from the raw IDATs

All five scripts write to `<outputs>/gse204943/` (`crucial_objects/`, `results/` and
`intermediate/`).

### 01_gse204943/01_minfi_preprocess.R

1. **Reading.** All IDAT pairs in the IDAT folder are read with
   `minfi::read.metharray.exp(force = TRUE)`.
   - IDAT names are matched to the GEO sample sheet by array barcode, after removing the GEO
     `GSM…_` prefix.
   - An IDAT without a sample-sheet row stops the script.
2. **Detection P.** `detectionP` (methylated + unmethylated vs background). A probe "fails" in
   a sample when P > 0.01.
3. **Sample filter.** Samples with ≥ 1% failed probes are removed. In the current data none
   are.
4. **Normalisation.** `preprocessNoob` (offset 15, single-sample dye-bias correction) on the
   unfiltered array.
5. **Probe filter.** Probes failing in ≥ 1% of the remaining samples are removed. This is
   applied to the MethylSet, whose rows are probes. An earlier pipeline applied it to the
   bead-address rows of the RGChannelSet by mistake (`docs/CHANGELOG.md`).
6. **Outputs.**
   - β = M/(M + U) for all retained probes, sex chromosomes included. `ann_auto.rds` holds the
     hg19 annotation (ilm10b4) of the autosomal probes.
   - No SNP, cross-reactive or sex-chromosome filter is applied here. Downstream analyses use
     the CpGs shared with the pilot or the autosomal annotation, as stated per script.
7. **Alignment value.** The mean β of probes annotated exactly "ZNF471" in TSS200/TSS1500 is
   stored in the sample sheet. Stage 03 uses it to check that β columns and sample rows are
   aligned.

### 01_gse204943/02_cell_composition.R

EpiDISH robust partial correlation (`epidish(method = "RPC")`, default `maxit = 50`) against
`centEpiFibIC.m`. This reference has 716 CpGs for epithelial, fibroblast and immune cells;
those present in the β matrix are used. Negative coefficients are set to zero and the
fractions rescaled to sum to one.

### 01_gse204943/03_pmd_blocks.R

1. **Data.** β values (`ratioConvert`, `mapToGenome`, hg19), all chromosomes.
2. **Clustering.** `cpgCollapse` with defaults: clusters of CpGs ≤ 500 bp apart and ≤ 1,500 bp
   wide; block groups within 250 kb; cluster value = mean β.
3. **Blocks.** `blockFinder` on the open-sea clusters:
   - design `~ case`, where case = OPL or OSCC and control = healthy mucosa;
   - β scale, `cutoff = 0.10`, `smooth = TRUE` (locfit by cluster, 250-kb span);
   - blocks with a negative smoothed case − control difference are kept;
   - no minimum length.
   - Permutations (`PMD_PERMUTATIONS`, default 0) only produce block P values, which are not
     used. They are unseeded.
4. **Per-sample erosion score** (individual CpGs, all array CpGs). Used only as a consistency
   check of the PMD index in stage 03.
   - `delta_in_minus_out` = mean β inside blocks − mean β of all other CpGs.
   - A block is eroded when its mean β is < 0.50 and ≤ (the sample's outside mean − 0.05).
   - Load = total width of the eroded blocks.
   - Score = mean of z(load) and z(−delta), z-scored across all 118 samples.

### 01_gse204943/04_copy_number_conumee2.R

1. **Input.** Total intensity (methylated + unmethylated) after `preprocessIllumina` with its
   defaults: Illumina background correction and normalisation to the internal controls, as in
   the original analysis.
   - Probes are not detection-filtered here; only the sample filter applies.
2. **Annotation.** `CNV.create_anno(array_type = "EPIC", detail_regions)`, hg19, autosomes,
   with the defaults:
   - bins of ≥ 15 probes, 50 kb–5 Mb;
   - no excluded regions.
3. **Fit.** Each OPL and OSCC with `CNV.fit` against all healthy mucosae: multiple regression of
   the query's log2 intensity on each reference sample. A reference correlating ≥ 0.99 with the
   query is dropped. The log ratio is query − fitted.
4. **Binning.** `CNV.bin` takes the median log ratio per bin and removes a per-sample shift.
5. **Segmentation.** `CNV.segment` (DNAcopy) with the defaults:
   - α = 0.001;
   - 50,000 permutations;
   - minimum width 5 bins;
   - undo splits closer than 2.2 SD.
   The seed is fixed (1) for each sample.
6. **Output.** Segment table (log2 ratio = segment mean − shift) for every OPL and OSCC.

### 01_gse204943/05_locus_cpg_sets.R

ZNF471 CpG table:
- CpGs whose RefGene names include ZNF471.
- Region = the RefGene group of the CpG's first ZNF471 entry.
- Class:
  - Promoter = TSS1500, TSS200, 5'UTR or 1stExon, in any island context;
  - Gene Body = Body or ExonBnd;
  - 3'UTR;
  - Other.
- Reported per CpG: the mean β per group, the OSCC − healthy difference, and whether the OPL
  change has the same sign.

## Stage 02 – Pilot methylation (EPIC v2)

### 02_pilot_methylation/01_preprocess_sesame.R

1. **Sample sheet.** Parsed from the GenomeStudio sheet (`R/lib/sample_sheet.R`). Group comes
   from the sample name: Tumor, Lesion or Normal (adjacent).
2. **Preprocessing.** `readIDATpair(platform = "EPICv2")`, then `prepSesame(prep = "QCDPB")`
   **once**. The steps, in order:
   - Q: quality mask (recommended masks);
   - C: Infinium-I channel inference;
   - D: non-linear dye-bias correction;
   - P: pOOBAH detection masking (P > 0.05);
   - B: noob (offset 15).
   - `SESAME_PREP_PASSES=2` reproduces the earlier double application, for comparison only.
3. **β values.** `getBetas` (masked probes as NA).
4. **QC table.** Unmasked fraction, pOOBAH pass rate (raw signal), mean intensity, and sex
   (`inferSex`).
5. **Filters.**
   - Probes masked in more than 25% of samples are removed.
   - Only chr1–22 probes are kept (sesame EPIC v2 manifest, hg38).
   - Uncollapsed matrix (used by DMRcate): complete probes only.
   - Collapsed matrix: replicate probes are averaged to one value per CpG
     (`betasCollapseToPfx`), then complete CpGs only.
   - All probe types are retained (cg, plus a few ch, nv and rs probes).
6. **M values.** log2(β/(1 − β)), with β clamped to [10⁻⁶, 1 − 10⁻⁶].
7. **PCA.** `prcomp` (centred, unscaled) of the 5,000 most variable collapsed CpGs, on all eight
   samples. A second PCA uses the seven non-tumour samples.

### 02_pilot_methylation/02_identity_sex_copy_number.R

Raw-signal analyses of the pilot arrays.

**Sex**
- `inferSex` on QCDPB β values.
- log2 ratio of mean chrY to mean autosomal total intensity, from dye-bias-corrected and
  noob-corrected raw signal.
- Mean chrX β.
- These are reported only; no thresholds are applied.

**Detection.** pOOBAH P on raw signal. The rate is the fraction of probes with P < 0.05.

**Identity**
- Genotype calls on the rs SNP probes (65 on EPIC v2), from unmasked QCDPB β:
  - β < 0.3 → 0;
  - β > 0.7 → 2;
  - otherwise 1.
- Concordance = fraction of identical calls for each pair of samples.

**Copy number**
1. **Intensity.**
   - Total intensity after quality mask, dye-bias correction and noob, with masked probes
     dropped.
   - Probes: autosomal cg probes with total intensity > 500 in all eight samples.
   - Transform: log2, centred on each sample's median.
2. **Exclusions.**
   - Probes in the HOXA, HOXB, HOXC, HOXD and ZIC1/ZIC4 clusters (hg38 windows in the
     script), where intensity tracks methylation, are removed.
   - The first and last 5 Mb of each chromosome are flagged as terminal and not used for calls.
3. **References.**
   - **Leave-patient-out (LPO):** the median of the other patients' adjacent tissues.
   - **Paired:** the patient's own adjacent tissue.
   - The log ratio is the residual of a least-squares fit of the sample on the reference, with
     an intercept and the reference as the single predictor. This is a conumee-style fit;
     conumee itself regresses on every reference sample.
4. **Bins.** Consecutive probes are grouped into bins of 50 within each chromosome. Bin value =
   median log ratio. Each sample is re-centred on the median of its non-terminal bins.
5. **Segmentation.** DNAcopy `smooth.CNA` and `segment`:
   - α = 0.001, minimum width 3 bins, undo splits closer than 2 SD;
   - 10,000 permutations, seed 1.
6. **Arm calls.**
   - Arm value = median of the non-terminal bins on the arm; arms with ≥ 10 bins.
   - Gain if LPO ≥ 0.07, loss if LPO ≤ −0.07.
   - A call is "paired-confirmed" when the paired value has the same sign and |value| ≥ 0.07.
   - It is "field-shared" when the patient's adjacent tissue has the same call.
7. **9p21.**
   - Region: chr9:21.6–22.25 Mb (hg38). Value = median probe-level LPO log ratio.
   - Probe-failure ratio = fraction of region probes with pOOBAH P ≥ 0.05, divided by the
     genome-wide fraction.
8. **Gene and band summaries.** Median bin value over each gene's transcript range or cytoband.
   Regions under 1 Mb are widened to 1 Mb.
9. **Burden.** Fraction of a 2.88-Gb genome in segments that are:
   - at least 3 Mb long;
   - |log2| > 0.12;
   - not entirely terminal.
   Focal events: < 3 Mb, |log2| > 0.3, ≥ 4 bins.

### 02_pilot_methylation/03_exploratory_dmp_dmr.R

Not used in the manuscript.

- **Samples.** Paired lesion vs adjacent analysis in patients 3–5 (three pairs).
- **Model.** limma on M values, design `~ 0 + group + patient`, contrast lesion − normal:
  - `eBayes(robust = TRUE, trend = TRUE)`;
  - BH;
  - Δβ = difference of the mean β values.
- **Candidate DMRs.**
  - `DMRcate::cpg.annotate`: EPIC v2, uncollapsed M values, default remapping of
    cross-hybridising probes, plain `eBayes`, permissive FDR 0.5.
  - `dmrcate(lambda = 1000, C = 2)`. If that fails, it is retried with `pcutoff = 1e-3`.

## Stage 03 – Methylation cohorts

All four scripts work on the CpGs measured on both arrays. These are the pilot EPIC v2 CpGs
after collapsing, intersected with the GSE204943 EPIC v1 CpGs (the "shared CpGs"). The
exceptions are stated where they apply.

### 03_methylation_cohorts/01_composition_signature_scores.R

**Purpose.** Place the eight pilot samples on an axis defined in GSE204943 that is adjusted for
composition.

**Inputs and data used**
- Pilot β values and PCA from `processed.rds`.
- GSE204943 β values, sample sheet, annotation, cell fractions, PMD blocks and ZNF471 CpG
  table (stage 01).
- EPIC v2 annotation `IlluminaHumanMethylationEPICv2anno.20a1.hg38`.

**Procedure**

1. **Alignment check.** The mean β of ZNF471 TSS200/TSS1500 CpGs must correlate at r > 0.999
   with the value stored in the sample sheet. This confirms that β columns and sample rows are
   aligned.
2. **Platform concordance.** Pearson r between per-CpG mean β of the four pilot adjacent
   tissues and of the 22 GSE204943 healthy mucosae, over the shared CpGs.
3. **Pilot composition (EpiDISH, RPC).**
   - First tier: `epidish(β, centEpiFibIC.m, method = "RPC")` on the shared CpGs, giving
     epithelial, fibroblast and immune fractions. **These first-tier fractions are used
     everywhere downstream.**
   - Second tier: `hepidish` with `centBloodSub.m` splits the immune fraction into seven blood
     subtypes. It is reported in `deconvolution_epidish.csv` only.
   - PC1 of the pilot-only PCA is correlated (Pearson) with the epithelial fraction.
4. **Tumour signature in GSE204943.**
   - Model: limma on M values (β clamped to [10⁻⁶, 1 − 10⁻⁶]) with design
     `~ group + Epi + IC + plate`, fitted to healthy mucosa and OSCC only. Plate is a factor;
     a blank plate becomes "none". `eBayes` uses default settings.
   - Statistics: the moderated t for OSCC vs healthy and its BH FDR across all shared CpGs.
     The same model is fitted for OPL vs healthy and reported in the per-CpG table.
   - Δβ is the difference of the unadjusted group mean β values.
   - Signature:
     - **hypermethylated:** FDR < 10⁻⁶ and Δβ > 0.2, the 2,000 with the largest t;
     - **hypomethylated:** FDR < 10⁻⁶ and Δβ < −0.2, the 2,000 with the smallest t.
5. **Tumour-likeness score.**
   - Each signature CpG is scaled as (β − μ_healthy)/(μ_OSCC − μ_healthy), using the
     GSE204943 group means.
   - The medians over the hypermethylated and over the hypomethylated CpGs are averaged.
   - A score of 0 is the healthy-mucosa mean and 1 is the OSCC mean.
   - Scored samples: all GSE204943 samples and the pilot.
     - GSE204943 healthy mucosa and OSCC define the CpGs and anchors, so they are scored in
       sample.
     - OPLs and the pilot are scored out of sample.
6. **Pilot summaries.**
   - Each lesion's shift = score − score of the same patient's adjacent tissue.
   - Percentile of each pilot score among the 22 GSE204943 OPLs, from the empirical CDF.
7. **Purity-scaled score (supplementary).**
   - A loess fit of score on epithelial fraction in GSE204943 OSCC (span 0.9, degree 2) gives
     the expected tumour score at a given epithelial content.
   - Each score is divided by that prediction, floored at 0.15.
8. **Composition independence.** Spearman ρ of score with epithelial fraction within healthy
   mucosa, and across all GSE204943 samples.
9. **PCA projection.**
   - PCA (centred, unscaled β) of the 5,000 most variable shared CpGs in GSE204943.
   - Pilot samples are projected with `predict()`.
10. **Replication of the pilot changes in GSE204943.**
    - For each pilot pair, the lesion − adjacent Δβ at the signature CpGs is compared with the
      GSE204943 OSCC − healthy Δβ:
      - Pearson r;
      - fraction of CpGs with concordant sign and |Δβ| > 0.05;
      - fraction with discordant sign and |Δβ| > 0.05.
    - Background: the same number of random non-signature CpGs, with seed 1. The same random
      set is used for every patient.
11. **CpGs that change consistently in the three premalignant pairs (P3, P4, P5).**
    - Selection: all three Δβ > 0.1 (gain) or all three < −0.1 (loss).
    - Null: the same count under the other six sign patterns of the three pairs. These form
      three distinct values, because a pattern and its complete reversal give the same total.
    - Replication in GSE204943:
      - same direction as the pilot change, with nominal limma P < 0.05, for OPL vs healthy and
        for OSCC vs healthy;
      - compared with the same rate among all other shared CpGs;
      - no formal test.
12. **PMD hypomethylation index.**
    - CpG sets, both restricted to open-sea CpGs:
      - inside: open-sea CpGs within the GSE204943 hypomethylated blocks (stage 01);
      - outside: all other open-sea CpGs.
    - Index = mean β inside − mean β outside. More negative means more hypomethylation of PMDs.
    - The manuscript reports the index with the opposite sign (outside − inside), so that
      higher means more hypomethylation.
    - Agreement with the per-sample PMD erosion score from stage 01 is checked by Spearman ρ.
      It is negative by construction, because the two measures run in opposite directions.
13. **ZNF471 promoter methylation.** Mean β of the ZNF471 CpGs classed "Promoter" in stage 01:
    TSS1500, TSS200, 5'UTR or 1stExon, in any island context.
14. **Candidate-gene promoter table.**
    - Gene list: 34 genes covering copy-number and tumour-suppressor genes, frequently mutated
      genes, the spatial-profiling programme genes and the GSE204943/pilot programme. Each is
      listed in the script.
    - Promoter CpGs (EPIC v2 annotation):
      - TSS200, TSS1500, 1stExon or 5'UTR in an island or shore;
      - if a gene has none, its TSS200/TSS1500 CpGs in any context.
    - Reported per gene: group means, minimum nominal limma P, and pilot pair Δβ.

**Main outputs.** `gse_signature.rds`, `external_methylation.rds` and tables
`tumor_likeness_*.csv`, `deconvolution_epidish.csv`, `pmd_znf471_pilot.csv`,
`replication_vs_gse204943.csv`, `pilot_consistent_cpgs_*.csv` and
`candidate_gene_promoter_methylation.csv`.

### 03_methylation_cohorts/02_opl_groups.R

**Purpose.** Define the two OPL groups and test their association with copy number. Also build
arm-level expression dosage in GSE227919.

**GSE204943 copy number by arm**
- Arms: from the hg19 UCSC cytobands.
- Arm value: the conumee2 segment log2 ratios (stage 01) averaged over each arm, weighted by
  overlap length.
- Each sample is then centred on its median arm value.
- Derived features:
  - core loss = −(3p + 9p)/2;
  - number of the seven oral LOH arms (3p, 9p, 13q, 4q, 8p, 11q, 17p) with log2 ratio ≤ −0.1.

**OPL groups**
- Density: Gaussian kernel density estimate of the 22 OPL tumour-likeness scores, bandwidth
  0.15 (`density()` defaults otherwise).
- Boundary: the density minimum between 0.1 and 0.8.
- OPLs above the boundary are **tumour-like**, the rest **normal-like**. The boundary value is
  in the data-derived values table.
- The scores used here are read from the stage-03 table, rounded to three decimals.

**Tests within the OPLs**
- Arm features vs tumour-likeness:
  - Spearman ρ, without and with adjustment for epithelial fraction (rank-based partial
    correlation: Pearson r of the residuals of rank(x) and rank(y) on rank(Epi)).
  - No correction for multiple testing.
- Group comparisons (two-sided Wilcoxon rank-sum, normal approximation) of:
  - the PMD index;
  - ZNF471 promoter methylation;
  - epithelial fraction;
  - LOH-arm count and core loss.

**GSE227919 expression-inferred arm dosage**
- Input transform: the deposited matrix is counts per million (each sample sums to 10⁶),
  analysed as log₂(CPM + 1).
- Gene filter: mean > 3 and > 1 in ≥ 90% of samples.
- Gene position: the hg38 transcript start (minimum over transcripts). The HLA region
  (chr6:28–34 Mb) is excluded.
- Gene z scores relative to the control biopsies: (x − control median) / max(control MAD,
  0.1), clipped to ±4.
- Arms: those with ≥ 40 genes are used.
- Arm dosage: mean z per arm, minus the sample's mean over arms.
- Validation, cancer vs control:
  - two-sided Wilcoxon per arm, BH across arms;
  - expected pattern: loss of 3p, 4q, 8p, 13q and 17p; gain of 3q, 8q and 20.
- Progression, among non-cancer biopsies with known outcome (progressed = any status other
  than "Stable"):
  - AUC (Mann–Whitney);
  - two-sided Wilcoxon;
  - logistic regression adjusted for histological class.

### 03_methylation_cohorts/03_mumbai_crosscheck.R

**Purpose.** Same-patient checks in the Mumbai cohort: GSE204943 methylation, GSE85514 array
CGH and GSE85195 expression.

**Sample matching.** Samples are matched by patient identifier: the sample title after the last
underscore, with spaces removed.

**Array CGH**
- Probes: GPL22311 probes with chromosome coordinates, each assigned to an arm.
- Each sample is median-centred over probes.
- Arm value = median probe log ratio.
- A focal 9p21 value is the median over chr9:21.75–22.15 Mb.

**Cross-platform agreement (conumee vs array CGH)**
- Pooled Pearson r over all sample × arm pairs.
- Per-lesion Pearson r across arms; the median is reported by group.
- Null: aCGH arm values permuted across samples within each arm, 500 iterations with seeds
  1–500. The median per-lesion r of each iteration is kept.
- Direction agreement: the fraction of altered arms (|aCGH| > 0.1) with the same sign.

**Array CGH by OPL group.** For the seven LOH arms, focal 9p21, core loss and instability (the
fraction of arms with |log2 ratio| > 0.1):
- two-sided Wilcoxon between the groups;
- Spearman ρ with tumour-likeness.
- BH adjustment across the seven arms is done in the audit (stage 06).

**Methylation vs expression (GSE85195)**
- Probes: one probe per gene, the one with the highest mean.
- Promoter methylation: mean β of island CpGs annotated TSS200, TSS1500 or 1stExon (EPIC v1
  annotation, all GSE204943 CpGs).
- Per-gene Spearman ρ across the matched samples, for genes whose promoter methylation varies
  (SD > 0.05).
- Signature-hypermethylated genes are compared with the other genes by a two-sided Wilcoxon
  test.

**Expression programmes by OPL group.**
- Genes are z-scored across all GSE85195 samples.
- Programme score = mean z.
- The groups are compared by a two-sided Wilcoxon test.

### 03_methylation_cohorts/04_gse234379_kzfp.R

**Purpose.** Transfer the signature, PMD index and KZFP promoter methylation to GSE234379, and
summarise KZFP methylation in GSE204943 and the pilot.

**KZFP promoter CpGs (this script and the audit).** Island CpGs annotated TSS200, TSS1500 or
1stExon to ZNF471, ZNF382, ZNF529, ZNF582 or ZNF671 (EPIC v1 annotation). The KZFP value is the
mean β of each gene's CpGs; the composite is the mean over the five genes.

**GSE234379 (450K)**
- Input: the deposited normalised β values.
- CpG filter: CpGs with < 10% missing values.
- Pairs: patients with exactly two samples.
- Scores computed with the unchanged GSE204943 anchors, on the signature CpGs present on the
  450K:
  - tumour-likeness;
  - the PMD index;
  - KZFP methylation.
- Composition: EpiDISH (RPC, centEpiFibIC) on complete-case CpGs.
- Statistics:
  - AUC of tumour vs adjacent tissue, over all samples;
  - two-sided Wilcoxon signed-rank test on the within-patient differences;
  - fraction of pairs changing in the expected direction;
  - adjacent tissue scoring above the highest GSE204943 healthy-mucosa score, or above 0.35.
- HPV status is parsed here. The HPV-negative subset is analysed in
  `06_audit/02_cohort_verification.R`.

**GSE204943 and pilot**
- GSE204943: median KZFP β per group (healthy, normal-like OPL, tumour-like OPL, OSCC).
  Normal-like and tumour-like OPLs are compared by a two-sided Wilcoxon test.
- Pilot: KZFP promoter Δβ of each lesion vs its own adjacent tissue.

## Stage 04 – Expression cohorts and outcome

### 04_expression_outcome/01_gse227919_programmes.R

**Purpose.** Score transcriptional programmes and methylation-linked gene sets in GSE227919
(66 oral biopsies) against histological grade and progression.

**Expression**
- Transform: log₂(CPM + 1).
- Gene filter: > 1 in ≥ 20% of samples.
- Each gene is z-scored across all samples.

**Gene sets**
- The seven spatial-profiling programmes and PI16 fibroblast. The genes are listed in the
  script.
- KZFP_chr19_hyper: ZNF382, ZNF529, ZNF471.
- p16_9p21: CDKN2A, CDKN2B, MTAP.
- Two sets from the pilot and GSE204943:
  - **Methylation-silenced promoters:** genes with promoter island CpGs (TSS200, TSS1500 or
    1stExon) among the signature-hypermethylated CpGs that gained > 0.1 in at least two of
    the three pilot pairs P2, P4 and P5.
  - **Methylation-activated promoters:** the analogue for hypomethylated CpGs (Δβ < −0.1),
    with promoter CpGs in any island context.
- Score = mean z of the set genes present.

**Tests**
- Trend with histological class (control < hyperkeratosis < dysplasia < cancer, coded 1–4):
  Spearman ρ, BH across sets.
- Progression, among non-cancer biopsies with known outcome (progressed = any status other
  than "Stable"):
  - AUC with a bootstrap 95% percentile interval (2,000 resamples, seed 7, so every set uses
    the same resamples);
  - two-sided Wilcoxon, BH across sets;
  - logistic regression on score + histological class (odds ratio per SD, Wald P).
- Reference: the AUC of histological class alone.
- A per-gene limma-trend table (dysplasia and cancer vs control) is written for the programme
  genes.

### 04_expression_outcome/02_gse26549_transformation.R

**Purpose.** Test whether the methylation-defined state, through an expression proxy, predicts
oral cancer in GSE26549 (86 leukoplakias from a chemoprevention trial).

**Data**
- Expression: the deposited series-matrix values, which are already on the log2 scale.
- Probe sets are mapped to symbols with `hugene10sttranscriptcluster.db`, keeping the probe set
  with the highest mean per gene. Genes are z-scored across lesions.
- Outcome: oral cancer-free survival (`ocfs_years`); event = "oral cancer development".

**Features**
- The eleven GSE227919 gene-set scores.
- Expression-inferred arm loss, built the same way as in GSE227919 with these differences:
  - genes above the 30th percentile of mean expression;
  - z scores clipped to ±4;
  - arm means centred on the sample's mean over arms;
  - loss = −dosage.
  Features: 3p, 9p, 13q, 4q, core 3p + 9p, and the mean of the four arms.
- **Methylome-state proxy:**
  - Training: the 15 Mumbai OPLs with both GSE204943 methylation and GSE85195 expression.
  - Training genes: those also measured in GSE26549 with SD > 0.3 in the training lesions.
  - Each training gene is correlated (Spearman) with methylation tumour-likeness.
  - The 50 most positively and 50 most negatively correlated genes form the proxy.
  - Proxy score in GSE26549 = mean z(up) − mean z(down). It is applied without refitting.

**Cox models (`survival`)**
- Per feature, per SD:
  - unadjusted;
  - adjusted for baseline histology (dysplasia vs hyperplasia).
- Likelihood-ratio test for adding the feature to histology.
- Harrell's C, with and without histology.
- BH across the 18 pre-specified features.
- Kaplan–Meier curves for a median split of the proxy.

**Robustness of the proxy**
1. **Multivariable model:** adjusted for histology, ever smoking, age, trial arm and biopsy
   timing.
2. **Specificity:** the proxy's correlation with the programme scores and arm losses.
   Adjustment for the cycling-epithelial and dysplastic-epithelial programmes.
3. **Composition:** adjustment for stromal, immune and keratinocyte marker scores (genes listed
   in the script), with histology. Also the likelihood-ratio test of the proxy over these.
4. **Size sensitivity:** proxies of 10, 25, 50, 100 and 200 genes per direction, and the
   up-only and down-only halves.
5. **Random-signature null:** 1,000 random, non-overlapping 50 + 50 gene sets drawn from the
   training genes with a defined correlation (seed 11).
6. **Label-permutation null:** 1,000 iterations. In each, the training tumour-likeness labels
   are shuffled once, every gene's correlation is recomputed against the shuffled labels, and
   the 50 + 50 proxy is rebuilt and tested.

   For both nulls, empirical P = (number of null Cox P ≤ observed P + 1) / 1,001.
7. **KZFP expression:** Cox per SD of ZNF471, ZNF382, ZNF529, ZNF582, ZNF154 and ZNF671
   expression. In the text, HRs are given per SD of *lower* expression, i.e. 1/HR.
8. **Replication in GSE227919:** the same proxy genes, scored as above. Tested for trend with
   histological class and for progression (AUC, Wilcoxon, logistic regression adjusted for
   class).

Further checks are in the audit: proportional hazards, bootstrap and leave-one-out (stage 06).

## Stage 05 – Pilot bulk RNA-seq

### 05_pilot_rnaseq/00_strandedness_check.sh

The first 300,000 read pairs of library S1 (lane 1) are pseudoaligned with kallisto in
unstranded, `--fr-stranded` and `--rf-stranded` modes. The library was dUTP (reverse-stranded):
the `--rf-stranded` rate (77.0%) was close to the unstranded rate (78.2%), and `--fr-stranded`
was 3.6%. The recorded result is `strandedness_check.tsv`, shown in Supplementary Fig. 3b.

### 05_pilot_rnaseq/01_kallisto_sequential.sh

- **Index:** kallisto 0.52.0 index of the GENCODE v50 transcriptome (k = 31, the default).
- **Quantification:** `kallisto quant --rf-stranded`, with both lanes of each library as one
  paired-end run, 8 threads, default bias and fragment-length estimation, no bootstraps.
- **Order and storage:** libraries are processed one at a time. When downloading, each FASTQ is
  verified by size and MD5 before use and deleted after quantification.

### 05_pilot_rnaseq/02_gene_matrices.R

- **Transcript to gene:** the GENCODE header fields: gene ID (field 2), gene name (field 6) and
  gene type (last field).
- **Summarisation:** tximport (kallisto) to gene level. It writes:
  - estimated counts;
  - length-scaled-TPM counts (`countsFromAbundance = "lengthScaledTPM"`);
  - TPM;
  - mean transcript length.
- **Sample sheet:** patient, tissue, kallisto QC, and the matching methylation sample.

### 05_pilot_rnaseq/03_rna_methylation_validation.R

**Purpose.** Check the methylation readouts against RNA from the same specimens. The analysis
is descriptive.

**Expression.** log₂(TPM + 1). For genes with duplicate symbols, the gene ID with the highest
mean TPM is kept.

**Checks**
- **Sex:** XIST and Y-linked genes.
- **Composition:**
  - Marker scores (mean z across libraries, all genes) for epithelial, immune, lymphoid,
    fibroblast and endothelial markers (listed in the script).
  - Compared with EpiDISH by Spearman ρ. The non-epithelial score is −(immune + fibroblast).
- **9p21:** log2 fold change of MTAP, CDKN2A, CDKN2B, DMRTA1 and IFNB1 in each lesion vs its own
  adjacent tissue.
- **Arm dosage vs array copy number:**
  - Genes with TPM > 3 in every library; HLA region excluded; arms with ≥ 30 genes.
  - Paired version: lesion − own adjacent, arm mean, centred on the mean over arms.
  - Leave-patient-out version: each library against the median of the other patients' adjacent
    tissues.
  - Compared with array copy number by per-sample Pearson r across arms.
- **KZFP methylation vs expression:** Spearman ρ across gene × lesion pairs.
- **Methylome-state proxy in the pilot:**
  - The GSE26549 proxy genes, scored as mean z(up) − mean z(down).
  - Correlated (Spearman) with methylation tumour-likeness.
  - Also a rank regression of proxy on tumour-likeness adjusted for the immune score.

### 05_pilot_rnaseq/04_rna_analyses.R

**Purpose.** The RNA analyses reported in the manuscript.

**Expression and pairs**
- log₂(TPM + 1). z scores across libraries for genes with TPM ≥ 1 in ≥ 2 libraries.
- Five lesion/tumour–adjacent pairs: P2, P4, P5, P6, P7. P3 has no adjacent library.

**Sections**
- **A. Programme and cell-type scores.**
  - Mean z of the seven spatial-profiling programmes and 13 cell-type and process sets (listed
    in the script and in Supplementary Table 3).
  - Lesion − adjacent change for each pair.
- **B. Paired differential expression.**
  - Data: edgeR on length-scaled-TPM counts, with `filterByExpr` (design-aware, default
    thresholds) and TMM normalisation.
  - Model: design `~ patient + lesion`; `estimateDisp`; `glmQLFit(robust = TRUE)`;
    quasi-likelihood F-test for the lesion term; BH FDR.
  - Hallmark gene sets (msigdbr, collection H) are tested with `camera` on the same design.
- **C. Promoter methylation vs expression.**
  - Promoter: island CpGs annotated TSS200, TSS1500 or 1stExon.
  - Gene promoter β = mean over its CpGs.
  - For P2, P4 and P5, each lesion is compared with its own adjacent tissue, using genes with
    TPM ≥ 1 in either sample.
  - Classes: gain (Δβ ≥ 0.2), stable (|Δβ| < 0.05) and loss (Δβ ≤ −0.2).
  - The expression change of each class is compared with stable by a two-sided Wilcoxon test.
- **D. RNA arm dosage for all 11 libraries (leave-patient-out).**
  - Noise floor = the 99th percentile of |dosage| over all arms of the adjacent libraries.
  - An array-called event counts as detected in RNA when the arm dosage passes the floor in the
    expected direction. For 9p21, detection is at least twofold lower expression of all three
    genes.
- **E. The generic lesion programme.**
  - Genes: the 50 protein-coding genes most significantly up and the 50 most significantly down
    in section B.
  - GSE26549: scored as up − down mean z and tested by Cox (unadjusted and histology-adjusted).
  - Null: 1,000 random gene sets from the protein-coding genes tested in B. Empirical
    P = (k + 1)/1,001.
  - GSE227919: tested for trend with grade and for progression (AUC; two-sided Wilcoxon,
    computed in the figure and number scripts).
- **F. Robustness.**
  - Composition-tracking genes: genes whose expression correlates (Spearman ρ > 0.5) across
    libraries with the immune, fibroblast (CAF) or endothelial score. They are removed and
    section C is repeated.
  - Skeletal-muscle score: mean log₂(TPM + 1) of ACTA1, MYH2, MYH7, CKM and DES.
- **G. 19q13 KZFP enrichment.**
  - KZFP definition: genes on chromosome 19 with symbols starting ZNF, ZFP, ZSCAN, ZIK, ZKSCAN
    or ZIM. This is a name rule, not a KRAB-domain annotation.
  - Universe: CpG-island promoters measured in both premalignant lesions (P4 and P5).
  - Tests: two-sided Fisher's exact test for KZFPs among promoters gaining ≥ 0.2 in both
    lesions, and among those that also have lower expression in both.

## Stage 06 – Audit

### 06_audit/01_results_audit.R

**Purpose.** Re-derive the headline results and run the checks that could overturn them. Each
check is written to `audit_checks.csv` with the reported value, the recomputed value and a
status (`confirmed` / `CHECK`).

**A. Circularity**
- The score is recomputed after removing two sets of CpGs from the signature:
  - every chromosome 19 CpG annotated to a gene named ZNF*, ZFP*, ZSCAN* or ZIK*;
  - every open-sea CpG, which includes all open-sea CpGs in PMD blocks (the CpGs the PMD index
    uses).
- The OPLs are re-classified with the same boundary and compared with the original groups.
- KZFP and PMD differences between the new groups are tested (two-sided Wilcoxon).

**B. Composition**
- In healthy mucosa: Spearman ρ of KZFP methylation and of the PMD index with epithelial
  fraction.
- In OPLs: linear models of KZFP methylation, and of the PMD index, on group + epithelial +
  immune fraction.
- Rank-based partial correlation of KZFP methylation with PMD hypomethylation, given
  composition.

**C. Multiplicity.** BH across the seven array-CGH arms, and across the KZFP genes.

**D. Outcome model (GSE26549)**
- Proxy per SD: Wald 95% CI.
- Proportional hazards: `cox.zph`.
- Histology-adjusted HR.
- Bootstrap (1,000 resamples, seed 42) of Harrell's C and of the HR, with percentile
  intervals.
- Leave-one-out refits.
- Unadjusted and histology-adjusted HRs with exact CIs for every feature (the forest table).

**E. Methylation vs expression within groups.** Spearman ρ of KZFP promoter methylation with
GSE85195 expression of the same gene, in all matched samples, within OSCC and within OPL.

**F. Numbers and the shape of the OPL score distribution**
- Gap between the two groups.
- Gaussian mixture models with 1–3 components (mclust, unequal variance), compared by BIC.
- Hartigan's dip test with a simulated P value (20,000 simulations, seed 20261009).

### 06_audit/02_cohort_verification.R

- **Mumbai OPL outcome by group.**
  - Histology and status at last follow-up are transcribed from the published supplementary
    tables.
  - Two-sided Fisher's exact tests for transformation, and for transformation or recurrence.
  - Supplementary Table 2.
- **Sex-adjusted signature.** The limma model with sex added: CpGs retained, OPL score
  correlation, OPL group agreement and change in pilot shifts.
- **GSE234379 by HPV status.** AUC, plus a signed-rank test of within-patient KZFP change in
  the HPV-negative and HPV-positive patients.
- **PMD blocks hypomethylated in OSCC alone.**
  - Blocks are kept if mean OSCC − healthy β < −0.10.
  - The index is recomputed on them and compared with the original index.
  - The OPL group difference is retested.
- **GSE227919 proxy AUC by subset.**
  - all biopsies;
  - hyperkeratosis and dysplasia only;
  - without three possible repeat patients;
  - cancer events vs stable.
- **KZFP CpG accounting.** CpGs per gene, overlap between ZNF382 and ZNF529, and the number of
  KZFP CpGs in the signature.

### 06_audit/03_manuscript_numbers.R

Writes every number quoted in the manuscript, formatted as in the text, to
`manuscript_numbers.csv` (section, statement, value). It reads stage 01–06 outputs and computes
a few summary statistics itself, such as correlations of saved scores and PCA distances. After
any change, compare this file before and after (`tools/compare_tables.R`).

### 06_audit/04_gse204943_integrity.R

Checks the stage-01 objects. Results go to `integrity_gse204943_summary.csv` with status pass /
note / FAIL / info.

- **Files:** local IDATs vs the GEO file list, and the RAW archive size.
- **Labels:** sample name, group and sex vs the GEO series matrix, joined by IDAT barcode.
- **Detection failures:** must be < 1% per array.
- **Intensity:** median log2 intensity; arrays below 10.5 are a note.
- **Sex from the data:**
  - two-cluster split of median log2 chrY − chrX intensity;
  - X CpG-island β > 0.3 for female.
  - A sample is discordant only if both signals disagree with the recorded sex.
- **Duplicated individuals:** concordance of genotype calls on the 59 SNP probes, flagged if
  > 0.85.
- **Profile vs group:** each sample's correlation with each group's median profile, over 50,000
  random CpGs (seed 1).
- **Robustness:** the signature re-derived without sex-discordant OSCCs, and without OSCCs
  scoring below 0.35. Pass if all OPLs keep their group and pilot shifts change by < 0.05.

## Stage 07 – Figures

### 07_figures/01_figures_1_to_4.R and 07_figures/02_figures_5_6.R

The figure scripts draw saved results. The computations made inside them are:

- **Fig. 1f:** Euclidean distance in the pilot PC1–PC2 plane from each lesion to the P2 tumour.
  Lesions are ranked by this distance and by composition-adjusted shift.
- **Fig. 2e:** paired Wilcoxon signed-rank test of KZFP methylation in GSE234379.
- **Fig. 4b:** Kaplan–Meier by median split of the proxy. The P value shown is the Cox Wald P
  for the continuous score.
- **Fig. 4c:** gene-level KZFP HRs, recomputed from the series matrix as Cox per SD of lower
  expression.
- **Fig. 4e:** each feature's value relative to the median of the four adjacent tissues,
  oriented towards tumour and scaled by its maximum.
- **Fig. 5f:** RNA calls for array copy-number events.
  - The arm events are computed with the noise floor.
  - The exome calls and the 9p21 RNA calls come from `resources/cnv_events_curated.csv`.
  - The exome copy-number calls are from the multimodal manuscript. The exome data are not
    analysed here.
- **Fig. 6f:** the GSE227919 progression P for the lesion programme (two-sided Wilcoxon).

**Format**
- Size: 183 mm wide, ≤ 170 mm tall.
- Text: Helvetica 4.5–7 pt.
- Files: Quartz PDF on macOS (`cairo_pdf` elsewhere) and 600-dpi PNG.
- **P values:**
  - two significant figures;
  - below 0.001 in scientific notation;
  - at or above 0.1, rounded to two decimals.

## Shared code

- **`R/lib/config.R`**
  - Finds the repository root and reads `config/paths.yml`, merged with the git-ignored
    `config/paths.local.yml`.
  - Applies `OPMD_<ENTRY>` environment overrides.
  - Defines every path used by the scripts, and the helpers `TAB()`, `FIG()` and
    `require_files()`.
- **`R/lib/sample_sheet.R`** parses the GenomeStudio sample sheet of the pilot chip:
  - the `[Data]` block;
  - patient, and group from the sample name;
  - IDAT basenames.
  It stops if an IDAT is missing.
- **`R/lib/utils.R`** provides `beta_to_m()`: β clamped to [10⁻⁶, 1 − 10⁻⁶], then
  log2(β/(1 − β)).
- **`R/lib/make_paths.R`** writes the configured output paths for the Makefile. It has no
  analysis role.

## Definitions that differ between analyses

These differences are deliberate or inherited from the original analysis. Each is stated here
so that a reader does not assume one definition throughout.

1. **Promoter CpGs.**
   - KZFP and gene-level promoter methylation (stages 03–06): island CpGs annotated TSS200,
     TSS1500 or 1stExon.
   - The ZNF471 index in `03_methylation_cohorts/01`: TSS1500, TSS200, 5'UTR or 1stExon, in any
     island context. The alignment check uses an exact "ZNF471" TSS200/TSS1500 match.
   - The candidate-gene table: island or shore, and also 5'UTR.
   - The GSE227919 methylation-activated set: any island context.
   - In all of these, the TSS/1stExon test is applied to the CpG's full RefGene group list. A
     CpG that is promoter-proximal for one gene is therefore also counted for any other gene
     named on that CpG.
2. **Group boundary.**
   - OPL groups use the data-derived boundary (KDE minimum).
   - Some sensitivity checks use the fixed value 0.35, which equals the current boundary to
     two decimals:
     - the sex-adjusted signature;
     - the GSE234379 field effect;
     - the integrity audit.
3. **Replication criteria.** Replication of pilot changes in GSE204943 uses nominal limma
   P < 0.05, not FDR.
4. **Annotation and genome build.**
   - GSE204943-anchored analyses use the EPIC v1 annotation, hg19.
   - Pilot-only annotation uses EPIC v2, hg38.
   - Pilot copy number and RNA gene positions are hg38; GSE204943 copy number is hg19.
5. **Copy-number methods.** GSE204943 uses conumee2: multiple regression on the reference
   samples, binning and DNAcopy. The pilot uses a single-reference conumee-style fit with
   50-probe bins. Arm values are summarised differently on each platform:
   - conumee: length-weighted segment mean;
   - array CGH: median of probes;
   - pilot: median of bins.
6. **Expression scaling.**
   - `05_pilot_rnaseq/03` z-scores all genes across libraries.
   - `05_pilot_rnaseq/04` z-scores genes with TPM ≥ 1 in ≥ 2 libraries.
   - Programme scores therefore differ slightly between the two scripts.
7. **KZFP gene definition.** In the RNA enrichment test, KZFPs are chromosome 19 genes named
   ZNF*, ZFP*, ZSCAN*, ZIK*, ZKSCAN* or ZIM*. The rule includes an antisense lncRNA
   (ZNF667-AS1).
8. **In-sample anchors.** GSE204943 healthy mucosa and OSCC both define the signature and are
   scored by it, so their own scores are in sample.

## Values that come from the data

These values appear in this document or in the manuscript Methods. Each is computed by the
pipeline and listed in `outputs/extended/tables/manuscript_numbers.csv` under the statement
shown. `make methods-check` fails when a value here no longer matches.

| Statement in `manuscript_numbers.csv` | Value |
|---|---|
| `pilot autosomal CpGs after filtering and collapsing (Methods)` | `762,517` |
| `GSE204943 CpGs and samples after detection-P filtering` | `856,001 CpGs; 118 samples` |
| `CpGs shared by pilot and GSE204943; r adjacent vs healthy` | `614,156; r = 0.97` |
| `group boundary (KDE minimum)` | `0.35` |
| `GSE204943 hypomethylated PMD blocks: number; total Mb` | `804; 100.3` |
| `RNA arm-dosage noise floor (99th percentile of absolute dosage in adjacent tissues)` | `0.25` |
| `genes excluded as composition-tracking` | `3,681` |
| `CpG-island promoters measured in P4 and P5; shared gain >= 0.2; of which KZFP` | `9,528; 50; 15` |
| `KZFP promoter CpGs in signature / unique; ZNF382 island CpGs also ZNF529 (Methods)` | `17 of 48; 13` |
| `sex-adjusted signature: CpGs retained; max pilot shift change` | `3,678 of 4,000; 0.007` |
| `tumour-only PMD blocks; Wilcoxon P (both block sets)` | `618 of 692; P = 3.1 x 10^-6 / 3.1 x 10^-6` |

## Keeping this document current

Update this document in the same commit as any change to an analysis, parameter, input or
result. The steps:

1. Run `make all`.
2. Compare `manuscript_numbers.csv` with the previous version (`tools/compare_tables.R`).
3. Update this document: procedure, parameters and the values table.
4. Update the manuscript text and `docs/CHANGELOG.md`.
5. Check that `make methods-check` passes. It runs at the end of `make all`.

The commit-msg hook (`make hooks`) and the GitHub workflow reject commits that change
`scripts/`, `R/lib/`, `resources/` or `env/renv.lock` without updating this file. A commit that
does not change the methods (a comment, plot styling, a refactor) must say so in its message:
`Methods-unchanged: <reason>`.
