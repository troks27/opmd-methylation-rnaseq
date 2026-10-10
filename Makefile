# Pipeline for "DNA methylation and transcriptome analysis of oral potentially malignant
# lesions". Run from the repository root.
#
#   make help        list targets
#   make all         every analysis, audit and figure (stages 02-07)
#   make figures     Figs 1-6 and Supplementary Figs 1-3
#
# Each target is a file written by one script. Its prerequisites are the script itself,
# the shared R/lib code and the outputs of upstream scripts. After you edit a script,
# `make all` reruns that script and everything downstream of it, and nothing else.
# `make -n all` shows what would run without running it.
#
# Data that Make does not track (fixed inputs): GEO downloads (make download), the
# pilot IDATs and sample sheet, the GSE227919 files and the pilot RNA-seq FASTQs.
# Locations are set in config/paths.yml (+ config/paths.local.yml).

SHELL := /bin/bash
R     := Rscript --no-save --no-restore
# Output locations are read from the config (and OPMD_* environment overrides) on every run.
$(shell Rscript R/lib/make_paths.R .make_paths.mk)
include .make_paths.mk

# Shared code. The path files (config/paths*.yml) are deliberately not prerequisites: moving data
# should not trigger a rebuild. After pointing the config at different input data, force the
# affected stages, e.g. `make -B methylation`.
LIB := R/lib/config.R R/lib/utils.R R/lib/sample_sheet.R
S   := scripts

# ------------------------------------------------------------------ sentinel outputs
G01          := $(OUT)/gse204943
GSE_BETA     := $(G01)/crucial_objects/beta.rds
GSE_CELLS    := $(G01)/results/cell_type_proportions.csv
GSE_PMD      := $(G01)/results/hypo_PMD_blocks.csv
GSE_CNV      := $(G01)/results/CNV_segments_all.csv
GSE_ZNF      := $(G01)/results/ZNF471_CpG_details.csv
# Downstream rules depend on these only when stage 01 builds the objects in $(OUT).
GSE_OBJ      := $(if $(GSE_BUILT),$(GSE_BETA) $(GSE_CELLS) $(GSE_PMD) $(GSE_ZNF))
GSE_CNV_DEP  := $(if $(GSE_BUILT),$(GSE_CNV))
PROCESSED    := $(OUT)/processed.rds
CNV          := $(EXT)/cnv_objects.rds
DMP          := $(EXT)/tables/dmp_lesion_vs_normal.csv
EXTMETH      := $(EXT)/external_methylation.rds
STORY        := $(EXT)/story_analyses.rds
MUMBAI       := $(EXT)/mumbai_crosscheck.rds
KZFP3        := $(EXT)/gse234379_kzfp.rds
RNA227919    := $(EXT)/external_rnaseq.rds
G26549       := $(EXT)/gse26549_progression.rds
MATRICES     := $(MAT)/gene_tpm.tsv.gz
RNASEQ_STRAND := $(RNASEQ)/strandedness_check.tsv
RNAVALID     := $(EXT)/tables/rnaseq_proxy_programmes.csv
RNASTORY     := $(EXT)/rnaseq_story.rds
AUDIT        := $(EXT)/audit.rds
VERIFY       := $(EXT)/tables/verification_kzfp_cpg_accounting.csv
NUMBERS      := $(EXT)/tables/manuscript_numbers.csv
INTEGRITY    := $(EXT)/tables/integrity_gse204943_summary.csv
FIG1_4       := $(FIGDIR)/Fig1.pdf
FIG5_6       := $(FIGDIR)/Fig5.pdf

.PHONY: help all methylation cohorts outcome rnaseq audit figures gse204943 download kallisto exploratory clean-outputs

help:
	@echo "Stages (see README.md and docs/pipeline.md):"
	@echo "  make download     00  download public GEO data into data/geo (~2.5 GB)"
	@echo "  make gse204943    01  GSE204943: minfi betas, cell fractions, PMD blocks, conumee2 copy number"
	@echo "  make methylation  02  pilot EPIC v2: sesame preprocessing, sex, genotype identity, copy number"
	@echo "  make cohorts      03  tumour signature/score, OPL groups, Mumbai cross-checks, GSE234379"
	@echo "  make outcome      04  GSE227919 programmes, GSE26549 transformation analysis"
	@echo "  make kallisto     05a pilot RNA-seq quantification (needs FASTQ access; see script header)"
	@echo "  make rnaseq       05  gene matrices, RNA vs methylation validation, RNA analyses"
	@echo "  make audit        06  results audit and cohort verification (assert reported numbers)"
	@echo "  make figures      07  Figs 1-6, Supplementary Figs 1-3 (PDF + PNG)"
	@echo "  make all          stages 02-07"
	@echo "  make exploratory  pilot paired DMP/DMR candidates (not in the manuscript)"
	@echo "Outputs: $(OUT)   (set in config/paths.yml)"

all: figures audit

# ------------------------------------------------------------------ 00 download
download:
	bash $(S)/00_download/download_geo.sh

# ------------------------------------------------------------------ 01 GSE204943
# Run explicitly (`make gse204943`); with the default config, `make all` also builds it.
gse204943: $(GSE_BETA) $(GSE_CELLS) $(GSE_PMD) $(GSE_CNV) $(GSE_ZNF)
$(GSE_BETA): $(S)/01_gse204943/01_minfi_preprocess.R $(LIB)
	$(R) $<
$(GSE_CELLS): $(S)/01_gse204943/02_cell_composition.R $(GSE_BETA) $(LIB)
	$(R) $<
$(GSE_PMD): $(S)/01_gse204943/03_pmd_blocks.R $(GSE_BETA) $(LIB)
	$(R) $<
$(GSE_CNV): $(S)/01_gse204943/04_copy_number_conumee2.R $(GSE_BETA) $(LIB)
	$(R) $<
$(GSE_ZNF): $(S)/01_gse204943/05_locus_cpg_sets.R $(GSE_BETA) $(LIB)
	$(R) $<

# ------------------------------------------------------------------ 02 pilot methylation
methylation: $(PROCESSED) $(CNV)
$(PROCESSED): $(S)/02_pilot_methylation/01_preprocess_sesame.R $(LIB)
	$(R) $<
$(CNV): $(S)/02_pilot_methylation/02_identity_sex_copy_number.R $(LIB)
	$(R) $<
exploratory: $(DMP)
$(DMP): $(S)/02_pilot_methylation/03_exploratory_dmp_dmr.R $(PROCESSED) $(LIB)
	$(R) $<

# ------------------------------------------------------------------ 03 methylation cohorts
cohorts: $(EXTMETH) $(STORY) $(MUMBAI) $(KZFP3)
$(EXTMETH): $(S)/03_methylation_cohorts/01_composition_signature_scores.R $(PROCESSED) $(GSE_OBJ) $(LIB)
	$(R) $<
$(STORY): $(S)/03_methylation_cohorts/02_opl_groups.R $(EXTMETH) $(GSE_CNV_DEP) $(LIB)
	$(R) $<
$(MUMBAI): $(S)/03_methylation_cohorts/03_mumbai_crosscheck.R $(EXTMETH) $(STORY) $(LIB)
	$(R) $<
$(KZFP3): $(S)/03_methylation_cohorts/04_gse234379_kzfp.R $(EXTMETH) $(STORY) $(PROCESSED) $(LIB)
	$(R) $<

# ------------------------------------------------------------------ 04 expression cohorts and outcome
outcome: $(RNA227919) $(G26549)
$(RNA227919): $(S)/04_expression_outcome/01_gse227919_programmes.R $(EXTMETH) $(PROCESSED) $(LIB)
	$(R) $<
$(G26549): $(S)/04_expression_outcome/02_gse26549_transformation.R $(RNA227919) $(MUMBAI) $(LIB)
	$(R) $<

# ------------------------------------------------------------------ 05 pilot RNA-seq
kallisto:
	bash $(S)/05_pilot_rnaseq/00_strandedness_check.sh
	bash $(S)/05_pilot_rnaseq/01_kallisto_sequential.sh
rnaseq: $(MATRICES) $(RNAVALID) $(RNASTORY)
$(MATRICES): $(S)/05_pilot_rnaseq/02_gene_matrices.R $(LIB)
	$(R) $<
$(RNAVALID): $(S)/05_pilot_rnaseq/03_rna_methylation_validation.R $(MATRICES) $(CNV) $(EXTMETH) $(G26549) $(KZFP3) $(LIB)
	$(R) $<
$(RNASTORY): $(S)/05_pilot_rnaseq/04_rna_analyses.R $(MATRICES) $(PROCESSED) $(EXTMETH) $(RNA227919) $(G26549) $(LIB)
	$(R) $<

# ------------------------------------------------------------------ 06 audit
audit: $(AUDIT) $(VERIFY) $(NUMBERS) $(INTEGRITY)
$(AUDIT): $(S)/06_audit/01_results_audit.R $(EXTMETH) $(STORY) $(MUMBAI) $(KZFP3) $(G26549) $(CNV) $(LIB)
	$(R) $<
$(VERIFY): $(S)/06_audit/02_cohort_verification.R $(EXTMETH) $(STORY) $(KZFP3) $(G26549) $(PROCESSED) $(LIB)
	$(R) $<

$(NUMBERS): $(S)/06_audit/03_manuscript_numbers.R $(AUDIT) $(VERIFY) $(RNASTORY) $(RNAVALID) $(LIB)
	$(R) $<

$(INTEGRITY): $(S)/06_audit/04_gse204943_integrity.R $(EXTMETH) $(PROCESSED) $(GSE_OBJ) $(GSE_CNV_DEP) $(LIB)
	$(R) $<

# ------------------------------------------------------------------ 07 figures
figures: $(FIG1_4) $(FIG5_6)
$(FIG1_4): $(S)/07_figures/01_figures_1_to_4.R $(EXTMETH) $(CNV) $(STORY) $(MUMBAI) $(G26549) $(KZFP3) $(RNA227919) $(AUDIT) $(PROCESSED) $(LIB)
	$(R) $<
$(FIG5_6): $(S)/07_figures/02_figures_5_6.R $(RNASTORY) $(RNAVALID) $(G26549) $(AUDIT) $(RNASEQ_STRAND) resources/cnv_events_curated.csv $(LIB)
	$(R) $<

# ------------------------------------------------------------------ housekeeping
clean-outputs:
	@echo "This deletes $(OUT)/extended and $(OUT)/figures_nature (not stage 01 or RNA-seq outputs)."
	@read -p "Continue? [y/N] " a; [ "$$a" = y ] && rm -rf $(EXT) $(FIGDIR) || echo "Cancelled."
