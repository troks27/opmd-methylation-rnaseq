# Curated inputs

Small tables that record judgements or results from outside this code. Each one has a
`note` or `source` column; scripts read them rather than hard-coding the values.

| File | Used by | Content |
|---|---|---|
| `cnv_events_curated.csv` | `scripts/07_figures/02_figures_5_6.R` (Fig. 5f) | Copy-number events called on the pilot methylation arrays, whether exome sequencing found them (from the multimodal manuscript; the exome data are not analysed here), and 9p21 calls from RNA. An empty `rna` cell is computed from arm-level expression dosage: detected when the mean dosage of `arms` passes the RNA noise floor in the stated `direction`. |
