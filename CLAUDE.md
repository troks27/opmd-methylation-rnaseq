# Working in this repository

Analysis code for the OPMD methylation + RNA-seq manuscript. Start with `README.md`. The
methods are described in `docs/methods.md`, the scripts' inputs and outputs in
`docs/pipeline.md`, and the history of changes in `docs/CHANGELOG.md`.

## Rule: the methods document follows the analysis

`docs/methods.md` must always describe what the code does and quote the current results.
Whenever an analysis, parameter, input or result changes:

1. Edit the script(s). Never edit a script while `make` is running it: Rscript reads scripts
   incrementally, so the running job would parse a corrupted file.
2. `make all` reruns what changed. Run `make -n all` first to see the list.
3. Compare `outputs/extended/tables/manuscript_numbers.csv` with the previous version
   (`Rscript tools/compare_tables.R <old> <new>`) to find the sentences that change.
4. Update `docs/methods.md`. This covers the procedure, the parameters, and the table
   "Values that come from the data". Then update the manuscript text and `docs/CHANGELOG.md`.
5. `make methods-check` must pass. It runs at the end of `make all`.
6. Commit methods and code together. The commit-msg hook (`make hooks`) and the GitHub
   workflow reject a commit that changes `scripts/`, `R/lib/`, `resources/` or
   `env/renv.lock` without `docs/methods.md`. They accept it only if the commit message has a
   `Methods-unchanged: <reason>` line, which is for comments, plot styling or refactors.

## Other rules

- The repository is public. Never commit data, outputs, local paths, specimen identifiers
  (the pilot's 99xxxxx numbers) or library file names. `config/paths.local.yml` is git-ignored.
- Every number in the manuscript must come from a script. `06_audit/03_manuscript_numbers.R`
  lists them, and `06_audit/01_results_audit.R` checks the key claims against their reported
  values.
- Fix random seeds for anything stochastic.
- Figures: Nature format (183 mm wide, ≤ 170 mm tall, Helvetica 5-7 pt, Quartz PDF on macOS).
  Check rendered PNGs for overlapping labels after every change.
