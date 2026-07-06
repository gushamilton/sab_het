# SAB HTE Simulation Study

This repository contains the code and final materials for the SAB heterogeneous treatment effect simulation study.

## What is included

- `paper/`: manuscript files, supplementary methods, and paper figure-generation helpers.
- `code/`: analysis scripts used to generate the paper results.
- `results/paper/final/`: final figure manifest and tables used in the manuscript.
- `scripts/`: local helper scripts.
- `slurm/`: batch submission scripts.

## Key outputs

- Supplementary methods: `paper/supplementary_methods.md`
- Final figure manifest: `results/paper/final/figures/figure_order.tsv`
- Final tables: `results/paper/final/tables/`
- Publication asset manifest: `results/paper/final/README.md`

## Reproducibility

The final paper assets are built from the scripts in `code/` and the manuscript helpers in `paper/`.
Figure and table ordering is defined in `results/paper/final/figures/figure_order.tsv` and `results/paper/final/tables/table_order.tsv`.

The committed final tables are the publication-facing assets used to check the manuscript values. Figure provenance is documented in the final manifest; one small supplementary figure binary still needs a normal credentialed Git push. Large raw simulation outputs are reproducible from the analysis scripts and are not all tracked in Git.

## Notes

The local archive under `old/` is kept for working history and is not part of the public release.
