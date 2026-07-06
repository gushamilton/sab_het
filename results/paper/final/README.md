# Final Paper Assets

This directory contains the publication-facing tables and figure manifest for the SAB heterogeneous treatment effect simulation paper.

## Scope

- `figures/`: final figure PDFs used for the manuscript where committed as repository binaries.
- `tables/`: final manuscript and supplementary tables.
- `figures/figure_order.tsv`: maps manuscript figure labels to generated source files.
- `tables/table_order.tsv`: maps manuscript table labels to final table files.
- `tables/raw_data_provenance.tsv`: records the intermediate analysis outputs used to generate the final tables.

## Manuscript Consistency Checks

The committed final tables match the submission manuscript values for the main numerical claims:

- Table 1 parameters: subgroup frequencies, ARREST 84-day mortality, ARREST ORs, and shrunk ORs.
- Table 2 power grid: subgroup B reaches power 1.000 at n=3,000; subgroups C, D, and E are 0.486, 0.527, and 0.708 at n=20,000.
- Enrichment results: subgroup B number-needed-to-screen is approximately 1,100 with a near-perfect classifier and 3,300 with a balanced classifier; corresponding NNR is approximately 140 to 930.
- Ordinal outcome comparison: the n=3,000, 100% accuracy proportional-odds and death-only non-proportional power values are in `tables/TableS8_ordinal_PO_vs_nonPO_power_bias.tsv`.

Large raw simulation outputs are reproducible from the scripts in `code/` and are not all tracked in Git. The final publication tables are tracked here so readers can audit the values reported in the paper. The small Figure S1 PDF is present in the prepared local publication commit and listed in `figures/figure_order.tsv`; uploading that binary to GitHub requires a normal credentialed Git push because the current connector path is text-oriented.
