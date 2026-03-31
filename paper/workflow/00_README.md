# Paper Workflow

Active paper-production entry points live here.

Recommended order:

1. `08_run_paper_analysis.sh`
   - runs the R paper-analysis pipeline and rebuilds final paper-facing figures/tables
2. `06_build_supp_tables_xlsx.py`
   - rebuilds the supplementary tables workbook from the final TSVs
3. `07_update_manuscript_2_1.sh`
   - updates the Word manuscript with the current figures, tables, text sync, and Arial font pass

Single-step full rebuild:

- `09_run_full_paper_build.sh`

Supporting manuscript-update steps:

1. `01_render_main_figures_for_doc.py`
2. `02_sync_manuscript_from_results.py`
3. `03_replace_main_figures_in_doc.py`
4. `04_patch_doc_consistency.py`
5. `05_polish_manuscript_2_1.py`

Primary outputs:

- manuscript: `paper/SAB_HET_2.1_JU.docx`
- supplementary workbook: `results/paper/final/tables/supplementary_tables.xlsx`
