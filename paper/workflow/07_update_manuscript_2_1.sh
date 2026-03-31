#!/usr/bin/env bash
set -euo pipefail

cd /user/work/fh6520/bc4/sab_het

export MANUSCRIPT_DOC="${MANUSCRIPT_DOC:-paper/SAB_HET_2.1_JU.docx}"

python3 paper/workflow/01_render_main_figures_for_doc.py
python3 paper/workflow/02_sync_manuscript_from_results.py
python3 paper/workflow/03_replace_main_figures_in_doc.py
python3 paper/workflow/04_patch_doc_consistency.py
python3 paper/workflow/05_polish_manuscript_2_1.py

echo "Updated manuscript: ${MANUSCRIPT_DOC}"
