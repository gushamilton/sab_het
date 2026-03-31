#!/usr/bin/env bash
set -euo pipefail

cd /user/work/fh6520/bc4/sab_het

bash paper/workflow/08_run_paper_analysis.sh
python3 paper/workflow/06_build_supp_tables_xlsx.py
bash paper/workflow/07_update_manuscript_2_1.sh

echo "Full paper build complete."
