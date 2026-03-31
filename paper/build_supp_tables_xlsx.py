#!/usr/bin/env python3
import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "workflow" / "06_build_supp_tables_xlsx.py"), run_name="__main__")
