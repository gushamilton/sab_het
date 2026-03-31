#!/usr/bin/env python3
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "results" / "paper" / "final" / "figures"
OUT_DIR = ROOT / "paper" / "fig_tmp"

FIG_MAP = {
    "Figure1_misclassification_accuracy.pdf": "Figure1.png",
    "Figure2_cohort_conditionality.pdf": "Figure2.png",
    "Figure3_enrichment_feasibility.pdf": "Figure3.png",
    "Figure4_ordinal_category_shift.pdf": "Figure4.png",
    "Figure5_ordinal_binary_summary.pdf": "Figure5.png",
}


def render_pdf_to_png(src: Path, dst: Path) -> None:
    cmd = [
        "gs",
        "-dSAFER",
        "-dBATCH",
        "-dNOPAUSE",
        "-sDEVICE=pngalpha",
        "-r300",
        f"-sOutputFile={dst}",
        str(src),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in FIG_MAP.items():
        src = FIG_DIR / src_name
        dst = OUT_DIR / dst_name
        if not src.exists():
            raise FileNotFoundError(f"Missing source figure: {src}")
        render_pdf_to_png(src, dst)
        print(f"Wrote {dst}")


if __name__ == "__main__":
    main()
