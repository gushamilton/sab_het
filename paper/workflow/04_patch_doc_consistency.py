#!/usr/bin/env python3
import shutil
from datetime import datetime
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[2]


def resolve_doc_path() -> Path:
    paper_dir = ROOT / "paper"
    env_doc = __import__("os").environ.get("MANUSCRIPT_DOC")
    if env_doc:
        return Path(env_doc)
    preferred = paper_dir / "SAB_HET_2.1_JU.docx"
    if preferred.exists():
        return preferred
    raise FileNotFoundError(f"No manuscript doc found in {paper_dir}")


DOC = resolve_doc_path()

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
backup = DOC.with_name(f"{DOC.stem}.backup_{stamp}.docx")
shutil.copy2(DOC, backup)

d = Document(str(DOC))

for p in d.paragraphs:
    t = p.text or ""
    if "impact classification accuracy (50%-100%)" in t:
        p.text = t.replace("impact classification accuracy (50%-100%)", "impact of classification accuracy (70%-100%)")
    if t.startswith("Table 3: Statistical power to detect"):
        p.text = (
            "Table 2: Statistical power to detect a subphenotype-specific effect at each sample size n. "
            "Power was estimated from 2,000 simulated replicates. Subphenotype A has no treatment effect, "
            "so power is not the appropriate metric."
        )
    if "Table 3 and Figure 1" in t:
        p.text = t.replace("Table 3 and Figure 1", "Table 2 and Figure 1")

d.save(str(DOC))
print("Updated", DOC)
print("Backup", backup)
