#!/usr/bin/env python3
import shutil
from datetime import datetime
from pathlib import Path
from docx import Document

DOC = Path('paper/SAB_HET_2.0_edited.docx')

stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
backup = DOC.with_name(f'{DOC.stem}.backup_{stamp}.docx')
shutil.copy2(DOC, backup)

d = Document(str(DOC))

for p in d.paragraphs:
    t = p.text or ''
    if 'impact classification accuracy (50%-100%)' in t:
        p.text = t.replace('impact classification accuracy (50%-100%)', 'impact classification accuracy (70%-100%)')
        t = p.text
    if t.startswith('Table 3: Statistical power to detect a subphenotype specific effect given a sample size n in each subphenotype. Power estimated from 2,000 replicates.'):
        p.text = 'Table 3: Statistical power to detect a subphenotype specific effect given a sample size n in each subphenotype. Power estimated from 2,000 replicates where simulated.'

d.save(str(DOC))
print('Updated', DOC)
print('Backup', backup)
