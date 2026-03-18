#!/usr/bin/env python3
import csv
from pathlib import Path
from openpyxl import Workbook

ROOT = Path(__file__).resolve().parents[1]
TAB_DIR = ROOT / 'results' / 'paper' / 'final' / 'tables'
OUT_XLSX = TAB_DIR / 'supplementary_tables.xlsx'
ORDER = TAB_DIR / 'table_order.tsv'

# Only supplementary tables
TABLES = [
    ('Table S1', TAB_DIR / 'TableS1_main_binary_operating_characteristics.tsv'),
    ('Table S2', TAB_DIR / 'TableS2_cohort_required_n.tsv'),
    ('Table S3', TAB_DIR / 'TableS3_enrichment_summary.tsv'),
    ('Table S4', TAB_DIR / 'TableS4_ordinal_binary_metrics_bias_key_n.tsv'),
    ('Table S5', TAB_DIR / 'TableS5_ordinal_type1.tsv'),
    ('Table S6', TAB_DIR / 'TableS6_ordinal_power_bias_tradeoff_key_n.tsv'),
    ('Table S7', TAB_DIR / 'TableS7_ordinal_PO_nonPO_calibration.tsv'),
]


def read_tsv(path: Path):
    with path.open('r', encoding='utf-8', newline='') as f:
        r = csv.reader(f, delimiter='\t')
        rows = list(r)
    return rows


def main():
    wb = Workbook()
    # README sheet with labels/purpose
    ws = wb.active
    ws.title = 'README'
    ws.append(['table_id', 'path', 'purpose'])
    if ORDER.exists():
        with ORDER.open('r', encoding='utf-8', newline='') as f:
            r = csv.reader(f, delimiter='\t')
            next(r, None)
            for row in r:
                if row and row[0].startswith('Table S'):
                    ws.append(row)
    else:
        ws.append(['(missing)', str(ORDER), 'table_order.tsv not found'])

    for name, path in TABLES:
        if not path.exists():
            continue
        # Excel sheet names max 31 chars
        sheet_name = name if len(name) <= 31 else name[:31]
        ws = wb.create_sheet(title=sheet_name)
        rows = read_tsv(path)
        for row in rows:
            ws.append(row)

    wb.save(OUT_XLSX)
    print(f'Wrote {OUT_XLSX}')


if __name__ == '__main__':
    main()
