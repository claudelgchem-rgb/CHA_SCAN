#!/usr/bin/env python3
"""Render an epiSCANEER shortlist as a Top-N workbook.

Reads the annotated shortlist written by annotate_ph20.py, sorts it by
SCI_norm and writes one sheet per requested size plus a README sheet whose
rows come from a key/value TSV, so the run's facts live in data rather than
in this script.

Values are written as literals, not formulas: this environment has no
working spreadsheet engine to cache formula results with, and a formula
without a cached value reads back empty in every tool that is not Excel.

    python scripts/make_report_xlsx.py \
        --shortlist results/PH20_paper500/PH20.shortlist.tsv \
        --readme results/PH20_paper500/README.tsv \
        --out results/PH20_paper500/PH20_top100_top50.xlsx \
        --sizes 100 50
"""
import argparse
import csv

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

FONT = "Arial"
HDR_FILL = PatternFill("solid", fgColor="1F3864")
HDR_FONT = Font(name=FONT, size=10, bold=True, color="FFFFFF")
BODY = Font(name=FONT, size=10)
NEAR_FILL = PatternFill("solid", fgColor="FFF2CC")
BORDER = Border(bottom=Side(style="thin", color="BFBFBF"))
NEAR_ANGSTROM = 10.0

COLS = [
    ("Rank", 6, None), ("Mutation (P38567)", 22, None),
    ("Mutation (mature / rHuPH20)", 26, None),
    ("pos1", 7, "0"), ("pos2", 7, "0"),
    ("mature_pos1", 12, "0"), ("mature_pos2", 12, "0"),
    ("WT_AA1", 8, None), ("MUT_AA1", 9, None),
    ("WT_AA2", 8, None), ("MUT_AA2", 9, None),
    ("WT_freq", 10, "0.0000"), ("MUT_SCI", 11, "0.000000"),
    ("SCI_norm", 10, "0.0000"),
    ("dist1_to_E148 (A)", 16, "0.00"), ("dist2_to_E148 (A)", 16, "0.00"),
    ("min_dist_to_E148 (A)", 19, "0.00"),
    ("Near active site (<=10 A)", 22, None),
]
CENTERED = {1, 4, 5, 6, 7, 8, 9, 10, 11, 18}


def write_sheet(ws, data):
    for i, (name, width, _fmt) in enumerate(COLS, start=1):
        c = ws.cell(row=1, column=i, value=name)
        c.font, c.fill = HDR_FONT, HDR_FILL
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.column_dimensions[get_column_letter(i)].width = width
    ws.row_dimensions[1].height = 30

    for n, r in enumerate(data, start=1):
        row = n + 1
        p1, p2 = int(r["pos1"]), int(r["pos2"])
        d1 = float(r["dist1_to_E148"]) if r["dist1_to_E148"] else None
        d2 = float(r["dist2_to_E148"]) if r["dist2_to_E148"] else None
        dists = [d for d in (d1, d2) if d is not None]
        dmin = min(dists) if dists else None
        values = [
            n,
            "%s%d%s + %s%d%s" % (r["WT_AA1"], p1, r["MUT_AA1"], r["WT_AA2"], p2, r["MUT_AA2"]),
            "%s%d%s + %s%d%s" % (r["WT_AA1"], p1 - 35, r["MUT_AA1"], r["WT_AA2"], p2 - 35, r["MUT_AA2"]),
            p1, p2, p1 - 35, p2 - 35,
            r["WT_AA1"], r["MUT_AA1"], r["WT_AA2"], r["MUT_AA2"],
            float(r["WT_freq"]), float(r["MUT_SCI"]), float(r["SCI_norm"]),
            d1, d2, dmin,
            "Y" if (dmin is not None and dmin <= NEAR_ANGSTROM) else "",
        ]
        near = dmin is not None and dmin <= NEAR_ANGSTROM
        for i, v in enumerate(values, start=1):
            c = ws.cell(row=row, column=i, value=v)
            c.font, c.border = BODY, BORDER
            fmt = COLS[i - 1][2]
            if fmt:
                c.number_format = fmt
            if i in CENTERED:
                c.alignment = Alignment(horizontal="center")
            if near:
                c.fill = NEAR_FILL

    ws.freeze_panes = "C2"
    ws.auto_filter.ref = "A1:%s%d" % (get_column_letter(len(COLS)), len(data) + 1)


def write_readme(ws, title, rows):
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 104
    c = ws.cell(row=1, column=1, value=title)
    c.font = Font(name=FONT, size=14, bold=True)
    ws.merge_cells("A1:B1")
    for i, (k, v) in enumerate(rows, start=3):
        a = ws.cell(row=i, column=1, value=k)
        a.font = Font(name=FONT, size=10, bold=bool(k))
        a.alignment = Alignment(vertical="top")
        b = ws.cell(row=i, column=2, value=v)
        b.font = BODY
        b.alignment = Alignment(vertical="top", wrap_text=True)
        if len(v) > 110:
            ws.row_dimensions[i].height = 14 * (1 + len(v) // 110)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shortlist", required=True)
    p.add_argument("--readme", required=True, help="TSV of key<TAB>value rows")
    p.add_argument("--out", required=True)
    p.add_argument("--title", default="epiSCANEER report")
    p.add_argument("--sizes", nargs="+", type=int, default=[100, 50])
    args = p.parse_args()

    with open(args.shortlist) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows.sort(key=lambda r: -float(r["SCI_norm"]))

    readme = []
    with open(args.readme) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                readme.append(("", ""))
                continue
            k, _, v = line.partition("\t")
            readme.append((k, v))

    wb = Workbook()
    write_readme(wb.active, args.title, readme)
    wb.active.title = "README"
    for n in args.sizes:
        write_sheet(wb.create_sheet("Top%d" % n), rows[:n])
    wb.save(args.out)
    print("wrote %s | %s" % (args.out, ", ".join("Top%d: %d rows" % (n, min(n, len(rows)))
                                                 for n in args.sizes)))


if __name__ == "__main__":
    main()
