#!/usr/bin/env python3
"""Annotate epiSCANEER PH20 predictions with UniProt and structural context.

Adds, for every predicted residue pair:
  - mature-protein numbering (P38567 position - 35), the numbering used by
    the rHuPH20 / Halozyme literature
  - which part of the precursor the residue belongs to
  - distance to the catalytic proton donor Glu148, measured on PDB 9JUB
    chain A (PH20 residues 37-442, deposited in P38567 numbering)
  - disulfide / N-glycosylation flags

and writes a filtered shortlist that keeps only pairs usable in the
secreted enzyme.
"""
import argparse
import csv
import json
import os
import urllib.request
import warnings

warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE_TO_ONE

ACCESSION = "P38567"
SIGNAL_PEPTIDE = (1, 35)
RHUPH20 = (36, 482)          # the secreted drug substance (Hylenex)
GPI_TAIL = (483, 509)
ACTIVE_SITE = 148            # Glu148, proton donor (UniProt)
DISULFIDES = [(60, 351), (224, 238), (376, 387), (381, 435), (437, 464)]
NGLYC = [82, 166, 235, 254, 368, 393]
DISULFIDE_CYS = {p for pair in DISULFIDES for p in pair}


def region(pos):
    if SIGNAL_PEPTIDE[0] <= pos <= SIGNAL_PEPTIDE[1]:
        return "signal_peptide"
    if RHUPH20[0] <= pos <= RHUPH20[1]:
        return "rHuPH20"
    if GPI_TAIL[0] <= pos <= GPI_TAIL[1]:
        return "GPI_tail"
    return "?"


def load_structure(path, url):
    if not os.path.isfile(path):
        urllib.request.urlretrieve(url, path)
    chain = next(iter(MMCIFParser(QUIET=True).get_structure("s", path)))["A"]
    return {r.id[1]: r for r in chain if r.id[0] == " "}


def min_heavy_distance(res_a, res_b):
    return min(a - b for a in res_a if a.element != "H" for b in res_b if b.element != "H")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="output/PH20/PH20.txt")
    p.add_argument("--structure", default="9jub.cif")
    p.add_argument("--structure-url", default="https://files.rcsb.org/download/9JUB.cif")
    p.add_argument("--out", default="output/PH20/PH20.annotated.tsv")
    p.add_argument("--shortlist", default="output/PH20/PH20.shortlist.tsv")
    args = p.parse_args()

    residues = load_structure(args.structure, args.structure_url)
    active = residues.get(ACTIVE_SITE)
    if active is None:
        raise SystemExit("Glu%d missing from the structure" % ACTIVE_SITE)

    def distance_to_active_site(pos):
        res = residues.get(pos)
        if res is None:
            return ""
        return round(min_heavy_distance(res, active), 2)

    with open(args.results) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    out_fields = ["pos1", "pos2", "mature_pos1", "mature_pos2", "region1", "region2",
                  "WT_AA1", "WT_AA2", "MUT_AA1", "MUT_AA2", "mutation",
                  "WT_freq", "MUT_SCI", "SCI_norm",
                  "dist1_to_E148", "dist2_to_E148", "flags"]
    annotated = []
    for r in rows:
        p1, p2 = int(r["pos1"]), int(r["pos2"])
        flags = []
        for pos in (p1, p2):
            if pos in DISULFIDE_CYS:
                flags.append("disulfide_Cys%d" % pos)
            if pos in NGLYC:
                flags.append("N-glyc_site%d" % pos)
            if pos == ACTIVE_SITE:
                flags.append("catalytic_E148")
        annotated.append({
            "pos1": p1, "pos2": p2,
            "mature_pos1": p1 - 35, "mature_pos2": p2 - 35,
            "region1": region(p1), "region2": region(p2),
            "WT_AA1": r["WT_AA1"], "WT_AA2": r["WT_AA2"],
            "MUT_AA1": r["MUT_AA1"], "MUT_AA2": r["MUT_AA2"],
            "mutation": "%s%d%s + %s%d%s" % (r["WT_AA1"], p1, r["MUT_AA1"],
                                             r["WT_AA2"], p2, r["MUT_AA2"]),
            "WT_freq": r["WT_freq"], "MUT_SCI": r["MUT_SCI"], "SCI_norm": r["SCI_norm"],
            "dist1_to_E148": distance_to_active_site(p1),
            "dist2_to_E148": distance_to_active_site(p2),
            "flags": ";".join(flags),
        })

    annotated.sort(key=lambda r: -float(r["SCI_norm"]))
    for path, rows_out in ((args.out, annotated),
                           (args.shortlist, [r for r in annotated
                                             if r["region1"] == "rHuPH20"
                                             and r["region2"] == "rHuPH20"
                                             and not r["flags"]])):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows_out)
        print("wrote %s (%d rows)" % (path, len(rows_out)))


if __name__ == "__main__":
    main()
