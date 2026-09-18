#!/usr/bin/env python3
"""Build an epiSCANEER input alignment for a UniProt query protein.

Pipeline: fetch the query -> collect its Pfam family from UniProtKB ->
score every member against the query by local alignment -> filter on
coverage/identity -> subsample -> drop near-duplicates -> align with the
bundled MUSCLE -> write CLUSTAL with the query as the first sequence.

Example:
    python scripts/build_msa.py --accession P38567 --pfam PF01630 \
        --prefix PH20 --out input/PH20.aln
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile
import urllib.request

from Bio import Align, AlignIO, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
UNIPROT = "https://rest.uniprot.org"


def fetch_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode()), r.headers.get("Link", "")


def fetch_query(accession):
    data, _ = fetch_json("%s/uniprotkb/%s.json?fields=accession,sequence,protein_name"
                         % (UNIPROT, accession))
    return data["sequence"]["value"]


def fetch_family(pfam, cache=None):
    if cache and os.path.isfile(cache):
        return json.load(open(cache))
    url = ("%s/uniprotkb/search?query=xref:pfam-%s&format=json&size=500"
           "&fields=accession,protein_name,organism_name,length,sequence" % (UNIPROT, pfam))
    out = []
    while url:
        data, link = fetch_json(url)
        out.extend(data["results"])
        url = link.split("<")[1].split(">")[0] if 'rel="next"' in link else None
        print("  fetched %d sequences" % len(out), file=sys.stderr)
    if cache:
        json.dump(out, open(cache, "w"))
    return out


def make_aligner():
    return Align.PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
    )


def identity(aligner, a, b):
    """Identity over aligned residues, and coverage of `a`."""
    counts = aligner.align(a, b)[0].counts()
    aligned = counts.identities + counts.mismatches
    if aligned == 0:
        return 0.0, 0.0
    return counts.identities / aligned, aligned / len(a)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--accession", required=True, help="UniProt accession of the query")
    p.add_argument("--pfam", required=True, help="Pfam family to draw homologs from")
    p.add_argument("--prefix", required=True, help="name of the query in the alignment")
    p.add_argument("--out", required=True, help="output CLUSTAL file")
    p.add_argument("--target", type=int, default=400, help="number of homologs to keep")
    p.add_argument("--min-coverage", type=float, default=0.70)
    p.add_argument("--min-identity", type=float, default=0.30)
    p.add_argument("--max-identity", type=float, default=0.95)
    p.add_argument("--max-redundancy", type=float, default=0.90,
                   help="drop homologs more identical than this to a kept homolog")
    p.add_argument("--min-length", type=int, default=300)
    p.add_argument("--max-length", type=int, default=700)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--cache", help="json cache of the family download")
    p.add_argument("--report", help="write the selected homologs as TSV")
    args = p.parse_args()

    print("1. query %s" % args.accession, file=sys.stderr)
    query_seq = fetch_query(args.accession)
    print("   %d residues" % len(query_seq), file=sys.stderr)

    print("2. family %s" % args.pfam, file=sys.stderr)
    family = fetch_family(args.pfam, args.cache)

    aligner = make_aligner()
    print("3. scoring %d family members against the query" % len(family), file=sys.stderr)
    pool = []
    for entry in family:
        seq = entry["sequence"]["value"]
        if entry["primaryAccession"] == args.accession:
            continue
        if not args.min_length <= len(seq) <= args.max_length:
            continue
        if not set(seq) <= STANDARD_AA:
            continue
        ident, cov = identity(aligner, query_seq, seq)
        if cov < args.min_coverage:
            continue
        if not args.min_identity <= ident <= args.max_identity:
            continue
        pool.append({
            "acc": entry["primaryAccession"],
            "organism": entry["organism"]["scientificName"],
            "length": len(seq),
            "identity": ident,
            "coverage": cov,
            "seq": seq,
        })
    print("   %d homologs pass coverage/identity filters" % len(pool), file=sys.stderr)

    random.Random(args.seed).shuffle(pool)
    print("4. removing near-duplicates (>%.0f%% identity)" % (args.max_redundancy * 100),
          file=sys.stderr)
    kept = []
    for cand in pool:
        if len(kept) >= args.target:
            break
        if any(identity(aligner, k["seq"], cand["seq"])[0] > args.max_redundancy for k in kept):
            continue
        kept.append(cand)
        if len(kept) % 50 == 0:
            print("   kept %d" % len(kept), file=sys.stderr)
    kept.sort(key=lambda r: -r["identity"])
    print("   %d homologs selected" % len(kept), file=sys.stderr)

    if args.report:
        with open(args.report, "w") as f:
            print("accession\torganism\tlength\tidentity_to_query\tquery_coverage", file=f)
            for r in kept:
                print("%s\t%s\t%d\t%.4f\t%.4f"
                      % (r["acc"], r["organism"], r["length"], r["identity"], r["coverage"]), file=f)

    records = [SeqRecord(Seq(query_seq), id=args.prefix, description="")]
    records += [SeqRecord(Seq(r["seq"]), id=r["acc"], description="") for r in kept]

    print("5. aligning %d sequences with MUSCLE" % len(records), file=sys.stderr)
    muscle = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "epiSCANEER", "msa_exec", "linux", "muscle")
    with tempfile.TemporaryDirectory() as tmp:
        fasta_in = os.path.join(tmp, "in.fasta")
        fasta_out = os.path.join(tmp, "out.fasta")
        SeqIO.write(records, fasta_in, "fasta")
        proc = subprocess.run([muscle, "-in", fasta_in, "-out", fasta_out],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            sys.exit("MUSCLE failed:\n%s" % proc.stderr[-2000:])
        aligned = {r.id: r for r in SeqIO.parse(fasta_out, "fasta")}

    # MUSCLE reorders sequences; epiSCANEER requires the query first.
    ordered = [aligned[args.prefix]] + [aligned[r.id] for r in records[1:] if r.id in aligned]
    for rec in ordered:
        rec.description = ""
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    AlignIO.write(Align.MultipleSeqAlignment(ordered), args.out, "clustal")
    print("6. wrote %s (%d sequences, %d columns)"
          % (args.out, len(ordered), len(ordered[0].seq)), file=sys.stderr)


if __name__ == "__main__":
    main()
