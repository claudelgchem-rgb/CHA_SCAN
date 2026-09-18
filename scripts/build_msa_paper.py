#!/usr/bin/env python3
"""Build an epiSCANEER input MSA exactly as described in the epiSCANEER paper.

Reference: Nucleic Acids Research 2026, 54, gkag879, Methods, "MSA generation":

    Homologous protein sequences corresponding to the query were retrieved
    from the UniRef90 database, which consists of non-redundant protein
    sequences clustered at 90% sequence identity, using PSI-BLAST. The
    PSI-BLAST search was performed with a single iteration and an e-value
    threshold of 0.001. To ensure the quality of the resulting MSAs, only
    sequences whose lengths ranged from 0.7 to 1.3 times that of the query
    sequence were retained, and sequences with 90% or more identity to the
    query were excluded. MSAs were constructed using MUSCLE with the
    ClustalW output format and default parameter settings.

No other filter is applied: no coverage cut, no mutual-identity cut, no cap
on the number of homologs, no subsampling. The only departure from a literal
reading is that the aligned records are re-ordered so the query comes first,
which epiSCANEER requires (it reads the query as the first record) and which
upstream's own msa.ProcMsa does via its MUSCLE "stable" patch. Re-ordering
rows does not change the alignment columns.

PSI-BLAST runs on the EBI REST service because UniRef90 (32 GB compressed)
does not fit in this environment; the EBI service caps reporting at 5000
hits, which is recorded in the report.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

UNIPROT = "https://rest.uniprot.org"
EBI = "https://www.ebi.ac.uk/Tools/services/rest/psiblast"
MAX_HITS = "5000"          # the largest value the EBI service accepts
EVALUE = "1.0e-3"          # paper: e-value threshold of 0.001
LEN_LO, LEN_HI = 0.7, 1.3  # paper: 0.7 to 1.3 times the query length
MAX_IDENTITY = 90.0        # paper: sequences with >=90% identity excluded


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def http(url, data=None, timeout=300):
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()


def fetch_query(accession):
    data = json.loads(http("%s/uniprotkb/%s.json?fields=accession,sequence" % (UNIPROT, accession)))
    return data["sequence"]["value"]


def submit(seq, email, prefix):
    payload = urllib.parse.urlencode({
        "email": email,
        "database": "uniref90",
        "expthr": EVALUE,
        "scores": MAX_HITS,
        "alignments": MAX_HITS,
        "sequence": ">%s\n%s\n" % (prefix, seq),
    }).encode()
    job = http(EBI + "/run", data=payload)
    log("  submitted job %s" % job)
    return job


def wait(job, poll=30):
    while True:
        status = http(EBI + "/status/" + job, timeout=120).strip()
        if status != "RUNNING":
            log("  status %s" % status)
            if status != "FINISHED":
                sys.exit("PSI-BLAST job %s ended as %s" % (job, status))
            return
        time.sleep(poll)


def parse_hits(job, out_dir):
    """Hits of the PSI-BLAST run, from the service's JSON result.

    Returns a list of (cluster_id, subject_length, percent_identity) in hit
    order, plus the run's metadata for the report.
    """
    raw = http(EBI + "/result/%s/json" % job, timeout=1800)
    path = os.path.join(out_dir, "psiblast.json")
    open(path, "w").write(raw)
    data = json.loads(raw)

    meta = {
        "program": "%s %s" % (data["program"], data["version"]),
        "database": data["dbs"],
        "db_num": data["db_num"],
        "expect_upper": data["expect_upper"],
        "matrix": data["matrix"],
        "gap_open": data["gap_open"],
        "gap_extend": data["gap_extend"],
        "filter": data["filter"],
        "alignments_requested": data["alignments"],
        "query_len": data["query_len"],
    }
    hits = []
    iterations = set()
    worst_evalue = 0.0
    for h in data["hits"]:
        iterations.add(h.get("hit_iter"))
        ident = max(x["hsp_identity"] for x in h["hit_hsps"])
        worst_evalue = max(worst_evalue, max(x["hsp_expect"] for x in h["hit_hsps"]))
        hits.append((h["hit_id"], h["hit_len"], ident))
    meta["hit_count"] = len(hits)
    meta["iterations_present"] = sorted(i for i in iterations if i)
    meta["worst_evalue"] = worst_evalue
    meta["hit_cap_reached"] = len(hits) >= int(data["alignments"])
    return hits, meta


def fetch_hit_sequences(job, out_dir, wanted):
    """Full subject sequences, from the run's own database snapshot.

    The service's `preselected_seq` result carries the hit sequences as
    searched, so they need no second retrieval. Anything it omits is pulled
    from the UniProt UniRef endpoint.
    """
    raw = http(EBI + "/result/%s/preselected_seq" % job, timeout=1800)
    open(os.path.join(out_dir, "psiblast_hit_sequences.fasta"), "w").write(raw)

    seqs, cur = {}, None
    for line in raw.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            cur = name.split(":", 1)[1] if ":" in name else name
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    log("   %d sequences from the run's own result" % len(seqs))

    missing = [c for c in wanted if c not in seqs]
    if missing:
        log("   fetching %d missing clusters from UniProt" % len(missing))
        for i in range(0, len(missing), 80):
            chunk = missing[i:i + 80]
            q = " OR ".join("id:%s" % c for c in chunk)
            url = "%s/uniref/stream?query=%s&format=fasta" % (UNIPROT, urllib.parse.quote(q))
            cur, parts = None, {}
            for line in http(url, timeout=600).splitlines():
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    parts[cur] = []
                elif cur:
                    parts[cur].append(line.strip())
            for k, v in parts.items():
                seqs[k] = "".join(v)
    return seqs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--accession", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--email", required=True)
    p.add_argument("--jobid", help="reuse an existing EBI job instead of submitting")
    p.add_argument("--workdir", default=".")
    p.add_argument("--report")
    args = p.parse_args()

    os.makedirs(args.workdir, exist_ok=True)

    log("1. query %s" % args.accession)
    query = fetch_query(args.accession)
    qlen = len(query)
    lo, hi = LEN_LO * qlen, LEN_HI * qlen
    log("   %d residues -> length window %.1f .. %.1f" % (qlen, lo, hi))

    log("2. PSI-BLAST vs UniRef90 (single iteration, e-value %s, max %s hits)" % (EVALUE, MAX_HITS))
    job = args.jobid or submit(query, args.email, args.prefix)
    wait(job)

    log("3. parsing hits")
    hits, meta = parse_hits(job, args.workdir)
    log("   hits reported: %d (cap %s, reached: %s)"
        % (meta["hit_count"], meta["alignments_requested"], meta["hit_cap_reached"]))
    log("   iterations present: %s | worst e-value: %s"
        % (meta["iterations_present"], meta["worst_evalue"]))
    json.dump(meta, open(os.path.join(args.workdir, "psiblast_meta.json"), "w"), indent=1)

    log("4. applying the paper's filters")
    kept, drop_len, drop_id, drop_self = [], 0, 0, 0
    for cid, hit_len, ident in hits:
        if cid == "UniRef90_%s" % args.accession:
            drop_self += 1
            continue
        if not lo <= hit_len <= hi:
            drop_len += 1
            continue
        if ident >= MAX_IDENTITY:
            drop_id += 1
            continue
        kept.append((cid, hit_len, ident))
    log("   query's own cluster removed: %d" % drop_self)
    log("   outside the %.1f-%.1f length window: %d" % (lo, hi, drop_len))
    log("   with >=%.0f%% identity to the query: %d" % (MAX_IDENTITY, drop_id))
    log("   homologs kept: %d" % len(kept))

    log("5. collecting hit sequences")
    seqs = fetch_hit_sequences(job, args.workdir, [c for c, _, _ in kept])
    final, no_seq, len_mismatch = [], 0, 0
    for cid, hit_len, ident in kept:
        s = seqs.get(cid)
        if not s:
            no_seq += 1
            continue
        if len(s) != hit_len:
            len_mismatch += 1
        final.append((cid, s, ident))
    if no_seq:
        log("   WARNING: %d kept clusters had no retrievable sequence" % no_seq)
    if len_mismatch:
        log("   WARNING: %d sequences differ in length from the BLAST record" % len_mismatch)
    log("   sequences available for the MSA: %d" % len(final))

    if args.report:
        with open(args.report, "w") as f:
            print("uniref90_cluster\tlength\tpercent_identity_to_query", file=f)
            for cid, s, ident in final:
                print("%s\t%d\t%.1f" % (cid, len(s), ident), file=f)

    records = [SeqRecord(Seq(query), id=args.prefix, description="")]
    records += [SeqRecord(Seq(s), id=cid, description="") for cid, s, _ in final]
    fasta_in = os.path.join(args.workdir, "%s_psiblast_hits.fasta" % args.prefix)
    SeqIO.write(records, fasta_in, "fasta")

    log("6. MUSCLE (default parameters, ClustalW output) on %d sequences" % len(records))
    muscle = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "epiSCANEER", "msa_exec", "linux", "muscle")
    clw_out = os.path.join(args.workdir, "%s_muscle.aln" % args.prefix)
    t0 = time.time()
    proc = subprocess.run([muscle, "-in", fasta_in, "-clwout", clw_out],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit("MUSCLE failed:\n%s" % proc.stderr[-3000:])
    log("   MUSCLE finished in %.1f min" % ((time.time() - t0) / 60))

    aligned = {r.id: r for r in AlignIO.read(clw_out, "clustal")}
    ordered = [aligned[args.prefix]] + [aligned[r.id] for r in records[1:] if r.id in aligned]
    for rec in ordered:
        rec.description = ""
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    AlignIO.write(MultipleSeqAlignment(ordered), args.out, "clustal")
    log("7. wrote %s (%d sequences, %d columns)"
        % (args.out, len(ordered), len(ordered[0].seq)))


if __name__ == "__main__":
    main()
