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
    """Return {cluster_id: best percent identity to the query}, in hit order."""
    types = http(EBI + "/resulttypes/" + job)
    available = re.findall(r"<identifier>([^<]+)</identifier>", types)
    log("  result types: %s" % ", ".join(available))

    raw = http(EBI + "/result/%s/out" % job, timeout=600)
    open(os.path.join(out_dir, "psiblast_out.txt"), "w").write(raw)

    ids_text = ""
    if "ids" in available:
        ids_text = http(EBI + "/result/%s/ids" % job, timeout=600)
        open(os.path.join(out_dir, "psiblast_ids.txt"), "w").write(ids_text)

    # Identities come from the tab-delimited or textual report; both are kept
    # on disk so the parse can be audited.
    hits = {}
    order = []
    # "UniRef90_XXXX ... <bits> <evalue> ... <identity>%" lines of the
    # NCBI-style summary table are parsed generically: take every cluster id
    # seen in the pairwise section together with its "Identities = a/b (c%)".
    blocks = re.split(r"\n(?=>)", raw)
    for block in blocks[1:]:
        m = re.match(r">\s*(UniRef90_\S+)", block)
        if not m:
            continue
        cid = m.group(1)
        idents = [float(x) for x in re.findall(r"Identities\s*=\s*\d+/\d+\s*\((\d+)%\)", block)]
        if not idents:
            continue
        if cid not in hits:
            order.append(cid)
        hits[cid] = max(hits.get(cid, 0.0), max(idents))
    return hits, order


def fetch_uniref_sequences(cluster_ids, batch=80):
    """Representative sequence of each UniRef90 cluster."""
    seqs = {}
    for i in range(0, len(cluster_ids), batch):
        chunk = cluster_ids[i:i + batch]
        q = " OR ".join("id:%s" % c for c in chunk)
        url = "%s/uniref/stream?query=%s&format=fasta" % (UNIPROT, urllib.parse.quote(q))
        text = http(url, timeout=600)
        cur = None
        for line in text.splitlines():
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
        log("  fetched %d / %d clusters" % (min(i + batch, len(cluster_ids)), len(cluster_ids)))
    return {k: "".join(v) for k, v in seqs.items()}


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
    idents, order = parse_hits(job, args.workdir)
    log("   hits reported: %d" % len(order))

    log("4. fetching UniRef90 representative sequences")
    seqs = fetch_uniref_sequences(order)
    missing = [c for c in order if c not in seqs]
    if missing:
        log("   WARNING: %d clusters returned no sequence (%s...)" % (missing.__len__(), missing[:3]))

    log("5. applying the paper's filters")
    kept, drop_len, drop_id = [], 0, 0
    for cid in order:
        s = seqs.get(cid)
        if not s:
            continue
        if cid == "UniRef90_%s" % args.accession:
            continue                      # the query's own cluster
        if not lo <= len(s) <= hi:
            drop_len += 1
            continue
        if idents[cid] >= MAX_IDENTITY:
            drop_id += 1
            continue
        kept.append((cid, s, idents[cid]))
    log("   dropped outside length window: %d" % drop_len)
    log("   dropped with >=%.0f%% identity to query: %d" % (MAX_IDENTITY, drop_id))
    log("   homologs kept: %d" % len(kept))

    if args.report:
        with open(args.report, "w") as f:
            print("uniref90_cluster\tlength\tpercent_identity_to_query", file=f)
            for cid, s, ident in kept:
                print("%s\t%d\t%.0f" % (cid, len(s), ident), file=f)

    records = [SeqRecord(Seq(query), id=args.prefix, description="")]
    records += [SeqRecord(Seq(s), id=cid, description="") for cid, s, _ in kept]
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
