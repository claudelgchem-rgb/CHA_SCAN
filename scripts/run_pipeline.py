#!/usr/bin/env python3
"""Stage runner that survives this container being reclaimed.

Long epiSCANEER runs keep being lost: the container is reclaimed whenever
the session goes idle, and it restores only what git tracks. This runner
treats each expensive stage as a checkpoint - it skips a stage whose output
is already present and, with --commit, commits and pushes that output as
soon as it is produced, so an interrupted run resumes instead of restarting.

Stages: msa -> episcaneer -> annotate -> report

    python scripts/run_pipeline.py --prefix PH20full \
        --aln input/PH20_full.aln --outdir results/PH20_full --commit

Each stage declares the files that prove it finished. A stage whose outputs
all exist is reported as "cached" and not re-run; --force re-runs it.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(REPO, ".venv", "bin", "python")
EPI_OUTPUTS = ["PH20.aln_cn", "PH20.cn", "PH20.coe_out_mcbasc",
               "PH20.coenet", "PH20.doubleSCI_mp", "PH20.txt"]


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def have(paths):
    return all(os.path.isfile(p) and os.path.getsize(p) > 0 for p in paths)


def git(*args):
    return subprocess.run(["git"] + list(args), cwd=REPO,
                          capture_output=True, text=True)


def checkpoint(paths, message, push=True):
    """Commit the stage's outputs so a container restart cannot lose them."""
    existing = [os.path.relpath(p, REPO) for p in paths if os.path.exists(p)]
    if not existing:
        return
    git("add", "--", *existing)
    if not git("diff", "--cached", "--quiet").returncode:
        log("  checkpoint: nothing new to commit")
        return
    r = git("commit", "-q", "-m", message)
    if r.returncode:
        log("  checkpoint: commit failed: %s" % (r.stderr or r.stdout).strip()[:300])
        return
    log("  checkpoint: committed %d path(s)" % len(existing))
    if push:
        for attempt in range(4):
            r = git("push", "-u", "origin", "HEAD")
            if r.returncode == 0:
                log("  checkpoint: pushed")
                return
            time.sleep(2 ** attempt)
        log("  checkpoint: push failed, commit is local only")


def run(cmd, logfile):
    log("  $ %s" % " ".join(cmd[:4] + (["..."] if len(cmd) > 4 else [])))
    with open(logfile, "wb") as f:
        p = subprocess.run(cmd, cwd=REPO, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        sys.exit("stage failed (rc=%d); see %s" % (p.returncode, logfile))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="name used in logs and commits")
    ap.add_argument("--aln", required=True, help="CLUSTAL alignment, query first")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--structure", default="9jub.cif")
    ap.add_argument("--structure-url", default="https://files.rcsb.org/download/9JUB.cif")
    ap.add_argument("--workdir", default=None, help="scratch for the epiSCANEER run")
    ap.add_argument("--commit", action="store_true", help="checkpoint each stage into git")
    ap.add_argument("--no-push", action="store_true")
    ap.add_argument("--force", nargs="*", default=[],
                    help="stage names to re-run even if their outputs exist")
    args = ap.parse_args()

    outdir = os.path.join(REPO, args.outdir)
    os.makedirs(outdir, exist_ok=True)
    workdir = args.workdir or os.path.join(outdir, ".work")
    os.makedirs(os.path.join(workdir, "input"), exist_ok=True)
    push = args.commit and not args.no_push

    # --- stage: episcaneer -------------------------------------------------
    epi_paths = [os.path.join(outdir, f) for f in EPI_OUTPUTS]
    if have(epi_paths) and "episcaneer" not in args.force:
        log("episcaneer: cached (%d outputs present)" % len(epi_paths))
    else:
        log("episcaneer: running - this is the long stage")
        staged = os.path.join(workdir, "input", "PH20.aln")
        if not os.path.isfile(staged) or os.path.getsize(staged) != os.path.getsize(args.aln):
            shutil.copyfile(args.aln, staged)
        t0 = time.time()
        run([PY, os.path.join(REPO, "run_epiSCANEER.py"),
             os.path.join(workdir, "input"), os.path.join(workdir, "output")],
            os.path.join(workdir, "episcaneer.log"))
        produced = os.path.join(workdir, "output", "PH20")
        for f in EPI_OUTPUTS:
            shutil.copyfile(os.path.join(produced, f), os.path.join(outdir, f))
        log("episcaneer: done in %.1f min" % ((time.time() - t0) / 60))
        if args.commit:
            checkpoint(epi_paths, "%s: epiSCANEER outputs" % args.prefix, push)

    # --- stage: annotate ---------------------------------------------------
    ann = [os.path.join(outdir, f) for f in ("PH20.annotated.tsv", "PH20.shortlist.tsv")]
    if have(ann) and "annotate" not in args.force:
        log("annotate: cached")
    else:
        log("annotate: running")
        struct = args.structure
        if not os.path.isfile(struct):
            struct = os.path.join(workdir, "9jub.cif")
            if not os.path.isfile(struct):
                import urllib.request
                urllib.request.urlretrieve(args.structure_url, struct)
        run([PY, os.path.join(REPO, "scripts", "annotate_ph20.py"),
             "--results", os.path.join(outdir, "PH20.txt"),
             "--structure", struct,
             "--out", ann[0], "--shortlist", ann[1]],
            os.path.join(workdir, "annotate.log"))
        log("annotate: done")
        if args.commit:
            checkpoint(ann, "%s: structural annotation" % args.prefix, push)

    log("pipeline complete -> %s" % args.outdir)


if __name__ == "__main__":
    main()
