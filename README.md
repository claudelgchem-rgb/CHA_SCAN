# epiSCANEER
---

This repository contains the core computational components of **epiSCANEER**.
**epiSCANEER** identifies and recommends residue combinations that improve enzymatic activity by exploiting epistatic relationships inferred from coevolutionary networks.

**epiSCANEER** is based on **SCANEER**, incorporating parts of the original code. \
SCANEER github link: https://github.com/SBIlab/SCANEER \
Upstream epiSCANEER link: https://github.com/khsub87/epiSCANEER

This copy is a *runnable* version of the upstream code: the sources were ported
to Python 3 and the few blocking bugs were fixed (see **Changes made to run the
upstream code** at the bottom). The scoring formulas are unchanged.

## Requirements
---
- Python >= 3.9 (tested on 3.11)
- Java >= 1.8 (used by the covariation step, `covariance.algorithms.*`) - tested on OpenJDK 21
- Python dependencies are listed in `requirements.txt`

## Setup
---
```
git clone <repository-url>
cd CHA_SCAN
./setup.sh
```
`setup.sh` checks the Python and Java versions, creates a virtual environment in
`.venv/`, installs the dependencies and verifies that the package imports.

Manual equivalent:
```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

If Java is missing:
```
sudo apt-get install -y default-jre     # Ubuntu / Debian
brew install openjdk                    # macOS
```

## Usage
---
1. Prepare input files
+ Place input files in the ```input/``` directory.
+ Each input file should be a multiple sequence alignment (MSA) of the target enzyme.
+ Input files must be in CLUSTAL format (*.aln).
+ **The first sequence of the alignment is the query enzyme**; every score is reported against its residue numbering.
+ You may use the example multiple sequence alignments provided in the ```input/``` directory (`KGSADH.aln`, `LdhA.aln`).

2. Configure input/output paths (optional)
+ To change the input directory, modify ```input_path``` in ```run_epiSCANEER.py```.
+ To change the output directory, modify ```output_path``` in ```run_epiSCANEER.py```.
+ Both can also be given on the command line: ```python run_epiSCANEER.py <input_dir> <output_dir>```.

3. Run epiSCANEER
```
.venv/bin/python run_epiSCANEER.py
```
or, after ```source .venv/bin/activate```:
```
python run_epiSCANEER.py
```
A gene is skipped if its final ```output/<prefix>/<prefix>.txt``` already exists;
delete that directory to recompute.

Runtime for the bundled examples (250 sequences each) is roughly 1 minute per
enzyme on a modern laptop.

4. Output files
+ All results will be saved in the ```output/<prefix>/``` directory.
+ The output directory is automatically created if it does not exist.
+ The contents of each output file are as follows:
     + ```*.aln``` - A copy of the input multiple sequence alignment.
     + ```*.aln_cn``` - The filtered alignment (query-gap columns and columns with >= 20% gaps removed) used to calculate covarying strength.
     + ```*.coe_out_mcbasc``` - A text file containing the calculated covarying strengths of all combination of residue pairs.
     + ```*.cn``` - A text file containing the number of co-evolutionary relationships of residues.
     + ```*.coenet``` - A text file containing a residue-residue co-evolutionary network of the enzyme.
     + ```*.doubleSCI_mp``` - Raw epiSCI scores of every candidate residue substitution pair.
     + ```*.txt``` - A text file containing final epiSCI scores of mutations (best-scoring substitution per co-evolving residue pair, with a min-max normalised ```SCI_norm``` column).

## Helper scripts
---
Two scripts in `scripts/` prepare inputs and interpret outputs. They are
additions to this copy, not part of upstream epiSCANEER.

`scripts/build_msa.py` builds an input alignment for any UniProt protein:
it fetches the query, downloads its Pfam family from UniProtKB, scores every
member against the query by local alignment, filters on coverage/identity,
removes near-duplicates, aligns the survivors with the bundled MUSCLE and
writes CLUSTAL with the query first (as epiSCANEER requires).
```
python scripts/build_msa.py --accession P38567 --pfam PF01630 \
    --prefix PH20 --out input/PH20.aln --target 300 \
    --report input/PH20_homologs.tsv
```

`scripts/annotate_ph20.py` annotates the PH20 results with mature-protein
numbering, region (signal peptide / rHuPH20 / GPI tail), distance to the
catalytic Glu148 measured on PDB 9JUB, and disulfide / N-glycosylation
flags, and writes a shortlist restricted to the secreted enzyme.
```
python scripts/annotate_ph20.py
```

`input/PH20.aln` (301 sequences of glycoside hydrolase family 56, human
PH-20 as query) and `input/PH20_homologs.tsv` are the alignment produced
this way.

## Changes made to run the upstream code
---
The upstream sources are half-way through a Python 2 -> 3 migration and do not
execute as published. The following was fixed; no scoring formula was modified.

| File | Problem | Fix |
| --- | --- | --- |
| `epiSCANEER/__init__.py` | Python 2 `print` statement | `print(...)` |
| `epiSCANEER/msa.py` | `from __init__ import *` (Python 2 implicit relative import), `from Bio.Alphabet import ...` (removed in Biopython >= 1.78), mixed tabs/spaces (`TabError`) | relative import, dropped the unused `Bio.Alphabet` import, consistent indentation |
| `epiSCANEER/coe.py` | `from __init__ import ...`, Python 2 `print >> fo`, hard-coded `/usr/bin/java`, config paths relative to the current directory | relative import, `print(..., file=...)`, `JAVA_HOME`/`PATH` lookup, paths derived from the package location, `subprocess` call that reports Java failures instead of failing silently |
| `epiSCANEER/coe/Energetics.properties` | `HOME_DIRECTORY` points at the original author's machine, so the McBASC metric file is not found | a machine-specific copy is generated in `epiSCANEER/coe/generated_conf/` and put first on the Java classpath (the checked-in file is left untouched) |
| `epiSCANEER/iscalc.py` | `import coe`, Python 2 `print >> f` | relative import, `print(..., file=...)` |
| `epiSCANEER/epiSCANEER.py` | `from dhpylib import ...` (package does not exist), `return output_dicF` (typo -> `NameError`), `final_epiSCI_normalization()` used the undefined global `base_pth`, `dict.values()` sorted as a Python 2 list | correct package import, `output_dic`, `base_pth` passed as an argument, `list(...)` |
| `run_epiSCNEER.py` | filename typo, Python 2 `print` statement, `build_msa()` called with one argument instead of two, the MSA was never copied to the output directory although the rest of the code reads it from there | renamed to `run_epiSCANEER.py`, Python 3 syntax, correct call, input alignment copied into `output/<prefix>/` |
| `requirements.txt` | `biopython==1.77` cannot be installed on Python >= 3.10 and is incompatible with `numpy 2.x`; four unused packages pinned | minimum versions of the three packages the code actually imports |

Performance: `get_num()`/`get_frequency()` re-counted amino-acid pairs over the
whole alignment for every candidate mutation, which made a single enzyme take
many hours. The counts are now memoised per residue pair (`_PAIR_COUNT_CACHE`),
and the zero-padded `np.mean` in `get_multiSCI()` is computed as
`sum(SCI_list) / number_of_edges`, which is the same value. Outputs were
verified to be identical to the original formulas (max absolute difference
2.2e-16, i.e. floating point summation order).

Reproducibility: candidate amino acids were iterated over a Python `set`, whose
order changes between processes, so runs disagreed on which substitution was
reported for residue pairs whose scores are exactly tied (7 of 962 pairs for
KGSADH). The amino acids are now iterated in sorted order, which makes the
output byte-identical across runs without changing any score.
