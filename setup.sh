#!/usr/bin/env bash
# epiSCANEER setup: creates a local virtual environment and checks the
# external tools the pipeline needs (Java for the covariation step).
set -euo pipefail

cd "$(dirname "$0")"

PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-.venv}"

echo "==> Checking Python"
"$PYTHON_BIN" --version
"$PYTHON_BIN" - <<'PY'
import sys
if sys.version_info < (3, 9):
    sys.exit("Python 3.9 or newer is required (found %s)" % sys.version.split()[0])
PY

echo "==> Checking Java (required by covariance.algorithms.*)"
if command -v java >/dev/null 2>&1; then
    java -version 2>&1 | head -1
elif [ -n "${JAVA_HOME:-}" ] && [ -x "$JAVA_HOME/bin/java" ]; then
    "$JAVA_HOME/bin/java" -version 2>&1 | head -1
else
    echo "ERROR: no java found. Install a JRE/JDK (>= 1.8), e.g.:" >&2
    echo "  Ubuntu/Debian : sudo apt-get install -y default-jre" >&2
    echo "  macOS (brew)  : brew install openjdk" >&2
    exit 1
fi

echo "==> Creating virtual environment in $VENV_DIR"
if [ ! -d "$VENV_DIR" ]; then
    "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

echo "==> Installing Python dependencies"
"$VENV_DIR/bin/python" -m pip install --upgrade pip >/dev/null
"$VENV_DIR/bin/python" -m pip install -r requirements.txt

echo "==> Verifying imports"
"$VENV_DIR/bin/python" - <<'PY'
import Bio, networkx, numpy
from epiSCANEER import msa, coe, iscalc
from epiSCANEER.epiSCANEER import calc_epiSCI
print("biopython %s / networkx %s / numpy %s" % (Bio.__version__, networkx.__version__, numpy.__version__))
print("epiSCANEER package imports OK")
PY

cat <<'MSG'

Setup complete. Run the pipeline with:

    .venv/bin/python run_epiSCANEER.py

or activate the environment first:

    source .venv/bin/activate
    python run_epiSCANEER.py

Results are written to ./output/<prefix>/.
MSG
