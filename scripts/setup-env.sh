#!/usr/bin/env bash
# Create (or update) the conda `chemmodel` environment with micromamba and install
# the package into it in editable mode. Idempotent — safe to re-run.
#
# Usage:  ./scripts/setup-env.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/home/specter/micromamba}"
MICROMAMBA="${MICROMAMBA:-/home/specter/micromamba/micromamba}"
ENV_NAME="${CHEMMODEL_ENV:-chemmodel}"

if [ ! -x "$MICROMAMBA" ]; then
  echo "micromamba not found at $MICROMAMBA — set MICROMAMBA to its path." >&2
  exit 1
fi

echo "→ Solving environment '$ENV_NAME' from environment.yml (this can take a few minutes)…"
if "$MICROMAMBA" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  "$MICROMAMBA" install -y -n "$ENV_NAME" -f environment.yml
else
  "$MICROMAMBA" create -y -n "$ENV_NAME" -f environment.yml
fi

echo "→ Installing chemmodel into '$ENV_NAME' (editable, no deps — conda owns the stack)…"
"$MICROMAMBA" run -n "$ENV_NAME" python -m pip install -e . --no-deps

echo "→ Smoke check:"
"$MICROMAMBA" run -n "$ENV_NAME" python -c "import chemmodel; print('chemmodel', chemmodel.__version__, 'ready')"

cat <<'EOF'

Done. Try:
  ./run.sh test
  ./run.sh solubility --solvent EO
  ./run.sh app        # http://localhost:8501
  ./run.sh sandbox    # http://localhost:8502
EOF
