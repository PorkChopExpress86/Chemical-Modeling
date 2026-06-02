#!/usr/bin/env bash
# Unified launcher for the Chemical-Modeling project.
#
# Everything runs inside the conda `chemmodel` env (Psi4 + PCM, PySCF, PubChem,
# ASE, Streamlit) via micromamba. The package is imported straight from src/ on
# PYTHONPATH, so no install step is required to run.
#
# Usage:
#   ./run.sh solubility [args...]   # nitrogen-solubility CLI
#   ./run.sh app                    # Single-molecule / reaction simulator  (:8501)
#   ./run.sh sandbox                # Reaction sandbox + reactor engineering (:8502)
#   ./run.sh test [args...]         # pytest suite
#   ./run.sh shell                  # interactive shell in the env
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/home/specter/micromamba}"
export PYTHONPATH="$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
MICROMAMBA="${MICROMAMBA:-/home/specter/micromamba/micromamba}"
ENV_NAME="${CHEMMODEL_ENV:-chemmodel}"

run() { "$MICROMAMBA" run -n "$ENV_NAME" "$@"; }

cmd="${1:-help}"
shift || true

case "$cmd" in
  solubility)
    run python -m chemmodel.solubility "$@"
    ;;
  app|simulation)
    run streamlit run src/chemmodel/ui/simulation_app.py --server.port "${PORT:-8501}" "$@"
    ;;
  sandbox)
    run streamlit run src/chemmodel/ui/sandbox_app.py --server.port "${PORT:-8502}" "$@"
    ;;
  test)
    # stdlib unittest by default (no extra deps); `pip install -e .[dev]` adds pytest.
    if run python -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('pytest') else 1)"; then
      run python -m pytest -q "$@"
    else
      run python -m unittest discover -s tests -t "$REPO_ROOT" "$@"
    fi
    ;;
  shell)
    run bash "$@"
    ;;
  *)
    sed -n '2,16p' "${BASH_SOURCE[0]}"
    exit 0
    ;;
esac
