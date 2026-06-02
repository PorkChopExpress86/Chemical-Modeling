# Project Memory

## Current State
Restructured into an installable `chemmodel` package (src layout). Three products:
- **Nitrogen solubility** model (`chemmodel.solubility`) for N₂ in EO/PO, US
  customary units, stdlib-only. EO literature-backed; PO fails closed.
- **Simulation app** (`chemmodel.ui.simulation_app`, port 8501): single-molecule
  energies/forces + catalytic-reaction PES scans over the engine registry.
- **Reaction sandbox** (`chemmodel.ui.sandbox_app`, port 8502): 🧪 Molecular
  Sandbox (PubChem + ASE + Psi4-PCM + Langevin MD) and 🏭 Reactor Engineering
  (closed-form Batch/CSTR/PFR, Arrhenius, Levenspiel, adiabatic ΔT).

## Architecture (design pattern)
- **Strategy + Registry** for QC engines: `chemmodel.engines` (`QuantumEngine`
  base; `Psi4Engine`/`PySCFEngine`/`CP2KEngine`/`DemoEngine`; `get_engine`,
  `available_engines`). Callers select by name and never branch on the backend.
- **Layered**: pure domain (`engines`, `chemistry`, `reactions`, `reactors`, `md`,
  `solubility`) imports no Streamlit; all UI in `chemmodel.ui`.
- **Single source of truth** for constants/conversions in `chemmodel.constants`.
- Streamlit pages are thin: widget wiring + `ui.viz` Plotly dashboards only.

## Environments
- `.venv` (Python 3.14): pip deps (streamlit, ase, pyscf, plotly, numpy). **No
  psi4 / pubchempy** here — apps run in demo/analytic fallback.
- Conda **`chemmodel`** env via micromamba (Python 3.11): full stack incl.
  psi4 1.10.1 (+ pcmsolver), pyscf, pubchempy, pandas, scipy.
- **Correct micromamba invocation** (the path matters — see ERROR.md):
  `MAMBA_ROOT_PREFIX=/home/specter/micromamba /home/specter/micromamba/micromamba run -n chemmodel ...`
  The `chemmodel` env lives at `/home/specter/micromamba/envs/chemmodel`.
- Just use `./run.sh {solubility|app|sandbox|test|shell}` and
  `./scripts/setup-env.sh` — they set the prefix and `PYTHONPATH=src`.

## Key Decisions
- Full restructure (relocate compute into the package; apps became thin UI). All
  domain logic is now unit-tested (39 tests; psi4/pyscf cases `skipTest` when the
  backend is absent so the suite passes in either env).
- Engines return energy in **eV** (convert from Hartree at their boundary); the PES
  scan keeps its own far-separation reference subtraction.
- Psi4 CWD litter contained at the source: per-call temp workdir + an atexit-LIFO
  handler that relocates `timer.dat` (see `engines/psi4_engine.py`). Scratch
  patterns are git-ignored.
- `fetch_pubchem_xyz` is pure in `chemmodel.chemistry.pubchem`; the UI wraps it
  with `st.cache_data`. MD worker stays Streamlit-free and single-flight via
  `MD_PROCESS_LOCK`.
- Tests run via stdlib `unittest` (no pytest dependency); `run.sh test` uses pytest
  only if `pip install -e .[dev]` added it.

## Verification done
- `./run.sh test` → 39 passing in the conda env (real Psi4 + PySCF exercised).
- Streamlit `AppTest` smoke: both page scripts execute with 0 exceptions.
- `./run.sh solubility --solvent EO` prints the sweep; PO fails closed.
- Repo root stays clean after Psi4 runs.
