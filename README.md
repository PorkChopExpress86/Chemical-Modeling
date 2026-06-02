# Chemical Modeling

A small chemical-engineering toolkit packaged as `chemmodel`:

1. **Nitrogen solubility** in ethylene oxide (EO) / propylene oxide (PO) — a
   stdlib-only, literature-backed Henry's-law model with a CLI.
2. **Molecular simulation app** — single-molecule energies/forces and catalytic
   reaction potential-energy-surface (PES) scans over pluggable quantum-chemistry
   engines (Psi4 / PySCF / CP2K, with a synthetic demo fallback).
3. **Reaction sandbox app** — PubChem → ab-initio Langevin molecular dynamics
   (Psi4, optional PCM solvation) plus closed-form ideal-reactor design
   (Batch / CSTR / PFR).

## Project layout

```
src/chemmodel/
  constants.py        physical constants & unit conversions (single source of truth)
  solubility/         N2-in-EO/PO Henry's-law model + CLI   (stdlib only)
  engines/            QC engines as a Strategy + Registry (Psi4/PySCF/CP2K/Demo)
  chemistry/          ASE structures, PubChem lookup, preset molecule libraries
  reactions/          engine-agnostic PES scans
  reactors/           closed-form Batch/CSTR/PFR kinetics
  md/                 ab-initio Langevin MD worker (thread-safe, Streamlit-free)
  ui/                 thin Streamlit layers + Plotly dashboards (simulation_app, sandbox_app)
tests/                unit tests (Psi4/PySCF cases skip when the engine is absent)
data/                 reactor_specs_template.csv
scripts/setup-env.sh  micromamba environment bootstrap
run.sh                unified launcher
```

Domain modules never import Streamlit; all UI lives under `chemmodel.ui`. Quantum
backends are selected by name from a registry, so callers never branch on the
engine — see `chemmodel/engines/__init__.py`.

## Environments

Two layers of dependency:

- **Pure model + reactor kinetics** need only NumPy/ASE/Streamlit/Plotly (pip).
- **Real quantum chemistry** (Psi4 + PCM solvation, PySCF, PubChem) is conda-only
  — Psi4 has no build for Python 3.12+. Use the `chemmodel` conda env.

Create the conda env (micromamba) and install the package into it:

```bash
./scripts/setup-env.sh
```

This solves `environment.yml` into a `chemmodel` env (Python 3.11, psi4 1.10.1,
pcmsolver, pyscf, pubchempy, …) and `pip install -e . --no-deps` into it. When a
backend is missing the apps degrade to a clearly-labelled demo mode.

## Running

`run.sh` runs everything inside the conda env with `PYTHONPATH=src`:

```bash
./run.sh solubility --solvent EO          # nitrogen-solubility CLI
./run.sh app                              # simulation app      → http://localhost:8501
./run.sh sandbox                          # reaction sandbox    → http://localhost:8502
./run.sh test                             # unit tests
./run.sh shell                            # shell inside the env
```

## Nitrogen-solubility model

Models dissolved-nitrogen concentration (`mol/ft3`) and liquid-phase mole
fraction at fixed nitrogen pressure, over a temperature sweep (default `50 psig`,
`68–140 °F`). Equilibrium only — no closed-vessel pressure rise, solvent vapor
pressure, or reaction behavior.

EO has literature-backed Henry data (LyondellBasell EO Product Stewardship
Guidance, Table A4: N₂ Henry constants in liquid EO at 32/77/122 °F; primary
trail Olson, *J. Chem. Eng. Data* 1977, 22, 326-329), fit as `ln(H) = a + b/T`.
PO **fails closed**: no defensible nitrogen-in-PO Henry data was found, so the
model raises rather than guess from analogs.

```bash
./run.sh solubility                        # both solvents (PO reports "unavailable")
./run.sh solubility --solvent EO
./run.sh solubility --pressure-psig 50 --t-min-f 68 --t-max-f 140 --t-step-f 12
```

Python API:

```python
from chemmodel.solubility import nitrogen_solubility, temperature_sweep

single = nitrogen_solubility("EO", temperature_f=68.0, pressure_psig=50.0)
print(single.dissolved_mol_per_ft3, single.liquid_mole_fraction)

for row in temperature_sweep("EO", pressure_psig=50.0):
    print(row.temperature_f, row.dissolved_mol_per_ft3, row.liquid_mole_fraction)
```

## Tests

```bash
./run.sh test                              # conda env: exercises real Psi4/PySCF too
PYTHONPATH=src python -m unittest discover -s tests -t .   # any env (psi4 cases skip)
```

## Adding a quantum-chemistry engine

Subclass `chemmodel.engines.base.QuantumEngine` (implement `available` and
`single_point`), then register the instance in `chemmodel/engines/__init__.py`.
Nothing else changes — the apps and PES scans resolve it by name.

## Next steps (solubility)

- Add direct PO Henry-law parameters when a defensible source is available.
- Temperature-dependent density for wider-range accuracy.
- Uncertainty bounds around the Henry-law fit for design screening.
