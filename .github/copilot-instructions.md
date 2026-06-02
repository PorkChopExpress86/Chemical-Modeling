# Chemical Modeling — Copilot Instructions

## What this project is
The `chemmodel` package (under `src/`) bundles three things:
1. **Nitrogen solubility** in EO/PO — a stdlib-only Henry's-law model + CLI.
2. **Simulation app** — single-molecule energies/forces and catalytic-reaction
   PES scans over pluggable quantum-chemistry engines.
3. **Reaction sandbox app** — PubChem → ab-initio Langevin MD (Psi4, optional PCM)
   plus closed-form ideal-reactor (Batch/CSTR/PFR) design.

## Architecture (layered; domain never imports Streamlit)
- `chemmodel.constants` — physical constants / unit conversions. **Single source of
  truth**; never re-define `HARTREE_EV`, conversion factors, etc. locally.
- `chemmodel.solubility` — `model.py` (physics) + `cli.py` (argparse/printing). No
  third-party deps.
- `chemmodel.engines` — **Strategy + Registry**. `base.QuantumEngine` defines
  `available` and `single_point(atoms, *, basis, method, charge=0, with_forces=False, **opts) -> EngineResult`
  (energy in **eV**). Concrete: `Psi4Engine`, `PySCFEngine`, `CP2KEngine`,
  `DemoEngine`. Resolve via `get_engine(name)` / `available_engines()`; **never
  branch on the backend in callers**.
- `chemmodel.chemistry` — ASE `structures`, `pubchem` lookup (pure, no Streamlit
  cache), `presets` (molecule/catalyst/reaction libraries).
- `chemmodel.reactions` — engine-agnostic PES scans (`run_reaction_scan`,
  `decompose_energies`, `demo_pes`).
- `chemmodel.reactors` — closed-form kinetics (`batch_profile`, `tau_pfr`,
  `tau_cstr`, `arrhenius_k`, …).
- `chemmodel.md` — Langevin MD worker (`md_worker`, `MD_PROCESS_LOCK`); runs off
  the Streamlit thread, never calls `st.*`, streams to a queue.
- `chemmodel.ui` — thin Streamlit apps (`simulation_app`, `sandbox_app`), shared
  `theme`, and Plotly/Streamlit `viz`. UI files run as `__main__` under
  `streamlit run`, so they use **absolute** `chemmodel...` imports.

## Solubility conventions (US customary, explicit units)
- Pressure psia (psig inputs converted), temperature °R (= °F + 459.67), volume
  ft³, Henry `Hcc` in psia·ft³/mol with `c_N2 = P_N2 / Hcc`.
- Henry correlation is a least-squares fit of `ln(H) = a + b/T` to tabulated
  `HenryPoint`s on each `SolventParameters` in `SOLVENTS`.
- **Fail closed:** a solvent with no `henry_points` (PO) raises
  `MissingHenryParameterError`. Do not invent parameters from analogs — add a real
  literature/lab source with its citation in `henry_reference`.
- Add a solvent by appending a `SolventParameters` to `SOLVENTS` in
  `solubility/model.py` (molecular weight, density + reference, Henry points +
  reference, notes).

## Engines & quantum-chemistry runtime
- Psi4/PySCF/CP2K are **conda-only**; the apps fall back to demo mode when absent.
  Run real chemistry in the `chemmodel` conda env (see README / `run.sh`).
- Psi4 is a process-global singleton that litters the CWD. Keep it contained:
  per-call work happens in a temp workdir (`_isolated_workdir`), and
  `_isolate_psi4_exit_litter()` relocates its atexit `timer.dat`. Always preserve
  `no_reorient` / `no_com` / `symmetry c1` when feeding ASE positions to Psi4 so
  force ordering matches.
- Add an engine by subclassing `QuantumEngine` and registering it in
  `engines/__init__.py`.

## Workflows
- Setup: `./scripts/setup-env.sh`. Run: `./run.sh {solubility|app|sandbox|test|shell}`.
- Tests are `unittest` under `tests/`; engine cases use `skipTest` when the backend
  is unavailable, so they pass in any env.

## Common pitfalls
- Don't put compute in the UI layer or import `streamlit` from a domain module.
- Don't re-derive constants/conversions locally — import from `chemmodel.constants`.
- Don't let Henry's law run outside the dilute regime; keep N₂ partial pressures
  modest.
