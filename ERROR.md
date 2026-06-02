# Error Log

## Psi4 not installable via pip
- **Cause:** Psi4 ships through conda-forge only; also has no build for Python 3.12+ (the `.venv` is 3.14).
- **Fix:** micromamba `chemmodel` env (Python 3.11) with psi4 + pcmsolver + pyscf + pubchempy.
- **Avoid:** Run anything that imports `psi4` (or needs PCM/pubchempy) through the `chemmodel` env (`./run.sh ...`), not `.venv`. Code paths that touch these degrade to demo mode when the backend is absent.

## Wrong MAMBA_ROOT_PREFIX silently points micromamba at a missing env
- **Symptom:** `critical libmamba The given prefix does not exist: ".../micromamba/envs/envs/chemmodel"` — note the doubled `envs/envs`.
- **Cause:** Setting `MAMBA_ROOT_PREFIX=/home/specter/micromamba/envs` makes micromamba resolve `-n chemmodel` to `$PREFIX/envs/chemmodel`, i.e. `.../envs/envs/chemmodel`. The env actually lives at `/home/specter/micromamba/envs/chemmodel`.
- **Fix:** `MAMBA_ROOT_PREFIX=/home/specter/micromamba` (root prefix is the micromamba *home*, not its `envs/` dir). `run.sh` / `setup-env.sh` set this correctly; prefer them. (An older note said `psi4env` / the wrong prefix — superseded.)
- **Avoid:** Confirm with `micromamba env list`; the env name's path is `$ROOT/envs/<name>`.

## Psi4 litters the project directory (CWD), incl. timer.dat at process exit
- **Symptom:** `PEDRA.OUT*`, `cavity.off*`, `cavity.npz`, `psi4-calc.dat`, and especially `timer.dat` / `psi.<pid>.clean` appearing in the repo root.
- **Cause:** Psi4 is a process-global singleton that writes scratch to the *current* directory, not `PSI_SCRATCH`. `timer.dat` is written by its **atexit `finalize`**, after any per-call workdir has been restored — and even when psi4 was only imported to probe availability.
- **Fix:** In `chemmodel.engines.psi4_engine`: each call runs inside `_isolated_workdir()` (chdir into a temp dir, rmtree after); `_isolate_psi4_exit_litter()` registers an atexit handler *after* psi4's own (LIFO → runs first) that chdir's into throwaway scratch so the finalize litter lands there. It is armed on the availability probe too. All scratch patterns are git-ignored.
- **Avoid:** Never assume `PSI_SCRATCH` keeps the CWD clean; wrap psi4 entry points in the workdir helpers.

## PCM solvation is very slow (performance, not a crash)
- **Symptom:** With a PCM solvent selected, ab-initio MD jumped from ~0.3 s/step to ~22 s/step at HF/STO-3G.
- **Cause:** PCMSolver rebuilds the GePol cavity on every gradient evaluation.
- **Fix/Mitigation:** `chemmodel.ui.sandbox_app` shows a runtime estimate (×55 penalty for PCM) and warns before running; defaults are Vacuum + HF/STO-3G + small step counts.
- **Avoid:** Don't pair PCM with large molecules / many steps unless the long runtime is acceptable.

## Force/position mismatch risk during MD
- **Cause:** Psi4 reorients/recenters geometry by default, which scrambles gradient ordering relative to ASE's atom array → wrong forces in MD.
- **Fix:** Build the Psi4 geometry with `no_reorient`, `no_com`, and `symmetry c1` (see `make_psi4_calculator` / `Psi4Engine` in `chemmodel.engines.psi4_engine`).
- **Avoid:** Keep those keywords whenever feeding ASE positions to Psi4 for forces.

## ASE Langevin fixcm FutureWarning
- **Symptom:** `FutureWarning: fixcm=True ... deprecated since ASE 3.28.0`.
- **Cause:** Default `Langevin(fixcm=True)` will be removed; correct NVT sampling wants `fixcm=False` + `FixCom`.
- **Status:** Non-breaking warning; MD runs fine (worker in `chemmodel.md.langevin`). Revisit if upgrading ASE past the removal.
