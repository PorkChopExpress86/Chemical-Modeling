"""Psi4 backend — molecular single points and a PCM-solvated ASE calculator.

Psi4 is a process-global singleton that scatters scratch files (``PEDRA.OUT``,
``cavity.*``, ``psi4-calc.dat``, ``timer.dat``) into the *current working
directory*. Every entry point here runs inside a private temp workdir so that
litter never lands in the repo root.
"""

from __future__ import annotations

import atexit
import contextlib
import os
import shutil
import tempfile
from typing import Any

import numpy as np
from ase import Atoms

from chemmodel.constants import HA_BOHR_TO_EV_ANG, HARTREE_EV
from chemmodel.engines.base import EngineResult, QuantumEngine

# PubChem/UI solvent name → PCMSolver named solvent
PCM_SOLVENTS = {"Water": "Water", "Ethanol": "Ethanol", "Benzene": "Benzene"}

_PSI4_EXIT_ISOLATED = False


def _isolate_psi4_exit_litter() -> None:
    """Keep Psi4's process-finalize litter out of the project directory.

    Psi4 writes ``timer.dat`` / ``psi.<pid>.clean`` to the *current* directory when
    its atexit ``finalize`` runs — after any per-call workdir has been restored.
    Registered here (always *after* Psi4's own atexit handler, so it runs first
    under atexit's LIFO order), this moves the cwd into throwaway scratch so the
    finalize litter lands there instead of the repo root. Call once, after Psi4 is
    imported.
    """
    global _PSI4_EXIT_ISOLATED
    if _PSI4_EXIT_ISOLATED:
        return
    atexit.register(os.chdir, tempfile.mkdtemp(prefix="psi4_exit_"))
    _PSI4_EXIT_ISOLATED = True


@contextlib.contextmanager
def _isolated_workdir(prefix: str = "psi4_scratch_"):
    """Run Psi4 inside a throwaway directory, restoring cwd afterwards."""
    workdir = tempfile.mkdtemp(prefix=prefix)
    prev_cwd = os.getcwd()
    try:
        os.environ["PSI_SCRATCH"] = workdir
        os.chdir(workdir)
        yield workdir
    finally:
        os.chdir(prev_cwd)
        shutil.rmtree(workdir, ignore_errors=True)


class Psi4Engine(QuantumEngine):
    name = "psi4"
    label = "Psi4 — Quantum Chemistry"

    @property
    def available(self) -> bool:
        try:
            import psi4  # type: ignore  # noqa: F401
        except ImportError:
            return False
        # importing psi4 alone arms its atexit finalize (which writes timer.dat) —
        # neutralise it even when callers only probe availability.
        _isolate_psi4_exit_litter()
        return True

    def single_point(
        self,
        atoms: Atoms,
        *,
        basis: str,
        method: str,
        charge: int = 0,
        with_forces: bool = False,
        max_iter: int = 150,
        **options: Any,
    ) -> EngineResult:
        if with_forces:
            return self._with_forces(atoms, basis, method, max_iter)
        return self._energy_only(atoms, basis, method, charge)

    # ── energy-only (charge/multiplicity aware) — used by PES scans ────────────
    def _energy_only(self, atoms, basis, method, charge) -> EngineResult:
        import psi4  # type: ignore

        _isolate_psi4_exit_litter()
        with _isolated_workdir():
            psi4.core.be_quiet()
            psi4.set_memory("1 GB")
            psi4.set_num_threads(os.cpu_count() or 1)

            n_elec = sum(atoms.get_atomic_numbers()) - charge
            mult = 1 + (n_elec % 2)
            geom_lines = [
                f"{s}  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}"
                for s, p in zip(atoms.get_chemical_symbols(), atoms.get_positions())
            ]
            geom_str = "\n".join(
                [f"{charge} {mult}", *geom_lines, "units angstrom", "no_reorient", "symmetry c1"]
            )
            mol = psi4.geometry(geom_str)
            energy_ha = float(psi4.energy(f"{method}/{basis}", molecule=mol))
        return EngineResult(energy_ha * HARTREE_EV, None, f"Psi4 {method}/{basis}")

    # ── single point with forces (ASE Psi4 calculator) ────────────────────────
    def _with_forces(self, atoms, basis, method, max_iter) -> EngineResult:
        from ase.calculators.psi4 import Psi4  # type: ignore

        _isolate_psi4_exit_litter()
        with _isolated_workdir() as scratch:
            calc = Psi4(
                atoms=atoms,
                method=method,
                basis=basis,
                maxiter=max_iter,
                num_threads=max(1, os.cpu_count() or 1),
                memory="1 GB",
                PSI_SCRATCH=scratch,
                symmetry="c1",
            )
            atoms.calc = calc
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
        return EngineResult(energy, forces, f"Psi4 {method}/{basis}")


def make_psi4_calculator(method, basis, solvent, charge, multiplicity, mem="2 GB", nthreads=None):
    """Build an ASE calculator wrapping Psi4 with optional PCM solvation.

    Returned calculator keeps all of Psi4's CWD litter inside its own per-instance
    temp workdir and resets the global singleton in a ``finally`` so an interrupted
    gradient cannot poison the next step. Call ``.cleanup()`` when done.
    """
    import psi4
    from ase.calculators.calculator import Calculator, all_changes

    _isolate_psi4_exit_litter()
    if nthreads is None:
        nthreads = os.cpu_count() or 1

    class Psi4PCM(Calculator):
        """ASE calculator wrapping Psi4 with optional PCM solvation.

        Psi4 writes PCM cavity caches (cavity.npz, PEDRA.OUT, cavity.off) into the
        *current working directory* — not PSI_SCRATCH. So each instance gets its
        own temp workdir and chdir's into it only for the duration of a gradient,
        then restores. The global state (options/variables/active SCF) is reset in
        a finally so an interrupted step cannot poison the next run. Energy and
        gradient are read BEFORE that reset (clean() wipes them).
        """

        implemented_properties = ["energy", "forces"]

        def __init__(self):
            super().__init__()
            self._workdir = tempfile.mkdtemp(prefix="psi4_run_")
            os.environ["PSI_SCRATCH"] = self._workdir
            psi4.core.IOManager.shared_object().set_default_path(self._workdir)
            psi4.core.be_quiet()
            psi4.set_memory(mem)
            psi4.set_num_threads(nthreads)

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            lines = [f"{charge} {multiplicity}"]
            for s, p in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
                lines.append(f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}")
            lines += ["units angstrom", "symmetry c1", "no_reorient", "no_com"]
            use_pcm = solvent in PCM_SOLVENTS

            prev_cwd = os.getcwd()
            try:
                # keep all Psi4 CWD litter inside the per-instance workdir
                os.chdir(self._workdir)
                psi4.core.clean_options()
                opts = {"basis": basis, "scf_type": "pk"}
                if use_pcm:
                    opts.update({"pcm": True, "pcm_scf_type": "total"})
                psi4.set_options(opts)
                mol = psi4.geometry("\n".join(lines))
                if use_pcm:
                    # a truncated cavity cache from an interrupted run would
                    # poison this gradient — drop it first
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(os.path.join(self._workdir, "cavity.npz"))
                    psi4.pcm_helper(f"""
Units = Angstrom
Medium {{
  SolverType = IEFPCM
  Solvent = {PCM_SOLVENTS[solvent]}
}}
Cavity {{
  Type = GePol
  Area = 0.3
}}
""")
                grad = psi4.gradient(method, molecule=mol)
                energy = psi4.variable("CURRENT ENERGY")  # read BEFORE clean()
                g = np.asarray(grad)  # Hartree/Bohr
                self.results["energy"] = energy * HARTREE_EV
                self.results["forces"] = -g * HA_BOHR_TO_EV_ANG
            finally:
                # reset the global singleton so the next step / next run starts
                # clean even if this one was interrupted
                with contextlib.suppress(Exception):
                    psi4.core.clean()
                    psi4.core.clean_options()
                    psi4.core.clean_variables()
                os.chdir(prev_cwd)

        def cleanup(self):
            shutil.rmtree(self._workdir, ignore_errors=True)

    return Psi4PCM()
