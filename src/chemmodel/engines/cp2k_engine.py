"""CP2K backend — Quickstep DFT forces via the ASE CP2K calculator."""

from __future__ import annotations

from typing import Any

from ase import Atoms

from chemmodel.engines.base import EngineResult, QuantumEngine


class CP2KEngine(QuantumEngine):
    name = "cp2k"
    label = "CP2K — MD Forces"

    @property
    def available(self) -> bool:
        try:
            import cp2k  # type: ignore  # noqa: F401

            return True
        except ImportError:
            return False

    def single_point(
        self,
        atoms: Atoms,
        *,
        basis: str = "",
        method: str = "",
        charge: int = 0,
        with_forces: bool = False,
        cutoff: float = 400.0,
        rel_cutoff: float = 60.0,
        **options: Any,
    ) -> EngineResult:
        from ase.calculators.cp2k import CP2K  # type: ignore

        inp = (
            f"&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    &MGRID\n"
            f"      CUTOFF {cutoff}\n      REL_CUTOFF {rel_cutoff}\n"
            f"    &END MGRID\n  &END DFT\n&END FORCE_EVAL\n"
        )
        atoms.calc = CP2K(inp=inp, auto_write=True)
        forces = atoms.get_forces() if with_forces else None
        return EngineResult(
            atoms.get_potential_energy(),
            forces,
            f"CP2K cutoff {cutoff} Ry",
        )
