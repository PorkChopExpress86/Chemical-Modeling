"""Synthetic engine — always available, never touches a real QC backend.

Used as the graceful fallback when no real engine is importable, so the apps and
tests still produce plausible (clearly-labelled) numbers.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from ase import Atoms

from chemmodel.engines.base import EngineResult, QuantumEngine


class DemoEngine(QuantumEngine):
    name = "demo"
    label = "Demo (synthetic)"

    @property
    def available(self) -> bool:
        return True

    def single_point(
        self,
        atoms: Atoms,
        *,
        basis: str = "",
        method: str = "",
        charge: int = 0,
        with_forces: bool = False,
        label: str = "Demo",
        **options: Any,
    ) -> EngineResult:
        rng = np.random.default_rng(42)
        n = len(atoms)
        energy_eV = -n * 5.0 + rng.normal(0, 0.5)
        forces = rng.normal(0, 0.3, (n, 3)) if with_forces else None
        return EngineResult(
            energy_eV=float(energy_eV),
            forces=forces,
            note=f"⚠️ {label} — synthetic demo data",
        )
