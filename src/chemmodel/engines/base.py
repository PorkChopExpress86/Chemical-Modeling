"""Quantum-chemistry engine Strategy interface.

Every backend (Psi4, PySCF, CP2K, and a synthetic Demo) implements
:class:`QuantumEngine`, so callers select a strategy by name from the registry
(:mod:`chemmodel.engines`) and never branch on the backend themselves.

Energy convention: engines return energy in **eV** (converting from Hartree at
their boundary). Forces, when requested, are in eV/Å with the ASE sign
convention (force = -dE/dx).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
from ase import Atoms


@dataclass(frozen=True)
class EngineResult:
    """Outcome of a single-point evaluation.

    ``forces`` is ``None`` when forces were not requested (energy-only, e.g. a
    potential-energy-surface scan); otherwise an ``(n_atoms, 3)`` eV/Å array.
    """

    energy_eV: float
    forces: Optional[np.ndarray]
    note: str


class QuantumEngine(ABC):
    """A pluggable quantum-chemistry backend."""

    #: short registry key, e.g. ``"psi4"``
    name: str = "base"
    #: human-facing label for UI badges
    label: str = "Quantum engine"

    @property
    @abstractmethod
    def available(self) -> bool:
        """True when this backend's Python/binary dependency is importable."""

    @abstractmethod
    def single_point(
        self,
        atoms: Atoms,
        *,
        basis: str,
        method: str,
        charge: int = 0,
        with_forces: bool = False,
        **options: Any,
    ) -> EngineResult:
        """Compute a single-point energy (eV), optionally with forces.

        ``options`` carries backend-specific knobs (e.g. ``max_iter`` for Psi4;
        ``use_pbc``/``cell_size``/``kpts`` for PySCF; ``cutoff``/``rel_cutoff``
        for CP2K). Engines read what they recognise and ignore the rest.
        """

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<{type(self).__name__} name={self.name!r} available={self.available}>"
