"""Quantum-chemistry engine registry (Strategy + Registry pattern).

Callers resolve a backend by name and ask it for availability or a single point;
they never import a concrete engine or branch on the backend themselves::

    from chemmodel.engines import get_engine, available_engines

    engine = get_engine("pyscf")
    if engine.available:
        result = engine.single_point(atoms, basis="sto-3g", method="HF")
"""

from __future__ import annotations

from typing import Dict

from chemmodel.engines.base import EngineResult, QuantumEngine
from chemmodel.engines.cp2k_engine import CP2KEngine
from chemmodel.engines.demo import DemoEngine
from chemmodel.engines.psi4_engine import Psi4Engine, make_psi4_calculator
from chemmodel.engines.pyscf_engine import PySCFEngine

# Single shared instance per backend (engines are stateless strategies).
_REGISTRY: Dict[str, QuantumEngine] = {
    engine.name: engine
    for engine in (Psi4Engine(), PySCFEngine(), CP2KEngine(), DemoEngine())
}


def get_engine(name: str) -> QuantumEngine:
    """Return the engine registered under ``name`` (case-insensitive)."""
    try:
        return _REGISTRY[name.strip().lower()]
    except KeyError as exc:
        known = ", ".join(sorted(_REGISTRY))
        raise KeyError(f"Unknown engine {name!r}. Registered engines: {known}") from exc


def available_engines() -> Dict[str, bool]:
    """Map each registered engine name to whether its dependency is importable."""
    return {name: engine.available for name, engine in _REGISTRY.items()}


def probe(module: str) -> bool:
    """True if ``module`` can be imported — used for optional, non-engine deps."""
    try:
        __import__(module)
        return True
    except ImportError:
        return False


__all__ = [
    "EngineResult",
    "QuantumEngine",
    "Psi4Engine",
    "PySCFEngine",
    "CP2KEngine",
    "DemoEngine",
    "get_engine",
    "available_engines",
    "probe",
    "make_psi4_calculator",
]
