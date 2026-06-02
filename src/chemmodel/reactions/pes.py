"""Potential-energy-surface scans for a bimolecular (optionally catalysed) system.

Engine-agnostic: the backend is selected by name from :mod:`chemmodel.engines`.
When the chosen engine is unavailable, the scan falls back to a synthetic Morse
curve (``demo_pes``) so the UI still has something to plot.
"""

from __future__ import annotations

import numpy as np
from ase import Atoms

from chemmodel.chemistry.structures import build_reaction_system
from chemmodel.engines import get_engine


def demo_pes(distances: np.ndarray, seed: int = 7) -> np.ndarray:
    """Morse-like relative potential (eV) for demo mode."""
    rng = np.random.default_rng(seed)
    d0, De, a = 2.2, 0.15, 1.2
    morse = De * (1 - np.exp(-a * (distances - d0))) ** 2 - De
    noise = rng.normal(0, 0.002, len(distances))
    return morse + noise  # eV, relative


def run_reaction_scan(
    atoms_a: Atoms,
    atoms_b: Atoms,
    catalyst_atoms: Atoms | None,
    catalyst_charge: int,
    distances: np.ndarray,
    engine: str,
    basis: str,
    method: str,
    progress_cb=None,
) -> list[dict]:
    """For each separation d, build (A + B [+ catalyst]) and compute a single point.

    Returns a list of scan-point dicts with absolute/relative energies (eV), the
    geometry, and an optional error string.
    """
    eng = get_engine(engine)
    demo = not eng.available

    results: list[dict] = []
    ref_e = None  # reference at largest separation

    if demo:
        pes_rel = demo_pes(distances)

    for idx, d in enumerate(distances):
        combined = build_reaction_system(atoms_a, atoms_b, catalyst_atoms, d)

        if demo:
            e_rel = float(pes_rel[idx])
            e_abs = e_rel
        else:
            try:
                e_abs = eng.single_point(
                    combined, basis=basis, method=method, charge=catalyst_charge
                ).energy_eV
                if ref_e is None:
                    ref_e = e_abs
                e_rel = e_abs - ref_e
            except Exception as exc:
                results.append({
                    "distance": d, "energy_abs": None, "energy_rel": None,
                    "atoms": combined.copy(), "error": str(exc),
                })
                if progress_cb:
                    progress_cb(idx + 1, len(distances), d, None)
                continue

        if ref_e is None and not demo:
            ref_e = e_abs
        results.append({
            "distance": d, "energy_abs": e_abs, "energy_rel": e_rel,
            "atoms": combined.copy(), "error": None,
        })
        if progress_cb:
            progress_cb(idx + 1, len(distances), d, e_rel)

    return results


def decompose_energies(
    atoms_a, atoms_b, catalyst_atoms, catalyst_charge, engine, basis, method
) -> dict:
    """Single-point energies (eV) of each fragment separately."""
    eng = get_engine(engine)
    if not eng.available:
        rng = np.random.default_rng(0)

        def demo_e(a):
            return -sum(a.get_atomic_numbers()) * 0.5 + rng.normal(0, 0.05)

        eA = demo_e(atoms_a)
        eB = demo_e(atoms_b)
        eC = demo_e(catalyst_atoms) if catalyst_atoms else 0.0
        return {"E_A": eA, "E_B": eB, "E_cat": eC, "demo": True}

    def sp(a, charge):
        return eng.single_point(a, basis=basis, method=method, charge=charge).energy_eV

    eA = sp(atoms_a, 0)
    eB = sp(atoms_b, 0)
    eC = sp(catalyst_atoms, catalyst_charge) if catalyst_atoms else 0.0
    return {"E_A": eA, "E_B": eB, "E_cat": eC, "demo": False}
