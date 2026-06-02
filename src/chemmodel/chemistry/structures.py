"""Molecular structure helpers built on ASE ``Atoms`` (the single source of truth
for geometry). Pure geometry/assembly — no quantum chemistry, no Streamlit.
"""

from __future__ import annotations

import io

import numpy as np
from ase import Atoms
from ase.io import read as ase_read

# Placement grid (fractional box coords), one slot per molecule lane.
PLACEMENT = [
    (0.30, 0.30, 0.50),
    (0.70, 0.30, 0.50),
    (0.30, 0.70, 0.50),
    (0.70, 0.70, 0.50),
]


def xyz_to_atoms(xyz_text: str) -> Atoms:
    """Parse an XYZ string into an ASE ``Atoms`` object."""
    return ase_read(io.StringIO(xyz_text.strip()), format="xyz")


def atoms_to_pyscf(atoms: Atoms) -> list:
    """Convert ASE ``Atoms`` to PySCF's ``[(symbol, (x, y, z)), ...]`` format."""
    return [
        (s, tuple(p))
        for s, p in zip(atoms.get_chemical_symbols(), atoms.get_positions())
    ]


def center(atoms: Atoms) -> Atoms:
    """Return a copy translated so its centroid sits at the origin."""
    a = atoms.copy()
    a.positions -= a.positions.mean(axis=0)
    return a


def build_reaction_system(
    atoms_a: Atoms,
    atoms_b: Atoms,
    catalyst: Atoms | None,
    separation: float,
    cat_lateral: float = 4.0,
) -> Atoms:
    """Lay out a reaction system along the z-axis reaction coordinate.

    A is centred at the origin; B is centred at
    ``(0, 0, separation + half_extent_a + half_extent_b)``; the catalyst (if any)
    is centred at ``(cat_lateral, 0, 0)``.
    """
    a = center(atoms_a.copy())
    b = center(atoms_b.copy())

    ext_a = (a.positions[:, 2].max() - a.positions[:, 2].min()) / 2
    ext_b = (b.positions[:, 2].max() - b.positions[:, 2].min()) / 2
    b.positions += np.array([0.0, 0.0, ext_a + separation + ext_b])

    parts = [a, b]
    if catalyst is not None:
        c = center(catalyst.copy())
        c.positions += np.array([cat_lateral, 0.0, 0.0])
        parts.append(c)

    combined = parts[0]
    for p in parts[1:]:
        combined = combined + p
    return combined


def assemble_system(mol_atoms: list, lane_names: list[str], box_L: float):
    """Center each molecule and drop it on a placement-grid slot inside the box.

    Returns ``(combined_atoms, group_sizes, min_pair_distance)``.
    """
    combined = Atoms(cell=[box_L, box_L, box_L], pbc=False)
    group_sizes = []
    centers = []
    for i, m in enumerate(mol_atoms):
        a = m.copy()
        a.positions -= a.positions.mean(axis=0)  # center
        fx, fy, fz = PLACEMENT[i]
        shift = np.array([fx * box_L, fy * box_L, fz * box_L])
        a.positions += shift
        centers.append(shift)
        combined += a
        group_sizes.append(len(a))

    # smallest center-to-center distance (overlap heuristic)
    min_d = np.inf
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            min_d = min(min_d, np.linalg.norm(centers[i] - centers[j]))
    if len(centers) < 2:
        min_d = box_L
    return combined, group_sizes, min_d


def auto_multiplicity(atoms: Atoms) -> tuple[int, int]:
    """Neutral system; multiplicity from electron-count parity. Returns (charge, mult)."""
    n_elec = int(sum(atoms.get_atomic_numbers()))
    charge = 0
    mult = 1 if (n_elec - charge) % 2 == 0 else 2
    return charge, mult
