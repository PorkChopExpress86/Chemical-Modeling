"""PySCF backend — molecular and periodic (PBC) single points."""

from __future__ import annotations

from typing import Any

import numpy as np
from ase import Atoms

from chemmodel.chemistry.structures import atoms_to_pyscf
from chemmodel.constants import BOHR_ANG, HARTREE_EV
from chemmodel.engines.base import EngineResult, QuantumEngine


def _finite_diff_pbc(cell, kpts, basis):
    """Central-difference forces for a periodic cell (no analytic PBC grad)."""
    from pyscf.pbc import gto as pbcgto, scf as pbcscf  # type: ignore

    pos0 = np.array(cell.atom_coords())
    syms = [a[0] for a in cell.atom]
    forces = np.zeros_like(pos0)
    for i in range(len(pos0)):
        for j in range(3):
            results = []
            for sign in (+1, -1):
                p = pos0.copy()
                p[i, j] += sign * 0.01
                c = pbcgto.Cell()
                c.atom = list(zip(syms, p.tolist()))
                c.basis = basis
                c.a = cell.a
                c.build(verbose=0)
                m = pbcscf.RHF(c)
                m.kpts = c.make_kpts([kpts, kpts, kpts])
                results.append(m.kernel())
            forces[i, j] = -(results[0] - results[1]) / 0.02
    return forces * HARTREE_EV


class PySCFEngine(QuantumEngine):
    name = "pyscf"
    label = "PySCF — QC + Periodic BC"

    @property
    def available(self) -> bool:
        try:
            import pyscf  # type: ignore  # noqa: F401

            return True
        except ImportError:
            return False

    def single_point(
        self,
        atoms: Atoms,
        *,
        basis: str,
        method: str,
        charge: int = 0,
        with_forces: bool = False,
        use_pbc: bool = False,
        cell_size: float = 15.0,
        kpts: int = 1,
        **options: Any,
    ) -> EngineResult:
        if not with_forces:
            return self._energy_only(atoms, basis, method, charge)
        if use_pbc:
            return self._pbc_with_forces(atoms, basis, cell_size, kpts)
        return self._molecular_with_forces(atoms, basis, method)

    # ── energy-only (charge/spin aware) — used by PES scans ────────────────────
    def _energy_only(self, atoms, basis, method, charge) -> EngineResult:
        from pyscf import dft, gto, scf  # type: ignore

        n_elec = sum(atoms.get_atomic_numbers()) - charge
        spin = n_elec % 2  # 0 = singlet, 1 = doublet

        mol = gto.Mole(verbose=0)
        mol.atom = atoms_to_pyscf(atoms)
        mol.basis = basis
        mol.charge = charge
        mol.spin = spin
        mol.build()

        m_upper = method.upper()
        if m_upper in ("HF", "RHF") and spin == 0:
            mf = scf.RHF(mol)
        elif m_upper == "UHF" or spin != 0:
            mf = scf.UHF(mol)
        else:
            mf = dft.RKS(mol)
            mf.xc = method

        mf.max_cycle = 200
        e = mf.kernel()
        if not mf.converged:
            raise RuntimeError(f"SCF did not converge (charge={charge}, spin={spin})")
        return EngineResult(float(e) * HARTREE_EV, None, f"PySCF {method}/{basis}")

    # ── molecular single point with analytic forces ───────────────────────────
    def _molecular_with_forces(self, atoms, basis, method) -> EngineResult:
        from pyscf import dft, grad, gto, scf  # type: ignore

        mol = gto.Mole(verbose=0)
        mol.atom = atoms_to_pyscf(atoms)
        mol.basis = basis
        mol.build()
        m_upper = method.upper()
        if m_upper in ("HF", "RHF"):
            mf = scf.RHF(mol)
        elif m_upper == "UHF":
            mf = scf.UHF(mol)
        else:
            mf = dft.RKS(mol)
            mf.xc = method
        e_tot = mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")
        g_mat = grad.RHF(mf).kernel()
        forces = -g_mat * HARTREE_EV / BOHR_ANG
        return EngineResult(e_tot * HARTREE_EV, forces, f"PySCF {method}/{basis}")

    # ── periodic single point with finite-difference forces ────────────────────
    def _pbc_with_forces(self, atoms, basis, cell_size, kpts) -> EngineResult:
        from pyscf.pbc import gto as pbcgto, scf as pbcscf  # type: ignore

        cell = pbcgto.Cell()
        cell.atom = atoms_to_pyscf(atoms)
        cell.basis = basis
        cell.a = np.eye(3) * cell_size
        cell.build()
        mf = pbcscf.RHF(cell)
        mf.kpts = cell.make_kpts([kpts, kpts, kpts])
        e_tot = mf.kernel()
        forces = _finite_diff_pbc(cell, kpts, basis)
        note = f"PBC k-grid {kpts}³, cell {cell_size} Å"
        return EngineResult(e_tot * HARTREE_EV, forces, note)
