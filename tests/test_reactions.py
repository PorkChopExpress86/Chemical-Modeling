import unittest

import numpy as np

from chemmodel.chemistry import xyz_to_atoms
from chemmodel.engines import get_engine
from chemmodel.reactions import decompose_energies, demo_pes, run_reaction_scan

A = xyz_to_atoms("2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74")
B = xyz_to_atoms("1\nH\nH 0.0 0.0 0.0")


class DemoPesTests(unittest.TestCase):
    def test_shape_and_well(self) -> None:
        d = np.linspace(4.5, 0.9, 12)
        pes = demo_pes(d)
        self.assertEqual(pes.shape, d.shape)
        self.assertLess(pes.min(), 0.0)        # Morse well dips below zero


class DemoScanTests(unittest.TestCase):
    """Force the demo path by selecting an engine that is not installed (cp2k)."""

    def setUp(self) -> None:
        if get_engine("cp2k").available:
            self.skipTest("cp2k available — cannot exercise demo fallback")

    def test_scan_returns_points_without_errors(self) -> None:
        dists = np.linspace(4.5, 0.9, 6)
        scan = run_reaction_scan(A, B, None, 0, dists, "cp2k", "sto-3g", "HF")
        self.assertEqual(len(scan), 6)
        self.assertTrue(all(p["error"] is None for p in scan))
        self.assertTrue(all(np.isfinite(p["energy_rel"]) for p in scan))

    def test_decompose_demo_keys(self) -> None:
        dec = decompose_energies(A, B, None, 0, "cp2k", "sto-3g", "HF")
        self.assertEqual(set(dec), {"E_A", "E_B", "E_cat", "demo"})
        self.assertTrue(dec["demo"])


class RealScanTests(unittest.TestCase):
    def setUp(self) -> None:
        if not get_engine("pyscf").available:
            self.skipTest("pyscf not installed")

    def test_two_point_real_scan(self) -> None:
        dists = np.array([4.0, 2.0])
        scan = run_reaction_scan(A, B, None, 0, dists, "pyscf", "sto-3g", "HF")
        self.assertEqual(len(scan), 2)
        self.assertTrue(all(p["error"] is None for p in scan))
        # first point is the reference → relative energy exactly 0
        self.assertAlmostEqual(scan[0]["energy_rel"], 0.0)


if __name__ == "__main__":
    unittest.main()
