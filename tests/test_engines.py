import unittest

import numpy as np

from chemmodel.chemistry import xyz_to_atoms
from chemmodel.engines import available_engines, get_engine, probe

H2 = xyz_to_atoms("2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74")


class RegistryTests(unittest.TestCase):
    def test_available_engines_shape(self) -> None:
        avail = available_engines()
        self.assertEqual(set(avail), {"psi4", "pyscf", "cp2k", "demo"})
        self.assertTrue(avail["demo"])                      # demo is always available
        self.assertTrue(all(isinstance(v, bool) for v in avail.values()))

    def test_get_engine_is_case_insensitive(self) -> None:
        self.assertIs(get_engine("PySCF"), get_engine("pyscf"))

    def test_get_unknown_engine_raises(self) -> None:
        with self.assertRaises(KeyError):
            get_engine("nope")

    def test_probe(self) -> None:
        self.assertTrue(probe("numpy"))
        self.assertFalse(probe("a_module_that_does_not_exist_zzz"))


class DemoEngineTests(unittest.TestCase):
    def test_always_available(self) -> None:
        self.assertTrue(get_engine("demo").available)

    def test_energy_only(self) -> None:
        r = get_engine("demo").single_point(H2, basis="x", method="y")
        self.assertTrue(np.isfinite(r.energy_eV))
        self.assertIsNone(r.forces)

    def test_with_forces_shape(self) -> None:
        r = get_engine("demo").single_point(H2, with_forces=True)
        self.assertEqual(r.forces.shape, (len(H2), 3))


class PySCFEngineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.eng = get_engine("pyscf")
        if not self.eng.available:
            self.skipTest("pyscf not installed")

    def test_h2_energy_only_is_bound(self) -> None:
        r = self.eng.single_point(H2, basis="sto-3g", method="HF")
        self.assertLess(r.energy_eV, 0.0)
        self.assertIsNone(r.forces)

    def test_h2_with_forces(self) -> None:
        r = self.eng.single_point(H2, basis="sto-3g", method="HF", with_forces=True)
        self.assertLess(r.energy_eV, 0.0)
        self.assertEqual(r.forces.shape, (2, 3))


class Psi4EngineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.eng = get_engine("psi4")
        if not self.eng.available:
            self.skipTest("psi4 not installed")

    def test_h2_energy_only_is_bound(self) -> None:
        r = self.eng.single_point(H2, basis="sto-3g", method="hf")
        self.assertLess(r.energy_eV, 0.0)


if __name__ == "__main__":
    unittest.main()
