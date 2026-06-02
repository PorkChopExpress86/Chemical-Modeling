import unittest

import numpy as np

from chemmodel.chemistry import (
    assemble_system,
    atoms_to_pyscf,
    auto_multiplicity,
    build_reaction_system,
    center,
    xyz_to_atoms,
)

H2 = "2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74"
H_ATOM = "1\nH\nH 0.0 0.0 0.0"
WATER = "3\nWater\nO 0 0 0.12\nH 0 0.76 -0.48\nH 0 -0.76 -0.48"


class StructureTests(unittest.TestCase):
    def test_xyz_to_atoms_roundtrip(self) -> None:
        atoms = xyz_to_atoms(H2)
        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms.get_chemical_symbols(), ["H", "H"])

    def test_center_moves_centroid_to_origin(self) -> None:
        centered = center(xyz_to_atoms(WATER))
        self.assertTrue(np.allclose(centered.get_positions().mean(axis=0), 0.0, atol=1e-9))

    def test_atoms_to_pyscf_format(self) -> None:
        spec = atoms_to_pyscf(xyz_to_atoms(H2))
        self.assertEqual(len(spec), 2)
        sym, xyz = spec[0]
        self.assertEqual(sym, "H")
        self.assertEqual(len(xyz), 3)

    def test_build_reaction_system_places_b_beyond_a(self) -> None:
        a = xyz_to_atoms(H2)
        b = xyz_to_atoms(H_ATOM)
        sep = 3.0
        sys_ = build_reaction_system(a, b, None, sep)
        self.assertEqual(len(sys_), len(a) + len(b))
        # B is shifted along +z by at least the requested separation
        self.assertGreater(sys_.get_positions()[:, 2].max(), sep)

    def test_build_reaction_system_includes_catalyst(self) -> None:
        a = xyz_to_atoms(H2)
        b = xyz_to_atoms(H2)
        cat = xyz_to_atoms(H_ATOM)
        sys_ = build_reaction_system(a, b, cat, 3.0)
        self.assertEqual(len(sys_), len(a) + len(b) + len(cat))

    def test_assemble_system_group_sizes_and_distance(self) -> None:
        mols = [xyz_to_atoms(H2), xyz_to_atoms(WATER)]
        combined, group_sizes, min_d = assemble_system(mols, ["Reactant A", "Reactant B"], 16.0)
        self.assertEqual(group_sizes, [2, 3])
        self.assertEqual(len(combined), 5)
        self.assertTrue(np.isfinite(min_d))
        self.assertGreater(min_d, 0.0)

    def test_assemble_single_molecule_distance_is_box(self) -> None:
        _, _, min_d = assemble_system([xyz_to_atoms(H2)], ["Reactant A"], 12.0)
        self.assertEqual(min_d, 12.0)

    def test_auto_multiplicity_parity(self) -> None:
        # H2: 2 electrons → singlet (mult 1)
        charge, mult = auto_multiplicity(xyz_to_atoms(H2))
        self.assertEqual((charge, mult), (0, 1))
        # single H: 1 electron → doublet (mult 2)
        _, mult_h = auto_multiplicity(xyz_to_atoms(H_ATOM))
        self.assertEqual(mult_h, 2)


if __name__ == "__main__":
    unittest.main()
