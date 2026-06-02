import math
import unittest

from chemmodel.constants import ATM_TO_PSIA
from chemmodel.solubility import (
    SOLVENTS,
    MissingHenryParameterError,
    nitrogen_solubility,
    psig_to_psia,
    temperature_sweep,
)
from chemmodel.solubility.model import _frange_inclusive, solvent_parameters


class NitrogenSolubilityTests(unittest.TestCase):
    def test_psig_to_psia(self) -> None:
        self.assertAlmostEqual(psig_to_psia(50.0), 64.6959, places=4)

    def test_henry_calculation_matches_tabulated_reference_temperature(self) -> None:
        eo = SOLVENTS["EO"]
        result = nitrogen_solubility("EO", temperature_f=77.0, pressure_psig=50.0)
        expected_hx_psia = 2180.0 * ATM_TO_PSIA
        expected_hcc = expected_hx_psia / eo.solvent_molar_density_mol_per_ft3
        expected_concentration = psig_to_psia(50.0) / expected_hcc
        expected_x = expected_concentration / (
            eo.solvent_molar_density_mol_per_ft3 + expected_concentration
        )

        self.assertTrue(math.isclose(result.henry_cc_psia_ft3_per_mol, expected_hcc, rel_tol=0.03))
        self.assertTrue(math.isclose(result.dissolved_mol_per_ft3, expected_concentration, rel_tol=0.03))
        self.assertTrue(math.isclose(result.liquid_mole_fraction, expected_x, rel_tol=0.03))

    def test_temperature_sweep_row_count_and_monotonicity(self) -> None:
        rows = temperature_sweep("EO", pressure_psig=50.0, t_min_f=68.0, t_max_f=140.0, t_step_f=12.0)

        self.assertEqual(len(rows), 7)
        self.assertEqual(rows[0].temperature_f, 68.0)
        self.assertEqual(rows[-1].temperature_f, 140.0)

        concentrations = [row.dissolved_mol_per_ft3 for row in rows]
        self.assertEqual(concentrations, sorted(concentrations))

    def test_missing_po_parameters_fail_closed(self) -> None:
        with self.assertRaises(MissingHenryParameterError):
            nitrogen_solubility("PO", temperature_f=68.0, pressure_psig=50.0)

    # ── added edge cases ───────────────────────────────────────────────────────
    def test_unsupported_solvent_raises_keyerror(self) -> None:
        with self.assertRaises(KeyError):
            solvent_parameters("ZZ")

    def test_solvent_lookup_is_case_insensitive(self) -> None:
        self.assertIs(solvent_parameters("eo"), SOLVENTS["EO"])

    def test_frange_rejects_nonpositive_step(self) -> None:
        with self.assertRaises(ValueError):
            list(_frange_inclusive(0.0, 1.0, 0.0))

    def test_frange_rejects_reversed_bounds(self) -> None:
        with self.assertRaises(ValueError):
            list(_frange_inclusive(10.0, 1.0, 1.0))


if __name__ == "__main__":
    unittest.main()
