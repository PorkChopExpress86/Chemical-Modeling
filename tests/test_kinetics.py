import math
import unittest

import numpy as np

from chemmodel.reactors import (
    arrhenius_k,
    batch_profile,
    batch_time_for_X,
    inv_rate,
    k_units,
    tau_cstr,
    tau_pfr,
)


class BatchProfileTests(unittest.TestCase):
    def test_zero_order_linear_decay_and_clip(self) -> None:
        CA, X = batch_profile(1.0, 0.1, 0, np.array([0.0, 2.0, 20.0]))
        self.assertAlmostEqual(CA[0], 1.0)
        self.assertAlmostEqual(CA[1], 0.8)
        self.assertEqual(CA[2], 0.0)          # clipped, never negative
        self.assertAlmostEqual(X[1], 0.2)
        self.assertLessEqual(X.max(), 1.0)

    def test_first_order_half_life(self) -> None:
        k = 0.05
        t_half = math.log(2.0) / k
        CA, X = batch_profile(1.0, k, 1, np.array([t_half]))
        self.assertAlmostEqual(float(CA[0]), 0.5, places=6)
        self.assertAlmostEqual(float(X[0]), 0.5, places=6)

    def test_second_order_closed_form(self) -> None:
        CA, _ = batch_profile(2.0, 0.1, 2, np.array([5.0]))
        self.assertAlmostEqual(float(CA[0]), 2.0 / (1.0 + 2.0 * 0.1 * 5.0), places=9)


class BatchTimeTests(unittest.TestCase):
    def test_time_for_X_inverts_profile(self) -> None:
        for order in (0, 1, 2):
            t = batch_time_for_X(1.0, 0.05, order, 0.6)
            _, X = batch_profile(1.0, 0.05, order, np.array([t]))
            self.assertAlmostEqual(float(X[0]), 0.6, places=5, msg=f"order {order}")


class ReactorSizingTests(unittest.TestCase):
    def test_cstr_needs_more_volume_than_pfr(self) -> None:
        for order in (1, 2):
            tp = float(tau_pfr(1.0, 0.02, order, 0.8))
            tc = float(tau_cstr(1.0, 0.02, order, 0.8))
            self.assertGreaterEqual(tc, tp, msg=f"order {order}")

    def test_pfr_and_cstr_equal_for_zero_order(self) -> None:
        tp = float(tau_pfr(1.0, 0.1, 0, 0.5))
        tc = float(tau_cstr(1.0, 0.1, 0, 0.5))
        self.assertAlmostEqual(tp, tc, places=9)

    def test_inv_rate_positive(self) -> None:
        self.assertTrue(np.all(inv_rate(1.0, 0.05, 1, np.linspace(0.0, 0.9, 10)) > 0))


class ArrheniusTests(unittest.TestCase):
    def test_k_increases_with_temperature(self) -> None:
        cold = arrhenius_k(1e13, 50.0, 300.0)
        hot = arrhenius_k(1e13, 50.0, 500.0)
        self.assertGreater(hot, cold)

    def test_units_table(self) -> None:
        self.assertEqual(k_units(0), "mol·L⁻¹·s⁻¹")
        self.assertEqual(k_units(1), "s⁻¹")
        self.assertEqual(k_units(2), "L·mol⁻¹·s⁻¹")


if __name__ == "__main__":
    unittest.main()
