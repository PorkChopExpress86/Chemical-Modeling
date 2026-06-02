"""Closed-form ideal-reactor kinetics for a single nth-order reaction.

Rate law ``−rA = k · CA^n`` in an ideal Batch / CSTR / PFR. Pure NumPy — runs
instantly, no quantum chemistry. Concentrations in mol/L, time/space-time in s.
"""

from __future__ import annotations

import numpy as np

from chemmodel.constants import R_GAS


def k_units(order: int) -> str:
    """Units of the rate constant for a given reaction order."""
    return {0: "mol·L⁻¹·s⁻¹", 1: "s⁻¹", 2: "L·mol⁻¹·s⁻¹"}.get(order, "(mixed)")


def arrhenius_k(A: float, Ea_kJmol: float, T_K: float) -> float:
    """Arrhenius rate constant ``k = A · exp(−Ea / R T)``."""
    return A * np.exp(-Ea_kJmol * 1000.0 / (R_GAS * T_K))


def batch_profile(CA0, k, order, t):
    """Concentration and conversion vs time for a constant-volume batch reactor."""
    t = np.asarray(t, dtype=float)
    if order == 0:
        CA = np.clip(CA0 - k * t, 0.0, None)
    elif order == 1:
        CA = CA0 * np.exp(-k * t)
    else:  # order == 2
        CA = CA0 / (1.0 + CA0 * k * t)
    X = np.clip(1.0 - CA / CA0, 0.0, 1.0)
    return CA, X


def batch_time_for_X(CA0, k, order, X):
    """Batch time required to reach conversion X."""
    X = min(X, 0.999999)
    if order == 0:
        return CA0 * X / k
    if order == 1:
        return -np.log(1.0 - X) / k
    return X / (k * CA0 * (1.0 - X))  # order 2


def tau_pfr(CA0, k, order, X):
    """Plug-flow space-time to reach conversion X (array-safe)."""
    X = np.asarray(X, dtype=float)
    if order == 0:
        return CA0 * X / k
    if order == 1:
        return -np.log(1.0 - X) / k
    return X / (k * CA0 * (1.0 - X))  # order 2


def tau_cstr(CA0, k, order, X):
    """CSTR space-time to reach conversion X (array-safe)."""
    X = np.asarray(X, dtype=float)
    CA = CA0 * (1.0 - X)
    rate = k * CA ** order
    return CA0 * X / rate


def inv_rate(CA0, k, order, X):
    """1 / (−rA) as a function of conversion — the Levenspiel integrand."""
    X = np.asarray(X, dtype=float)
    CA = CA0 * (1.0 - X)
    return 1.0 / (k * CA ** order)
