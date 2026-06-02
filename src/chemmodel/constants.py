"""Physical constants and unit conversions — the single source of truth.

Every module imports the values it needs from here so that the quantum-chemistry
apps, the reactor kinetics, and the solubility model can never drift apart on a
conversion factor. Energies inside the quantum-chemistry layer are handled in
electron-volts (eV); engines convert from Hartree at their boundary.
"""

from __future__ import annotations

# ── quantum chemistry ─────────────────────────────────────────────────────────
HARTREE_EV = 27.211386                      # eV per Hartree
BOHR_ANG = 0.529177                         # Å per Bohr
HA_BOHR_TO_EV_ANG = HARTREE_EV / BOHR_ANG   # eV/Å per Hartree/Bohr (force unit)

# ── reactor / thermodynamics ──────────────────────────────────────────────────
R_GAS = 8.314462618          # J/mol/K — universal gas constant
EV_TO_KJMOL = 96.485332      # kJ/mol per eV
L_TO_FT3 = 0.0353146667      # ft³ per litre

# ── solubility model (US customary) ───────────────────────────────────────────
ATM_TO_PSIA = 14.6959
PSIA_PER_PSIG_OFFSET = 14.6959
LBMOL_TO_MOL = 453.59237
G_PER_ML_TO_LB_PER_FT3 = 62.42796
GAL_PER_FT3 = 7.48051948
