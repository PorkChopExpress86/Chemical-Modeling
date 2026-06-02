"""
Equilibrium solubility model for nitrogen in liquid ethylene oxide (EO)
and propylene oxide (PO).

The model is intentionally limited to dissolved nitrogen at equilibrium. It
does not calculate closed-vessel pressure rise or solvent vapor pressure.

Units used internally:
- Pressure: psia, unless an argument explicitly says psig
- Temperature: degrees Rankine
- Henry constant Hcc: psia*ft^3/mol, where c_N2 = P_N2 / Hcc
- Dissolved concentration: mol/ft^3 of liquid
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from chemmodel.constants import (
    ATM_TO_PSIA,
    G_PER_ML_TO_LB_PER_FT3,
    GAL_PER_FT3,
    LBMOL_TO_MOL,
    PSIA_PER_PSIG_OFFSET,
)


class MissingHenryParameterError(ValueError):
    """Raised when a solvent has no defensible nitrogen Henry correlation."""


@dataclass(frozen=True)
class HenryPoint:
    """One tabulated Henry point on a pressure/mole-fraction basis."""

    temperature_f: float
    henry_atm_per_mole_fraction: float


@dataclass(frozen=True)
class SolventParameters:
    """Physical and Henry-law parameters for one solvent."""

    key: str
    name: str
    molecular_weight_g_per_mol: float
    density_lb_per_ft3: float
    density_reference: str
    henry_points: Optional[Sequence[HenryPoint]]
    henry_reference: str
    notes: str = ""

    @property
    def solvent_molar_density_mol_per_ft3(self) -> float:
        """Return solvent molar density from liquid density and molecular weight."""
        return self.density_lb_per_ft3 / self.molecular_weight_g_per_mol * LBMOL_TO_MOL

    def _fit_ln_h_vs_inverse_t(self) -> tuple[float, float]:
        if not self.henry_points:
            raise MissingHenryParameterError(
                f"No direct nitrogen-in-{self.name} Henry-law parameters are available. "
                "Add literature or lab data before calculating solubility."
            )

        x_values = [1.0 / fahrenheit_to_rankine(point.temperature_f) for point in self.henry_points]
        y_values = [math.log(point.henry_atm_per_mole_fraction) for point in self.henry_points]
        x_mean = sum(x_values) / len(x_values)
        y_mean = sum(y_values) / len(y_values)
        numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
        denominator = sum((x - x_mean) ** 2 for x in x_values)
        if denominator == 0:
            raise ValueError(f"Henry data for {self.key} must span more than one temperature.")

        slope = numerator / denominator
        intercept = y_mean - slope * x_mean
        return intercept, slope

    def henry_x_atm_per_mole_fraction(self, temperature_f: float) -> float:
        """Return Hx where x_N2 = P_N2_atm / Hx."""
        intercept, slope = self._fit_ln_h_vs_inverse_t()
        temperature_r = fahrenheit_to_rankine(temperature_f)
        return math.exp(intercept + slope / temperature_r)

    def henry_cc_psia_ft3_per_mol(self, temperature_f: float) -> float:
        """Return Hcc where c_N2 = P_N2_psia / Hcc."""
        h_x_psia_per_mole_fraction = self.henry_x_atm_per_mole_fraction(temperature_f) * ATM_TO_PSIA
        return h_x_psia_per_mole_fraction / self.solvent_molar_density_mol_per_ft3


SOLVENTS: Dict[str, SolventParameters] = {
    "EO": SolventParameters(
        key="EO",
        name="ethylene oxide",
        molecular_weight_g_per_mol=44.0526,
        density_lb_per_ft3=0.8795 * G_PER_ML_TO_LB_PER_FT3,
        density_reference=(
            "IARC/NCBI Bookshelf, Some Industrial Chemicals: Ethylene Oxide; "
            "liquid density 0.8795 g/mL at 20 C."
        ),
        henry_points=(
            HenryPoint(32.0, 2800.0),
            HenryPoint(77.0, 2180.0),
            HenryPoint(122.0, 1820.0),
        ),
        henry_reference=(
            "LyondellBasell Ethylene Oxide Product Stewardship Guidance, Table A4; "
            "nitrogen Henry constants in liquid EO at 32, 77, and 122 F. "
            "Primary literature trail: Olson, J. Chem. Eng. Data 1977, 22, 326-329."
        ),
        notes="The 140 F sweep endpoint extrapolates beyond the highest tabulated EO point at 122 F.",
    ),
    "PO": SolventParameters(
        key="PO",
        name="propylene oxide",
        molecular_weight_g_per_mol=58.0791,
        density_lb_per_ft3=6.9 * GAL_PER_FT3,
        density_reference="PubChem/CAMEO Chemicals: propylene oxide density 6.9 lb/gal.",
        henry_points=None,
        henry_reference=(
            "No direct, authoritative nitrogen-in-liquid-propylene-oxide Henry data "
            "was found in the initial ACS/NIST/PubChem/handbook-oriented search."
        ),
        notes="Calculation intentionally fails closed until direct PO Henry data is supplied.",
    ),
}


@dataclass(frozen=True)
class SolubilityResult:
    """One equilibrium solubility result row."""

    solvent: str
    solvent_name: str
    temperature_f: float
    pressure_psig: float
    pressure_psia: float
    henry_cc_psia_ft3_per_mol: float
    dissolved_mol_per_ft3: float
    liquid_mole_fraction: float


def fahrenheit_to_rankine(temperature_f: float) -> float:
    return temperature_f + 459.67


def psig_to_psia(pressure_psig: float) -> float:
    return pressure_psig + PSIA_PER_PSIG_OFFSET


def _frange_inclusive(start: float, stop: float, step: float) -> Iterable[float]:
    if step <= 0:
        raise ValueError("Temperature step must be greater than zero.")
    if stop < start:
        raise ValueError("Temperature stop must be greater than or equal to start.")

    value = start
    epsilon = step * 1.0e-9
    while value <= stop + epsilon:
        yield round(value, 10)
        value += step


def solvent_parameters(solvent: str) -> SolventParameters:
    key = solvent.strip().upper()
    try:
        return SOLVENTS[key]
    except KeyError as exc:
        supported = ", ".join(sorted(SOLVENTS))
        raise KeyError(f"Unsupported solvent {solvent!r}. Supported solvents: {supported}") from exc


def nitrogen_solubility(
    solvent: str,
    temperature_f: float,
    pressure_psig: float = 50.0,
) -> SolubilityResult:
    """Calculate equilibrium nitrogen solubility for one solvent and temperature."""
    params = solvent_parameters(solvent)
    pressure_psia = psig_to_psia(pressure_psig)
    henry_cc = params.henry_cc_psia_ft3_per_mol(temperature_f)
    dissolved = pressure_psia / henry_cc
    solvent_molar_density = params.solvent_molar_density_mol_per_ft3
    mole_fraction = dissolved / (solvent_molar_density + dissolved)

    return SolubilityResult(
        solvent=params.key,
        solvent_name=params.name,
        temperature_f=temperature_f,
        pressure_psig=pressure_psig,
        pressure_psia=pressure_psia,
        henry_cc_psia_ft3_per_mol=henry_cc,
        dissolved_mol_per_ft3=dissolved,
        liquid_mole_fraction=mole_fraction,
    )


def temperature_sweep(
    solvent: str,
    pressure_psig: float = 50.0,
    t_min_f: float = 68.0,
    t_max_f: float = 140.0,
    t_step_f: float = 12.0,
) -> List[SolubilityResult]:
    """Return equilibrium nitrogen solubility rows across an inclusive F sweep."""
    return [
        nitrogen_solubility(solvent, temperature_f, pressure_psig)
        for temperature_f in _frange_inclusive(t_min_f, t_max_f, t_step_f)
    ]
