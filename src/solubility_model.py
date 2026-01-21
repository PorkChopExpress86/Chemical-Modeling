"""
Solubility and pressure-change model for nitrogen in liquid EO/PO using Henry's law.

This module keeps units explicit:
- Pressure: psia (pounds per square inch absolute)
- Temperature: °R (degrees Rankine)
- Volume: ft^3 (cubic feet)
- Henry's constant Hcp: psia·ft^3/mol (the proportionality P = Hcp * x, where x is mole fraction of gas in liquid)

Key assumptions
---------------
1) Nitrogen behaves ideally in the gas phase.
2) Henry's law is valid in the dilute regime for nitrogen in EO/PO.
3) Temperature-dependence of Henry's constant follows van't Hoff form:
      H(T) = H_ref * exp( (ΔH_sol/R) * (1/T_ref - 1/T) )
   where ΔH_sol is the enthalpy of solution (J/mol). A positive ΔH_sol typically means
   solubility decreases with increasing temperature.
4) Closed system: total moles of nitrogen are conserved when heating from storage to reaction temperature.
5) Only nitrogen is considered in the gas phase (other vapors ignored for simplicity).

You must supply solvent-specific Henry parameters. Default example values are illustrative only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, List
import math
import csv
from pathlib import Path

R_GAS = 10.7316  # psia·ft^3/(mol·°R)
BAR_TO_PSIA = 14.5038
PSIA_TO_BAR = 1.0 / 14.5038


def henry_constant(T: float, H_ref: float, T_ref: float, delta_h_sol: float) -> float:
    """Return Henry's constant at temperature T (°R).

    Parameters
    ----------
    T : float
        Temperature of interest in °R (Rankine).
    H_ref : float
        Henry's constant at reference temperature (psia·ft^3/mol).
    T_ref : float
        Reference temperature in °R (Rankine).
    delta_h_sol : float
        Enthalpy of solution (BTU/mol). Positive values usually mean decreasing solubility with higher T.
    """
    exponent = (delta_h_sol / R_GAS) * (1.0 / T_ref - 1.0 / T)
    return H_ref * math.exp(exponent)


def dissolved_moles(p_partial: float, v_liquid: float, henry_T: float) -> float:
    """Moles of gas dissolved at partial pressure p_partial (psia) in volume v_liquid (ft^3).

    Uses Henry's law: p = H * x  => n = (p/H) * n_liquid; for dilute gas we can simplify to
    n_gas = p * V_liquid / H, treating the gas mole fraction times total moles as p/H.
    This approximation is standard when the gas solubility is very low and the solvent is in large excess.
    """
    return p_partial * v_liquid / henry_T


def total_moles_closed_system(
    p_initial: float,
    v_headspace: float,
    v_liquid: float,
    t_initial: float,
    h_ref: float,
    t_ref: float,
    delta_h_sol: float,
) -> float:
    """Total moles of nitrogen in a closed tank at initial conditions.

    Includes dissolved plus headspace gas.
    """
    h_t0 = henry_constant(t_initial, h_ref, t_ref, delta_h_sol)
    n_gas0 = p_initial * v_headspace / (R_GAS * t_initial)
    n_dissolved0 = dissolved_moles(p_initial, v_liquid, h_t0)
    return n_gas0 + n_dissolved0


def final_pressure_closed(
    n_total: float,
    v_headspace: float,
    v_liquid: float,
    t_final: float,
    h_ref: float,
    t_ref: float,
    delta_h_sol: float,
) -> float:
    """Solve for final nitrogen partial pressure after heating a closed system.

    Mass balance: n_total = n_gas + n_dissolved = (P*V_g)/(R*T) + (P*V_l)/H(T)
    => P = n_total / ( V_g/(R*T) + V_l/H(T) )
    """
    h_tf = henry_constant(t_final, h_ref, t_ref, delta_h_sol)
    denominator = (v_headspace / (R_GAS * t_final)) + (v_liquid / h_tf)
    return n_total / denominator


@dataclass
class NitrogenSolubilityCase:
    """Configuration for a single tank scenario."""

    name: str
    volume_liquid_ft3: float
    volume_headspace_ft3: float
    p_initial_psia: float
    t_initial_f: float
    t_final_f: float
    henry_ref_psia_ft3_per_mol: float
    henry_tref_f: float
    delta_h_sol_btu_per_mol: float
    notes: Optional[str] = None

    def run(self) -> Dict[str, float]:
        """Compute dissolved moles at storage, total moles, and final pressure after heating.

        Returns values in US customary units, with pressure in psia.
        """
        t0_R = self.t_initial_f + 459.67
        tf_R = self.t_final_f + 459.67
        p0_psia = self.p_initial_psia

        n_total = total_moles_closed_system(
            p_initial=p0_psia,
            v_headspace=self.volume_headspace_ft3,
            v_liquid=self.volume_liquid_ft3,
            t_initial=t0_R,
            h_ref=self.henry_ref_psia_ft3_per_mol,
            t_ref=self.henry_tref_f + 459.67,
            delta_h_sol=self.delta_h_sol_btu_per_mol,
        )

        p_final_psia = final_pressure_closed(
            n_total=n_total,
            v_headspace=self.volume_headspace_ft3,
            v_liquid=self.volume_liquid_ft3,
            t_final=tf_R,
            h_ref=self.henry_ref_psia_ft3_per_mol,
            t_ref=self.henry_tref_f + 459.67,
            delta_h_sol=self.delta_h_sol_btu_per_mol,
        )

        h_t0 = henry_constant(t0_R, self.henry_ref_psia_ft3_per_mol, self.henry_tref_f + 459.67, self.delta_h_sol_btu_per_mol)
        n_dissolved0 = dissolved_moles(p0_psia, self.volume_liquid_ft3, h_t0)

        return {
            "initial_pressure_psia": p0_psia,
            "initial_henry_psia_ft3_per_mol": h_t0,
            "initial_dissolved_mol": n_dissolved0,
            "total_n_mol": n_total,
            "final_pressure_psia": p_final_psia,
            "final_pressure_bar": p_final_psia * PSIA_TO_BAR,
            "delta_pressure_psi": p_final_psia - p0_psia,
        }


# Example reference parameters (illustrative — replace with plant data!)
DEFAULT_PARAMS = {
    "N2_in_EO": {
        "H_ref": 1.16e6,  # psia·ft^3/mol at 77 °F (illustrative)
        "T_ref_F": 77.0,
        "delta_h_sol": 10.4,  # BTU/mol, typical for N2 in many organics (approx.)
    },
    "N2_in_PO": {
        "H_ref": 1.01e6,  # psia·ft^3/mol at 77 °F (illustrative)
        "T_ref_F": 77.0,
        "delta_h_sol": 10.4,
    },
}


def example_case() -> Dict[str, Dict[str, float]]:
    """Run simple EO and PO examples with placeholder parameters.

    The numbers here are *not* site-specific; update with lab/handbook data before using for decisions.
    """
    results = {}

    eo = NitrogenSolubilityCase(
        name="EO holding tank",
        volume_liquid_ft3=353.0,
        volume_headspace_ft3=70.6,
        p_initial_psia=87.0,
        t_initial_f=68.0,
        t_final_f=140.0,
        henry_ref_psia_ft3_per_mol=DEFAULT_PARAMS["N2_in_EO"]["H_ref"],
        henry_tref_f=DEFAULT_PARAMS["N2_in_EO"]["T_ref_F"],
        delta_h_sol_btu_per_mol=DEFAULT_PARAMS["N2_in_EO"]["delta_h_sol"],
    )
    results[eo.name] = eo.run()

    po = NitrogenSolubilityCase(
        name="PO holding tank",
        volume_liquid_ft3=353.0,
        volume_headspace_ft3=70.6,
        p_initial_psia=87.0,
        t_initial_f=68.0,
        t_final_f=140.0,
        henry_ref_psia_ft3_per_mol=DEFAULT_PARAMS["N2_in_PO"]["H_ref"],
        henry_tref_f=DEFAULT_PARAMS["N2_in_PO"]["T_ref_F"],
        delta_h_sol_btu_per_mol=DEFAULT_PARAMS["N2_in_PO"]["delta_h_sol"],
    )
    results[po.name] = po.run()

    return results


def run_reactors_from_csv(csv_path: str = "reactor_specs_template.csv") -> Dict[str, Dict[str, float]]:
    """Load reactor specifications from CSV and run nitrogen solubility calculations for each.

    CSV format:
        Reactor Name,Liquid Volume (ft3),Headspace Volume (ft3),Initial Pressure (psia),Storage Temperature (F),Reaction Temperature (F),Solvent

    Parameters
    ----------
    csv_path : str
        Path to the CSV file with reactor specifications.

    Returns
    -------
    Dict[str, Dict[str, float]]
        Results for each reactor keyed by reactor name.
    """
    results = {}
    csv_file = Path(csv_path)
    
    if not csv_file.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            reactor_name = row["Reactor Name"]
            solvent = row["Solvent"].strip().upper()
            
            if solvent not in ["EO", "PO"]:
                print(f"Warning: Unknown solvent '{solvent}' for {reactor_name}, skipping.")
                continue
            
            # Select appropriate Henry parameters
            params = DEFAULT_PARAMS[f"N2_in_{solvent}"]
            
            case = NitrogenSolubilityCase(
                name=reactor_name,
                volume_liquid_ft3=float(row["Liquid Volume (ft3)"]),
                volume_headspace_ft3=float(row["Headspace Volume (ft3)"]),
                p_initial_psia=float(row["Initial Pressure (psia)"]),
                t_initial_f=float(row["Storage Temperature (F)"]),
                t_final_f=float(row["Reaction Temperature (F)"]),
                henry_ref_psia_ft3_per_mol=params["H_ref"],
                henry_tref_f=params["T_ref_F"],
                delta_h_sol_btu_per_mol=params["delta_h_sol"],
            )
            results[reactor_name] = case.run()
    
    return results


def print_reactor_comparison(results: Dict[str, Dict[str, float]]) -> None:
    """Print a formatted comparison table of reactor pressure calculations.

    Parameters
    ----------
    results : Dict[str, Dict[str, float]]
        Results dictionary from run_reactors_from_csv or example_case.
    """
    print("\n" + "=" * 90)
    print("REACTOR NITROGEN PRESSURE ANALYSIS")
    print("=" * 90)
    print(f"{'Reactor':<12} {'Initial (psia)':<15} {'Final (psia)':<15} {'ΔP (psi)':<12} {'Final (bar)':<12}")
    print("-" * 90)
    
    for reactor_name, data in results.items():
        initial_psia = data["initial_pressure_psia"]
        final_psia = data["final_pressure_psia"]
        delta_p = data["delta_pressure_psi"]
        final_bar = data["final_pressure_bar"]
        
        print(f"{reactor_name:<12} {initial_psia:<15.2f} {final_psia:<15.2f} {delta_p:<12.2f} {final_bar:<12.2f}")
    
    print("=" * 90)
    print("\nNOTE: Using placeholder Henry's law parameters. Replace with measured data before decisions.\n")


if __name__ == "__main__":
    from pprint import pprint

    pprint(example_case())
    
    # Uncomment to run reactors from CSV:
    # results = run_reactors_from_csv("reactor_specs_template.csv")
    # print_reactor_comparison(results)
