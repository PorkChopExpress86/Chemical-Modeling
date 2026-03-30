"""
Batch sulfonation reaction model for olefin + sodium bisulfite with peroxide initiation.

This module models a semi-batch radical addition of sodium bisulfite (NaHSO₃)
to a terminal olefin (e.g., 1-octene), catalyzed by an organic peroxide initiator
(e.g., tert-butyl peroxyacetate).

Reaction (1:1 molar stoichiometry)::

    R-CH=CH₂  +  NaHSO₃  →  R-CH₂-CH₂-SO₃Na

Process overview
----------------
* The olefin is charged to the reactor initially.
* Aqueous sodium bisulfite solution is fed continuously at a constant mass rate.
* The organic peroxide catalyst is added in equal discrete shots at regular
  bisulfite-solution increments.
* The reactor is maintained at a constant temperature (e.g., 190 °F).

This module keeps units explicit:
- Mass: lb (pounds)
- Temperature: °F (degrees Fahrenheit)
- Time: min (minutes)
- Molecular weight: g/mol
- Moles: lb-mol (pound-moles)

Key assumptions
---------------
1) The radical-initiated addition follows 1:1 olefin-to-bisulfite stoichiometry.
2) Bisulfite solution is fed at a constant mass rate.
3) Peroxide catalyst is added in equal-mass shots at uniform bisulfite-solution intervals.
4) Constant reaction temperature (isothermal).
5) Kinetic parameters are placeholders — replace with measured data before making
   process decisions.

WARNING: Default parameters and rate constants are illustrative only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# Unit conversion constants
# ---------------------------------------------------------------------------
LB_TO_G = 453.592        # grams per pound
G_TO_LB = 1.0 / LB_TO_G

# ---------------------------------------------------------------------------
# Molecular weights (g/mol)
# ---------------------------------------------------------------------------
MW_1_OCTENE = 112.21           # C₈H₁₆
MW_SODIUM_BISULFITE = 104.06   # NaHSO₃
MW_SODIUM_1_OCTANESULFONATE = 216.28  # C₈H₁₇SO₃Na
MW_TERT_BUTYL_PEROXYACETATE = 132.16  # C₆H₁₂O₃
MW_WATER = 18.015              # H₂O


def lb_to_lbmol(mass_lb: float, mw: float) -> float:
    """Convert pounds to pound-moles.

    Parameters
    ----------
    mass_lb : float
        Mass in pounds.
    mw : float
        Molecular weight in g/mol (numerically equal to lb/lb-mol).
    """
    return mass_lb / mw


def lbmol_to_lb(moles_lbmol: float, mw: float) -> float:
    """Convert pound-moles to pounds.

    Parameters
    ----------
    moles_lbmol : float
        Moles in pound-moles.
    mw : float
        Molecular weight in g/mol (numerically equal to lb/lb-mol).
    """
    return moles_lbmol * mw


@dataclass
class FeedScheduleStep:
    """One interval of the semi-batch feed schedule."""

    step: int
    time_start_min: float
    time_end_min: float
    bisulfite_solution_added_lb: float
    bisulfite_added_lb: float
    water_added_lb: float
    peroxide_shot_lb: float
    cumulative_bisulfite_solution_lb: float
    cumulative_peroxide_lb: float


@dataclass
class SulfonationReactionCase:
    """Configuration and solver for a batch olefin sulfonation reaction.

    The reaction:  Olefin + NaHSO₃ → Alkyl sulfonate (1:1 molar)
    initiated by radical peroxide catalyst.

    Parameters
    ----------
    name : str
        Descriptive label for this case.
    olefin_name : str
        Name of the olefin reactant.
    olefin_charge_lb : float
        Initial olefin charge in pounds.
    olefin_mw : float
        Olefin molecular weight (g/mol).
    bisulfite_solution_total_lb : float
        Total aqueous sodium bisulfite solution in pounds.
    bisulfite_wt_fraction : float
        Weight fraction of NaHSO₃ in the aqueous solution (e.g., 0.40 for 40%).
    bisulfite_feed_rate_lb_per_min : float
        Feed rate of the bisulfite solution in lb/min.
    peroxide_name : str
        Name of the organic peroxide catalyst.
    peroxide_total_lb : float
        Total mass of organic peroxide catalyst in pounds.
    peroxide_mw : float
        Peroxide molecular weight (g/mol).
    peroxide_shot_interval_lb : float
        Pounds of bisulfite solution between each peroxide shot.
    product_name : str
        Name of the sulfonation product.
    product_mw : float
        Product molecular weight (g/mol).
    temperature_f : float
        Reaction temperature in °F.
    notes : str or None
        Optional notes or warnings.
    """

    name: str
    olefin_name: str
    olefin_charge_lb: float
    olefin_mw: float
    bisulfite_solution_total_lb: float
    bisulfite_wt_fraction: float
    bisulfite_feed_rate_lb_per_min: float
    peroxide_name: str
    peroxide_total_lb: float
    peroxide_mw: float
    peroxide_shot_interval_lb: float
    product_name: str
    product_mw: float
    temperature_f: float
    notes: Optional[str] = None

    # --- Derived quantities (computed in __post_init__) ---
    _bisulfite_pure_lb: float = field(init=False, repr=False)
    _water_lb: float = field(init=False, repr=False)
    _num_peroxide_shots: int = field(init=False, repr=False)
    _peroxide_per_shot_lb: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._bisulfite_pure_lb = self.bisulfite_solution_total_lb * self.bisulfite_wt_fraction
        self._water_lb = self.bisulfite_solution_total_lb * (1.0 - self.bisulfite_wt_fraction)
        self._num_peroxide_shots = int(
            round(self.bisulfite_solution_total_lb / self.peroxide_shot_interval_lb)
        )
        if self._num_peroxide_shots > 0:
            self._peroxide_per_shot_lb = self.peroxide_total_lb / self._num_peroxide_shots
        else:
            self._peroxide_per_shot_lb = self.peroxide_total_lb

    # ------------------------------------------------------------------
    # Stoichiometric analysis
    # ------------------------------------------------------------------
    def stoichiometry(self) -> Dict[str, float]:
        """Compute molar amounts, limiting reagent, and theoretical yield.

        Returns
        -------
        dict with keys:
            olefin_lbmol, bisulfite_lbmol, peroxide_lbmol,
            molar_ratio_bisulfite_to_olefin, limiting_reagent,
            excess_reagent, excess_lbmol, excess_lb,
            theoretical_product_lbmol, theoretical_product_lb
        """
        n_olefin = lb_to_lbmol(self.olefin_charge_lb, self.olefin_mw)
        n_bisulfite = lb_to_lbmol(self._bisulfite_pure_lb, MW_SODIUM_BISULFITE)
        n_peroxide = lb_to_lbmol(self.peroxide_total_lb, self.peroxide_mw)

        ratio = n_bisulfite / n_olefin if n_olefin > 0.0 else float("inf")

        if n_olefin <= n_bisulfite:
            limiting = self.olefin_name
            excess = "NaHSO3"
            excess_mol = n_bisulfite - n_olefin
            excess_lb = lbmol_to_lb(excess_mol, MW_SODIUM_BISULFITE)
            n_product = n_olefin
        else:
            limiting = "NaHSO3"
            excess = self.olefin_name
            excess_mol = n_olefin - n_bisulfite
            excess_lb = lbmol_to_lb(excess_mol, self.olefin_mw)
            n_product = n_bisulfite

        return {
            "olefin_lbmol": n_olefin,
            "bisulfite_lbmol": n_bisulfite,
            "peroxide_lbmol": n_peroxide,
            "molar_ratio_bisulfite_to_olefin": ratio,
            "limiting_reagent": limiting,
            "excess_reagent": excess,
            "excess_lbmol": excess_mol,
            "excess_lb": excess_lb,
            "theoretical_product_lbmol": n_product,
            "theoretical_product_lb": lbmol_to_lb(n_product, self.product_mw),
        }

    # ------------------------------------------------------------------
    # Feed schedule
    # ------------------------------------------------------------------
    def feed_schedule(self) -> List[FeedScheduleStep]:
        """Build the semi-batch feed schedule.

        Returns a list of FeedScheduleStep objects, one per peroxide-shot
        interval.  Each step covers the addition of ``peroxide_shot_interval_lb``
        pounds of bisulfite solution plus one peroxide shot at the start of
        the interval.
        """
        steps: List[FeedScheduleStep] = []
        cumul_solution = 0.0
        cumul_peroxide = 0.0
        time_cursor = 0.0

        for i in range(1, self._num_peroxide_shots + 1):
            solution_this_step = min(
                self.peroxide_shot_interval_lb,
                self.bisulfite_solution_total_lb - cumul_solution,
            )
            duration = solution_this_step / self.bisulfite_feed_rate_lb_per_min

            peroxide_this_step = self._peroxide_per_shot_lb
            cumul_solution += solution_this_step
            cumul_peroxide += peroxide_this_step

            steps.append(
                FeedScheduleStep(
                    step=i,
                    time_start_min=round(time_cursor, 2),
                    time_end_min=round(time_cursor + duration, 2),
                    bisulfite_solution_added_lb=round(solution_this_step, 2),
                    bisulfite_added_lb=round(solution_this_step * self.bisulfite_wt_fraction, 2),
                    water_added_lb=round(solution_this_step * (1.0 - self.bisulfite_wt_fraction), 2),
                    peroxide_shot_lb=round(peroxide_this_step, 2),
                    cumulative_bisulfite_solution_lb=round(cumul_solution, 2),
                    cumulative_peroxide_lb=round(cumul_peroxide, 2),
                )
            )
            time_cursor += duration

        return steps

    # ------------------------------------------------------------------
    # Mass balance
    # ------------------------------------------------------------------
    def mass_balance(self) -> Dict[str, float]:
        """Overall mass balance for the reactor at full addition.

        Returns
        -------
        dict with keys:
            olefin_charged_lb, bisulfite_solution_charged_lb,
            bisulfite_pure_charged_lb, water_charged_lb,
            peroxide_charged_lb, total_charged_lb,
            theoretical_product_lb, theoretical_product_plus_water_lb,
            total_feed_time_min, total_feed_time_hr,
            num_peroxide_shots, peroxide_per_shot_lb,
            temperature_f
        """
        stoich = self.stoichiometry()
        total_feed_time = (
            self.bisulfite_solution_total_lb / self.bisulfite_feed_rate_lb_per_min
        )
        total_charged = (
            self.olefin_charge_lb
            + self.bisulfite_solution_total_lb
            + self.peroxide_total_lb
        )

        return {
            "olefin_charged_lb": self.olefin_charge_lb,
            "bisulfite_solution_charged_lb": self.bisulfite_solution_total_lb,
            "bisulfite_pure_charged_lb": self._bisulfite_pure_lb,
            "water_charged_lb": self._water_lb,
            "peroxide_charged_lb": self.peroxide_total_lb,
            "total_charged_lb": total_charged,
            "theoretical_product_lb": stoich["theoretical_product_lb"],
            "theoretical_product_plus_water_lb": stoich["theoretical_product_lb"] + self._water_lb,
            "total_feed_time_min": round(total_feed_time, 2),
            "total_feed_time_hr": round(total_feed_time / 60.0, 2),
            "num_peroxide_shots": self._num_peroxide_shots,
            "peroxide_per_shot_lb": round(self._peroxide_per_shot_lb, 2),
            "temperature_f": self.temperature_f,
        }

    # ------------------------------------------------------------------
    # Full run
    # ------------------------------------------------------------------
    def run(self) -> Dict[str, object]:
        """Execute full analysis: stoichiometry, mass balance, and feed schedule.

        Returns
        -------
        dict with keys 'stoichiometry', 'mass_balance', 'feed_schedule'.
        """
        return {
            "stoichiometry": self.stoichiometry(),
            "mass_balance": self.mass_balance(),
            "feed_schedule": self.feed_schedule(),
        }


# ======================================================================
# Default example — 1-octene + 40 % sodium bisulfite + tert-butyl
#                   peroxyacetate
# WARNING: Illustrative parameters only. Replace with measured data!
# ======================================================================

DEFAULT_REACTION_PARAMS = {
    "olefin": {
        "name": "1-octene",
        "mw": MW_1_OCTENE,
        "charge_lb": 7700.0,
    },
    "bisulfite_solution": {
        "total_lb": 22000.0,
        "wt_fraction": 0.40,
        "feed_rate_lb_per_min": 90.0,
    },
    "peroxide": {
        "name": "tert-butyl peroxyacetate",
        "mw": MW_TERT_BUTYL_PEROXYACETATE,
        "total_lb": 435.0,
        "shot_interval_lb": 2200.0,
    },
    "product": {
        "name": "sodium 1-octanesulfonate",
        "mw": MW_SODIUM_1_OCTANESULFONATE,
    },
    "temperature_f": 190.0,
}


def example_sulfonation_case() -> Dict[str, object]:
    """Run the default 1-octene sulfonation example.

    Returns the full results dictionary from ``SulfonationReactionCase.run()``.
    """
    p = DEFAULT_REACTION_PARAMS
    case = SulfonationReactionCase(
        name="1-Octene sulfonation batch",
        olefin_name=p["olefin"]["name"],
        olefin_charge_lb=p["olefin"]["charge_lb"],
        olefin_mw=p["olefin"]["mw"],
        bisulfite_solution_total_lb=p["bisulfite_solution"]["total_lb"],
        bisulfite_wt_fraction=p["bisulfite_solution"]["wt_fraction"],
        bisulfite_feed_rate_lb_per_min=p["bisulfite_solution"]["feed_rate_lb_per_min"],
        peroxide_name=p["peroxide"]["name"],
        peroxide_total_lb=p["peroxide"]["total_lb"],
        peroxide_mw=p["peroxide"]["mw"],
        peroxide_shot_interval_lb=p["peroxide"]["shot_interval_lb"],
        product_name=p["product"]["name"],
        product_mw=p["product"]["mw"],
        temperature_f=p["temperature_f"],
        notes="Illustrative parameters — replace with plant/lab data before decisions.",
    )
    return case.run()


def print_sulfonation_summary(results: Dict[str, object]) -> None:
    """Pretty-print the sulfonation reaction results.

    Parameters
    ----------
    results : dict
        Output from ``SulfonationReactionCase.run()``.
    """
    stoich = results["stoichiometry"]
    mb = results["mass_balance"]
    schedule = results["feed_schedule"]

    print("\n" + "=" * 80)
    print("SULFONATION REACTION SUMMARY")
    print("=" * 80)

    print("\n--- Stoichiometry ---")
    print(f"  Olefin charged:         {stoich['olefin_lbmol']:>10.2f} lb-mol")
    print(f"  NaHSO₃ charged:         {stoich['bisulfite_lbmol']:>10.2f} lb-mol")
    print(f"  Peroxide charged:       {stoich['peroxide_lbmol']:>10.2f} lb-mol")
    print(f"  Molar ratio (NaHSO₃ / olefin): {stoich['molar_ratio_bisulfite_to_olefin']:.3f}")
    print(f"  Limiting reagent:       {stoich['limiting_reagent']}")
    print(f"  Excess reagent:         {stoich['excess_reagent']}")
    print(f"  Excess amount:          {stoich['excess_lbmol']:>10.2f} lb-mol  "
          f"({stoich['excess_lb']:>10.1f} lb)")
    print(f"  Theoretical product:    {stoich['theoretical_product_lbmol']:>10.2f} lb-mol  "
          f"({stoich['theoretical_product_lb']:>10.1f} lb)")

    print("\n--- Mass Balance ---")
    print(f"  Olefin charged:          {mb['olefin_charged_lb']:>10.1f} lb")
    print(f"  Bisulfite solution:      {mb['bisulfite_solution_charged_lb']:>10.1f} lb")
    print(f"    Pure NaHSO₃:           {mb['bisulfite_pure_charged_lb']:>10.1f} lb")
    print(f"    Water:                 {mb['water_charged_lb']:>10.1f} lb")
    print(f"  Peroxide charged:        {mb['peroxide_charged_lb']:>10.1f} lb")
    print(f"  Total material charged:  {mb['total_charged_lb']:>10.1f} lb")
    print(f"  Theoretical product:     {mb['theoretical_product_lb']:>10.1f} lb")
    print(f"  Reaction temperature:    {mb['temperature_f']:>10.1f} °F")

    print("\n--- Feed Schedule ---")
    print(f"  Total feed time:         {mb['total_feed_time_min']:.1f} min  "
          f"({mb['total_feed_time_hr']:.2f} hr)")
    print(f"  Peroxide shots:          {mb['num_peroxide_shots']}")
    print(f"  Peroxide per shot:       {mb['peroxide_per_shot_lb']:.1f} lb")

    print(f"\n  {'Step':>4}  {'Start':>8}  {'End':>8}  {'Bisulfite Soln':>15}  "
          f"{'NaHSO₃':>10}  {'Water':>10}  {'Peroxide':>10}  {'Cumul Soln':>12}  {'Cumul Cat':>10}")
    print(f"  {'':>4}  {'(min)':>8}  {'(min)':>8}  {'(lb)':>15}  "
          f"{'(lb)':>10}  {'(lb)':>10}  {'(lb)':>10}  {'(lb)':>12}  {'(lb)':>10}")
    print("  " + "-" * 98)

    for s in schedule:
        print(f"  {s.step:>4}  {s.time_start_min:>8.1f}  {s.time_end_min:>8.1f}  "
              f"{s.bisulfite_solution_added_lb:>15.1f}  "
              f"{s.bisulfite_added_lb:>10.1f}  {s.water_added_lb:>10.1f}  "
              f"{s.peroxide_shot_lb:>10.1f}  {s.cumulative_bisulfite_solution_lb:>12.1f}  "
              f"{s.cumulative_peroxide_lb:>10.1f}")

    print("\n" + "=" * 80)
    print("NOTE: Parameters are illustrative. Replace with measured data before decisions.")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    results = example_sulfonation_case()
    print_sulfonation_summary(results)
