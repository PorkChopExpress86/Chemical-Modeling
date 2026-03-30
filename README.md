# Chemical Modeling

Python utilities for chemical engineering calculations:

1. **Nitrogen Solubility** (`src/solubility_model.py`) — Estimate dissolved nitrogen in ethylene oxide (EO) and propylene oxide (PO) under storage, and predict pressure rise when heating a closed vessel.
2. **Sulfonation Reaction** (`src/reaction_model.py`) — Model a semi-batch radical sulfonation of a terminal olefin (e.g., 1-octene) with aqueous sodium bisulfite, catalyzed by an organic peroxide.

---

## Nitrogen Solubility Model

### What it models
- Dissolved nitrogen at storage conditions using Henry's law.
- Total moles of nitrogen in a closed tank (gas + dissolved).
- Final nitrogen partial pressure after heating the same closed tank, accounting for lower solubility at higher temperature.

### Key assumptions
- Ideal gas behavior for nitrogen in the headspace.
- Dilute nitrogen in liquid; Henry's law applies.
- Temperature dependence of Henry's constant follows a van't Hoff exponential with an enthalpy of solution, ΔH.
- Closed system; nitrogen moles are conserved while heating.
- Only nitrogen is considered in the headspace (other vapors not modeled here).

### Quick start
```bash
python -m src.solubility_model
```

### Using the API
```python
from src.solubility_model import NitrogenSolubilityCase

case = NitrogenSolubilityCase(
    name="My EO tank",
    volume_liquid_ft3=353.0,
    volume_headspace_ft3=70.6,
    p_initial_psia=87.0,
    t_initial_f=68.0,
    t_final_f=140.0,
    henry_ref_psia_ft3_per_mol=1.16e6,  # replace with data
    henry_tref_f=77.0,
    delta_h_sol_btu_per_mol=10.4,
)
results = case.run()
print(results["final_pressure_psia"], "psia")
```

---

## Sulfonation Reaction Model

### What it models
A semi-batch radical-initiated addition of sodium bisulfite to a terminal olefin:

    R-CH=CH₂  +  NaHSO₃  →  R-CH₂-CH₂-SO₃Na

The default example models 1-octene reacting with 40 % aqueous sodium bisulfite, catalyzed by tert-butyl peroxyacetate at 190 °F.

**Default process parameters:**
| Material | Amount |
|---|---|
| 1-Octene | 7,700 lb |
| 40 % Sodium bisulfite solution | 22,000 lb |
| tert-Butyl peroxyacetate | 435 lb |
| Bisulfite feed rate | 90 lb/min |
| Peroxide shot interval | every 2,200 lb bisulfite solution |
| Reaction temperature | 190 °F |

### Key assumptions
- 1:1 molar stoichiometry (olefin : bisulfite).
- Bisulfite solution is fed at a constant mass rate.
- Peroxide catalyst is added in equal-mass shots at uniform bisulfite-solution intervals.
- Isothermal operation at the specified reaction temperature.
- Kinetic parameters are not included — this is a stoichiometric / mass-balance model.

### Quick start
```bash
python -m src.reaction_model
```

### Using the API
```python
from src.reaction_model import SulfonationReactionCase

case = SulfonationReactionCase(
    name="My sulfonation batch",
    olefin_name="1-octene",
    olefin_charge_lb=7700.0,
    olefin_mw=112.21,
    bisulfite_solution_total_lb=22000.0,
    bisulfite_wt_fraction=0.40,
    bisulfite_feed_rate_lb_per_min=90.0,
    peroxide_name="tert-butyl peroxyacetate",
    peroxide_total_lb=435.0,
    peroxide_mw=132.16,
    peroxide_shot_interval_lb=2200.0,
    product_name="sodium 1-octanesulfonate",
    product_mw=216.28,
    temperature_f=190.0,
)
results = case.run()
# results contains 'stoichiometry', 'mass_balance', and 'feed_schedule'
```

---

## Setup
Create and activate a virtual environment before running:
```bash
# On Debian/Ubuntu, install python3-venv first if needed:
# sudo apt install python3.12-venv

# Create virtual environment
python3 -m venv venv

# Activate (Linux/Mac)
source venv/bin/activate

# Activate (Windows)
venv\Scripts\activate
```

No external packages are required — both models use only the Python standard library.

## Inputs you must supply
All placeholder numbers in the code are **illustrative only**. Replace them with plant or literature values before making process decisions.

## Next steps
- Swap in measured Henry constants for nitrogen in EO/PO at your temperatures.
- Extend to include vapor pressures of EO/PO if total vessel pressure is needed.
- Add reaction kinetics (rate constants, activation energy) to the sulfonation model.
- Include heat-of-reaction calculations for the sulfonation to size cooling systems.
- Add uncertainty bounds if lab measurements have error bars.
