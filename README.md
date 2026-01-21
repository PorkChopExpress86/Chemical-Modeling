# Chemical Modeling – Nitrogen Solubility in EO/PO

Python utilities to estimate how much nitrogen is dissolved in ethylene oxide (EO) and propylene oxide (PO) while stored under nitrogen, and what pressure rise to expect when heating a closed vessel.

## What this models
- Dissolved nitrogen at storage conditions using Henry's law.
- Total moles of nitrogen in a closed tank (gas + dissolved).
- Final nitrogen partial pressure after heating the same closed tank, accounting for lower solubility at higher temperature.

## Key assumptions
- Ideal gas behavior for nitrogen in the headspace.
- Dilute nitrogen in liquid; Henry's law applies.
- Temperature dependence of Henry's constant follows a van't Hoff exponential with an enthalpy of solution, ΔH.
- Closed system; nitrogen moles are conserved while heating.
- Only nitrogen is considered in the headspace (other vapors not modeled here).

## Inputs you must supply
Accurate Henry parameters for nitrogen in EO and PO:
- `H_ref` (Pa·m³/mol) at a reference temperature `T_ref` (°C).
- `ΔH_sol` (J/mol) for nitrogen dissolution.
- Tank geometry: liquid volume and headspace volume (m³).
- Pressures in bar (gauge or absolute — be consistent) and temperatures in °C.

The placeholder numbers in the code are illustrative. Replace them with plant or literature values before making decisions.

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

## Quick start
```bash
python -m src.solubility_model
```
This runs EO and PO example cases with placeholder parameters and prints results.

## Using the API
Create a `NitrogenSolubilityCase` and call `run()`:
```python
from src.solubility_model import NitrogenSolubilityCase

case = NitrogenSolubilityCase(
    name="My EO tank",
    volume_liquid_m3=10.0,
    volume_headspace_m3=2.0,
    p_initial_bar=6.0,
    t_initial_c=20.0,
    t_final_c=60.0,
    henry_ref_pa_m3_per_mol=8.0e7,  # replace with data
    henry_tref_c=25.0,
    delta_h_sol_j_per_mol=11000.0,
)
results = case.run()
print(results["final_pressure_bar"], "bar")
```

## Next steps
- Swap in measured Henry constants for nitrogen in EO/PO at your temperatures.
- Extend to include vapor pressures of EO/PO if total vessel pressure (not just N₂ partial pressure) is needed.
- Add uncertainty bounds if lab measurements have error bars.
