# Chemical Modeling Copilot Instructions

## Project Purpose
This is a specialized chemical engineering tool that models nitrogen solubility in ethylene oxide (EO) and propylene oxide (PO) under storage conditions, and predicts pressure rise when heating closed vessels. It uses Henry's law with temperature-dependent corrections.

## Core Architecture

### Single-Module Design
- All modeling logic lives in [src/solubility_model.py](../src/solubility_model.py)
- No external dependencies beyond Python standard library
- Runs as a script or imported as an API

### Key Physics
The model combines three thermodynamic relationships:
1. **Henry's law**: Gas solubility proportional to partial pressure (P = H·x)
2. **van't Hoff equation**: Temperature-dependent Henry constant with enthalpy of solution (ΔH_sol)
3. **Ideal gas law**: For nitrogen in headspace
4. **Mass conservation**: Closed system where nitrogen redistributes between dissolved and gas phases when heated

## Critical Conventions

### Units Must Stay Explicit
- **Never use ambiguous units**. All functions expect US customary:
  - Pressure: psia (pounds per square inch absolute)
  - Temperature: °R (degrees Rankine, convert from Fahrenheit by adding 459.67)
  - Volume: ft³ (cubic feet)
  - Henry constant: psia·ft³/mol
- User-facing `NitrogenSolubilityCase` accepts psia and Fahrenheit for convenience, converting °F to °R internally

### Henry Constant Convention
- Uses `Hcp` formulation: `P = Hcp * x` (partial pressure = Henry constant × mole fraction)
- Positive `delta_h_sol` → solubility decreases with temperature (typical for gases)
- `delta_h_sol` in BTU/mol, with gas constant R = 10.7316 psia·ft³/(mol·°R)
- See [solubility_model.py#L38-L47](../src/solubility_model.py#L38-L47) for implementation

### Placeholder Data Warning
- **`DEFAULT_PARAMS` are illustrative only**, not real chemical data
- Code must warn users to replace with measured values before decisions
- Any new features should maintain this cautionary approach

## Development Workflows

### Virtual Environment Setup
Always use a virtual environment for Python development:
```bash
# On Debian/Ubuntu systems, install python3-venv first if needed:
# sudo apt install python3.12-venv

# Create virtual environment (once)
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate  # Linux/Mac
venv\Scripts\activate     # Windows
```

### Running the Model
```bash
# Run examples with placeholder data (venv must be activated)
python -m src.solubility_model
```

No external packages need to be installed - the model uses only Python standard library.

### Using as API
Import `NitrogenSolubilityCase`, instantiate with tank geometry + thermodynamic parameters, call `.run()` to get results dictionary. See [README.md](../README.md#using-the-api) for example.

### Adding New Solvents
1. Get measured Henry parameters (H_ref, T_ref, ΔH_sol) from literature or lab
2. Add to `DEFAULT_PARAMS` dictionary
3. Create example case in `example_case()` function
4. Document assumptions in module docstring

## Common Pitfalls
- **Don't mix units**: Use conversion constants (`BAR_TO_PSIA`, `PSIA_TO_BAR`)
- **Temperature in Rankine**: All thermodynamic calculations use absolute temperature (°R = °F + 459.67)
- **Closed system assumption**: Model conserves nitrogen moles; doesn't apply to vented or actively purged scenarios
- **Dilute regime**: Henry's law breaks down for high solubility; keep nitrogen partial pressures reasonable

## Extension Points
- Add vapor pressure of EO/PO to compute total vessel pressure (currently only nitrogen partial pressure)
- Include uncertainty propagation if Henry parameters have error bars
- Support multi-component headspace (e.g., nitrogen + oxygen)
