"""chemmodel — chemical & reactor modelling toolkit.

Layered architecture:

- ``chemmodel.constants``  — physical constants / unit conversions (single source).
- ``chemmodel.solubility`` — equilibrium N₂ solubility in EO/PO (stdlib-only core).
- ``chemmodel.engines``    — quantum-chemistry engines as a Strategy + Registry
                             (Psi4 / PySCF / CP2K / Demo).
- ``chemmodel.chemistry``  — molecular structures, PubChem lookup, preset libraries.
- ``chemmodel.reactions``  — potential-energy-surface scans (engine-agnostic).
- ``chemmodel.reactors``   — closed-form ideal-reactor kinetics (Batch/CSTR/PFR).
- ``chemmodel.md``         — ab-initio Langevin molecular-dynamics worker.
- ``chemmodel.ui``         — thin Streamlit presentation layers (import the above).

Domain modules never import Streamlit; all UI lives under ``chemmodel.ui``.
"""

from __future__ import annotations

__version__ = "0.1.0"
