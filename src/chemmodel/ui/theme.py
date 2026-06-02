"""Shared presentation constants and CSS for the Streamlit apps.

Centralises the colour palettes, atom sizes, and page styling that used to be
copy-pasted across ``app.py`` and ``sandbox.py``.
"""

from __future__ import annotations

# ── atom rendering ────────────────────────────────────────────────────────────
ATOM_COLORS = {
    "H": "#e8e8e8", "C": "#505050", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "Cl": "#1FF01F", "Br": "#A62929", "S": "#FFFF30",
    "Pt": "#aaaaaa", "Pd": "#cccccc", "Al": "#f0c040",
}
ATOM_RADII = {
    "H": 7, "C": 13, "N": 12, "O": 12, "F": 11, "Cl": 15, "Br": 16,
    "S": 14, "P": 14, "Na": 16, "Pt": 18, "Pd": 17, "Al": 14,
}

# group colours for reaction systems (A = blue, B = green, catalyst = gold)
GROUP_COLORS = ["#5b9bd5", "#70ad47", "#ffc000"]

# ── molecular-sandbox input lanes ─────────────────────────────────────────────
LANES = ["Reactant A", "Reactant B", "Catalyst/Additive", "Secondary Molecule"]
LANE_COLORS = {
    "Reactant A": "#5b9bd5", "Reactant B": "#70ad47",
    "Catalyst/Additive": "#ffc000", "Secondary Molecule": "#c0504d",
}

# ── page CSS ──────────────────────────────────────────────────────────────────
SIMULATION_CSS = """
<style>
.stApp { background: #0f0f1a; }
.block-container { padding-top: 1.2rem; }
.engine-badge {
    display:inline-block; padding:3px 10px; border-radius:12px;
    font-size:0.78rem; font-weight:600;
    background:#1a3a5c; color:#7ec8e3; border:1px solid #7ec8e3;
}
</style>
"""

SANDBOX_CSS = """
<style>
.stApp { background: #0f0f1a; }
.block-container { padding-top: 1.1rem; max-width: 1400px; }
.lane-badge { display:inline-block; padding:2px 10px; border-radius:10px;
  font-size:0.75rem; font-weight:700; color:#0f0f1a; }
.engine-badge { display:inline-block; padding:3px 10px; border-radius:12px;
  font-size:0.78rem; font-weight:600; background:#1a3a5c; color:#7ec8e3;
  border:1px solid #7ec8e3; }
</style>
"""
